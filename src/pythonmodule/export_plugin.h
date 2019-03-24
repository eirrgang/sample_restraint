#ifndef GMXAPI_SAMPLE_RESTRAINT_EXPORT_PLUGIN_H
#define GMXAPI_SAMPLE_RESTRAINT_EXPORT_PLUGIN_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "gmxapi/exceptions.h"
#include "gmxapi/md.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/gmxapi.h"

#include "restraint.h"
#include "sessionresources.h"

namespace py = pybind11;

namespace plugin
{
////////////////////////////////
// Begin PyRestraint static code
/*!
 * \brief Templated wrapper to use in Python bindings.
 *
 * Boilerplate
 *
 * Mix-in from below. Adds a bind behavior, a getModule() method to get a gmxapi::MDModule adapter,
 * and a create() method that assures a single shared_ptr record for an object that may sometimes
 * be referred to by a raw pointer and/or have shared_from_this called.
 * \tparam T class implementing gmx::IRestraintPotential
 *
 */
template<class T>
class PyRestraint : public T, public std::enable_shared_from_this<PyRestraint<T>>
{
    public:
        void bind(py::object object);

        using T::name;

        /*!
         * \brief
         *
         * T must either derive from gmxapi::MDModule or provide a template specialization for
         * PyRestraint<T>::getModule(). If T derives from gmxapi::MDModule, we can keep a weak pointer
         * to ourself and generate a shared_ptr on request, but std::enable_shared_from_this already
         * does that, so we use it when we can.
         * \return
         */
        std::shared_ptr<gmxapi::MDModule> getModule();

        /*!
         * \brief Factory function to get a managed pointer to a new restraint.
         *
         * \tparam ArgsT
         * \param args
         * \return
         */
        template<typename ... ArgsT>
        static std::shared_ptr<PyRestraint<T>> create(ArgsT... args)
        {
            auto newRestraint = std::make_shared<PyRestraint<T>>(args...);
            return newRestraint;
        }

        template<typename ... ArgsT>
        explicit PyRestraint(ArgsT... args) :
            T{args...}
        {}

};

/*!
 * \brief Implement the gmxapi binding protocol for restraints.
 *
 * All restraints will use this same code automatically.
 *
 * \tparam T restraint class exported below.
 * \param object Python Capsule object to allow binding with a simple C API.
 */
template<class T>
void PyRestraint<T>::bind(py::object object)
{
    PyObject* capsule = object.ptr();
    if (PyCapsule_IsValid(capsule,
                          gmxapi::MDHolder::api_name))
    {
        auto holder = static_cast<gmxapi::MDHolder*>(PyCapsule_GetPointer(capsule,
                                                                          gmxapi::MDHolder::api_name));
        auto workSpec = holder->getSpec();
        std::cout << this->name() << " received " << holder->name();
        std::cout << " containing spec of size ";
        std::cout << workSpec->getModules().size();
        std::cout << std::endl;

        auto module = getModule();
        workSpec->addModule(module);
    }
    else
    {
        throw gmxapi::ProtocolError("bind method requires a python capsule as input");
    }
}
// end PyRestraint static code
//////////////////////////////


/*!
 * \brief Interact with the restraint framework and gmxapi when launching a simulation.
 *
 * This should be generalized and removed from here. Unfortunately, some things need to be
 * standardized first. If a potential follows the example of EnsembleRestraint or HarmonicRestraint,
 * the template specializations below can be mimicked to give GROMACS access to the potential.
 *
 * \tparam T class implementing the gmxapi::MDModule interface.
 * \return shared ownership of a T object via the gmxapi::MDModule interface.
 */
// If T is derived from gmxapi::MDModule, create a default-constructed std::shared_ptr<T>
// \todo Need a better default that can call a shared_from_this()
template<class T>
std::shared_ptr<gmxapi::MDModule> PyRestraint<T>::getModule()
{
    auto module = std::make_shared<typename std::enable_if<std::is_base_of<gmxapi::MDModule, T>::value, T>::type>();
    return module;
}

template<class PotentialT>
class RestraintBuilder
{
    public:
        explicit RestraintBuilder(py::object element)
        {
            name_ = py::cast<std::string>(element.attr("name"));
            assert(!name_.empty());

            // It looks like we need some boilerplate exceptions for plugins so we have something to
            // raise if the element is invalid.
            assert(py::hasattr(element,
                               "params"));

            // Params attribute should be a Python dict
            parameter_dict_ = element.attr("params");
            // \todo Check for the presence of these dictionary keys to avoid hard-to-diagnose error.

            // Get positional parameters.
            py::list sites = parameter_dict_["sites"];
            for (auto&& site : sites)
            {
                siteIndices_.emplace_back(py::cast<int>(site));
            }


            // Note that if we want to grab a reference to the Context or its communicator, we can get it
            // here through element.workspec._context. We need a more general API solution, but this code is
            // in the Python bindings code, so we know we are in a Python Context.
            assert(py::hasattr(element,
                               "workspec"));
            auto workspec = element.attr("workspec");
            assert(py::hasattr(workspec,
                               "_context"));
            context_ = workspec.attr("_context");
        }

        /*!
         * \brief Add node(s) to graph for the work element.
         *
         * \param graph networkx.DiGraph object still evolving in gmx.context.
         *
         * \todo This may not follow the latest graph building protocol as described.
         */
        void build(py::object graph)
        {
            // For each registered input, call the stored function object to set
            // the C++ data from the provided Python data.
            for(const auto& setter : setters_)
            {
                setter(&params_, parameter_dict_);
            }


            // Temporarily subvert things to get quick-and-dirty solution for testing.
            // Need to capture Python communicator and pybind syntax in closure so Resources
            // can just call with matrix arguments.

            // This can be replaced with a subscription and delayed until launch, if necessary.
            if (!py::hasattr(context_, "ensemble_update"))
            {
                throw gmxapi::ProtocolError("context does not have 'mdstring_update'.");
            }
            // make a local copy of the Python object so we can capture it in the lambda
            auto update = context_.attr("ensemble_update");
            // Make a callable with standardizeable signature.
            const std::string name{name_};
            auto functor = [update, name](const plugin::Matrix<double>& send,
                                          plugin::Matrix<double>* receive) {
                update(send,
                       receive,
                       py::str(name));
            };

            // To use a reduce function on the Python side, we need to provide it with a Python buffer-like object,
            // so we will create one here. Note: it looks like the SharedData element will be useful after all.
            auto resources = std::make_shared<plugin::Resources>(std::move(functor));

            auto potential = PyRestraint<plugin::RestraintModule<plugin::Restraint<PotentialT>>>::create(name_,
                                                                                                                        siteIndices_,
                                                                                                                        params_,
                                                                                                                        resources);

            auto subscriber = subscriber_;
            py::list potentialList = subscriber.attr("potential");
            potentialList.append(potential);

        };

        /*!
         * \brief Accept subscription of an MD task.
         *
         * \param subscriber Python object with a 'potential' attribute that is a Python list.
         *
         * During build, an object is added to the subscriber's self.potential, which is then bound with
         * system.add_potential(potential) during the subscriber's launch()
         */
        void addSubscriber(py::object subscriber)
        {
            assert(py::hasattr(subscriber,
                               "potential"));
            subscriber_ = subscriber;
        };

        /*!
         * \brief Register an input name and storage location.
         *
         * Example:
         *
         *      builder.add_setter("nbins", &input_param_type::nBins);
         */
        template<typename T>
        RestraintBuilder& add_input(const std::string& name,
                                    T PotentialT::input_param_type::* data_ptr)
        {
            auto setter =
                [=](typename PotentialT::input_param_type* p, const py::dict& d)->void
            {
                p->*data_ptr = py::cast<T>(d[name.c_str()]);
            };
            this->setters_.emplace_back(setter);
            return *this;
        };

        py::dict parameter_dict_;
        py::object subscriber_;
        py::object context_;
        std::vector<int> siteIndices_;

        typename PotentialT::input_param_type params_;

        std::string name_;

        std::vector<std::function<void(typename PotentialT::input_param_type*, const py::dict&)>> setters_;
};

} // end namespace plugin

#endif //GMXAPI_SAMPLE_RESTRAINT_EXPORT_PLUGIN_H
