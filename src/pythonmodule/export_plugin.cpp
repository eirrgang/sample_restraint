/*! \file
 * \brief Provide Python bindings and helper functions for setting up restraint potentials.
 *
 * There is currently a lot of boilerplate here that will be generalized and removed in a future version.
 * In the mean time, follow the example for EnsembleRestraint to create the proper helper functions
 * and instantiate the necessary templates.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#include "export_plugin.h"

#include <cassert>

#include <memory>

#include "gmxapi/exceptions.h"
#include "gmxapi/md.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/gmxapi.h"

#include "ensemblepotential.h"

// Make a convenient alias to save some typing...
namespace py = pybind11;

template<>
std::shared_ptr<gmxapi::MDModule> PyRestraint<plugin::RestraintModule<plugin::EnsembleRestraint>>::getModule()
{
    return shared_from_this();
}
//////////////////////////////////////////////////////////////////////////////////////////
// New restraints mimicking EnsembleRestraint should specialize getModule() here as above.
//////////////////////////////////////////////////////////////////////////////////////////

class EnsembleRestraintBuilder
{
    public:
        explicit EnsembleRestraintBuilder(py::object element)
        {
            name_ = py::cast<std::string>(element.attr("name"));
            assert(!name_.empty());

            // It looks like we need some boilerplate exceptions for plugins so we have something to
            // raise if the element is invalid.
            assert(py::hasattr(element,
                               "params"));

            // Params attribute should be a Python list
            py::dict parameter_dict = element.attr("params");
            // \todo Check for the presence of these dictionary keys to avoid hard-to-diagnose error.

            // Get positional parameters.
            py::list sites = parameter_dict["sites"];
            for (auto&& site : sites)
            {
                siteIndices_.emplace_back(py::cast<int>(site));
            }

            auto nbins = py::cast<size_t>(parameter_dict["nbins"]);
            auto binWidth = py::cast<double>(parameter_dict["binWidth"]);
            auto minDist = py::cast<double>(parameter_dict["min_dist"]);
            auto maxDist = pybind11::cast<double>(parameter_dict["max_dist"]);
            auto experimental = pybind11::cast<std::vector<double>>(parameter_dict["experimental"]);
            auto nSamples = pybind11::cast<unsigned int>(parameter_dict["nsamples"]);
            auto samplePeriod = pybind11::cast<double>(parameter_dict["sample_period"]);
            auto nWindows = pybind11::cast<unsigned int>(parameter_dict["nwindows"]);
            auto k = pybind11::cast<double>(parameter_dict["k"]);
            auto sigma = pybind11::cast<double>(parameter_dict["sigma"]);

            auto params = plugin::makeEnsembleParams(nbins,
                                                     binWidth,
                                                     minDist,
                                                     maxDist,
                                                     experimental,
                                                     nSamples,
                                                     samplePeriod,
                                                     nWindows,
                                                     k,
                                                     sigma);
            params_ = std::move(*params);

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
            if (!subscriber_)
            {
                return;
            }
            else
            {
                if (!py::hasattr(subscriber_, "potential")) throw gmxapi::ProtocolError("Invalid subscriber");
            }

            // Temporarily subvert things to get quick-and-dirty solution for testing.
            // Need to capture Python communicator and pybind syntax in closure so EnsembleResources
            // can just call with matrix arguments.

            // This can be replaced with a subscription and delayed until launch, if necessary.
            if (!py::hasattr(context_, "ensemble_update"))
            {
                throw gmxapi::ProtocolError("context does not have 'ensemble_update'.");
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
            auto resources = std::make_shared<plugin::EnsembleResources>(std::move(functor));

            auto potential = PyRestraint<plugin::RestraintModule<plugin::EnsembleRestraint>>::create(name_,
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

        py::object subscriber_;
        py::object context_;
        std::vector<int> siteIndices_;

        plugin::ensemble_input_param_type params_;

        std::string name_;
};

/*!
 * \brief Factory function to create a new builder for use during Session launch.
 *
 * \param element WorkElement provided through Context
 * \return ownership of new builder object
 */
std::unique_ptr<EnsembleRestraintBuilder> createEnsembleBuilder(const py::object& element)
{
    using std::make_unique;
    auto builder = make_unique<EnsembleRestraintBuilder>(element);
    return builder;
}


////////////////////////////////////////////////////////////////////////////////////////////
// New potentials modeled after EnsembleRestraint should define a Builder class and define a
// factory function here, following the previous two examples. The factory function should be
// exposed to Python following the examples near the end of the PYBIND11_MODULE block.
////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
// The PYBIND11_MODULE block uses the pybind11 framework (ref https://github.com/pybind/pybind11 )
// to generate Python bindings to the C++ code elsewhere in this repository. A copy of the pybind11
// source code is included with this repository. Use syntax from the examples below when exposing
// a new potential, along with its builder and parameters structure. In future releases, there will
// be less code to include elsewhere, but more syntax in the block below to define and export the
// interface to a plugin. pybind11 is not required to write a GROMACS extension module or for
// compatibility with the ``gmx`` module provided with gmxapi. It is sufficient to implement the
// various protocols, C API and Python function names, but we do not provide example code
// for other Python bindings frameworks.
//////////////////////////////////////////////////////////////////////////////////////////////////

// The first argument is the name of the module when importing to Python. This should be the same as the name specified
// as the OUTPUT_NAME for the shared object library in the CMakeLists.txt file. The second argument, 'm', can be anything
// but it might as well be short since we use it to refer to aspects of the module we are defining.
PYBIND11_MODULE(myplugin, m) {
    m.doc() = "sample plugin"; // This will be the text of the module's docstring.

    // Matrix utility class (temporary). Borrowed from http://pybind11.readthedocs.io/en/master/advanced/pycpp/numpy.html#arrays
    py::class_<plugin::Matrix<double>, std::shared_ptr<plugin::Matrix<double>>>(m,
                                                                                "Matrix",
                                                                                py::buffer_protocol())
        .def_buffer([](plugin::Matrix<double>& matrix) -> py::buffer_info {
            return py::buffer_info(
                matrix.data(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                2,                                      /* Number of dimensions */
                {matrix.rows(), matrix.cols()},                 /* Buffer dimensions */
                {sizeof(double) * matrix.cols(),             /* Strides (in bytes) for each index */
                 sizeof(double)}
            );
        });

    //////////////////////////////////////////////////////////////////////////
    // Begin EnsembleRestraint
    //
    // Define Builder to be returned from ensemble_restraint Python function defined further down.
    pybind11::class_<EnsembleRestraintBuilder> ensembleBuilder(m,
                                                               "EnsembleBuilder");
    ensembleBuilder.def("add_subscriber",
                        &EnsembleRestraintBuilder::addSubscriber);
    ensembleBuilder.def("build",
                        &EnsembleRestraintBuilder::build);

    // Get more concise name for the template instantiation...
    using PyEnsemble = PyRestraint<plugin::RestraintModule<plugin::EnsembleRestraint>>;

    // Export a Python class for our parameters struct
    py::class_<plugin::EnsembleRestraint::input_param_type> ensembleParams(m, "EnsembleRestraintParams");
    m.def("make_ensemble_params",
          &plugin::makeEnsembleParams);

    // API object to build.
    py::class_<PyEnsemble, std::shared_ptr<PyEnsemble>> ensemble(m, "EnsembleRestraint");
    // EnsembleRestraint can only be created via builder for now.
    ensemble.def("bind",
                 &PyEnsemble::bind,
                 "Implement binding protocol");
    /*
     * To implement gmxapi_workspec_1_0, the module needs a function that a Context can import that
     * produces a builder that translates workspec elements for session launching. The object returned
     * by our function needs to have an add_subscriber(other_builder) method and a build(graph) method.
     * The build() method returns None or a launcher. A launcher has a signature like launch(rank) and
     * returns None or a runner.
     */

    // Generate the name operation that will be used to specify elements of Work in gmxapi workflows.
    // WorkElements will then have namespace: "myplugin" and operation: "ensemble_restraint"
    m.def("ensemble_restraint",
          [](const py::object element) { return createEnsembleBuilder(element); });
    //
    // End EnsembleRestraint
    ///////////////////////////////////////////////////////////////////////////




}
