/*! \file
 * \brief Provide Python bindings and helper functions for setting up restraint potentials.
 *
 * There is currently a lot of boilerplate here that will be generalized and removed in a future version.
 * In the mean time, follow the example for MDStringRestraint to create the proper helper functions
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

#include "mdstring_potential.h"

// Make a convenient alias to save some typing...
namespace py = pybind11;

namespace plugin
{

//////////////////////////////////////////////////////////////////////////////////////////
// New restraints mimicking MDStringRestraint should specialize getModule() here.
//////////////////////////////////////////////////////////////////////////////////////////
template<>
std::shared_ptr<gmxapi::MDModule> PyRestraint<RestraintModule<Restraint<MDStringPotential>>>::getModule()
{
    return shared_from_this();
}


using MDStringRestraintBuilder = RestraintBuilder<MDStringPotential>;


/*!
 * \brief Factory function to create a new builder for use during Session launch.
 *
 * \param element WorkElement provided through Context
 * \return ownership of new builder object
 */
std::unique_ptr<MDStringRestraintBuilder> createMDStringBuilder(const py::object element)
{
    using std::make_unique;
    auto builder = make_unique<MDStringRestraintBuilder>(element);
    builder->add_input("nbins", &mdstring_data_t::nBins)
            .add_input("binWidth", &mdstring_data_t::binWidth)
            .add_input("min_dist", &mdstring_data_t::minDist)
            .add_input("max_dist", &mdstring_data_t::maxDist)
            .add_input("experimental", &mdstring_data_t::experimental)
            .add_input("nsamples", &mdstring_data_t::nSamples)
            .add_input("sample_period", &mdstring_data_t::samplePeriod)
            .add_input("nwindows", &mdstring_data_t::nWindows)
            .add_input("k", &mdstring_data_t::sigma)
            .add_input("sigma", &mdstring_data_t::sigma);
    return builder;
}


////////////////////////////////////////////////////////////////////////////////////////////
// New potentials modeled after MDStringRestraint should define a Builder class and define a
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
PYBIND11_MODULE(mdstring, m) {
    m.doc() = "String method for molecular dynamics"; // This will be the text of the module's docstring.

    // Matrix utility class (temporary). Borrowed from http://pybind11.readthedocs.io/en/master/advanced/pycpp/numpy.html#arrays
    py::class_<Matrix2D<double>, std::shared_ptr<Matrix2D<double>>>(m,
                                                                                "Matrix2D",
                                                                                py::buffer_protocol())
        .def_buffer([](Matrix2D<double>& matrix) -> py::buffer_info {
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
    // Begin MDStringRestraint
    //
    // Define Builder to be returned from mdstring_restraint Python function defined further down.
    pybind11::class_<MDStringRestraintBuilder> mdstringBuilder(m,
                                                               "MDStringBuilder");
    mdstringBuilder.def("add_subscriber",
                        &MDStringRestraintBuilder::addSubscriber);
    mdstringBuilder.def("build",
                        &MDStringRestraintBuilder::build);

    // Define a more concise name for the template instantiation...
    using PyMDString = PyRestraint<RestraintModule<Restraint<MDStringPotential>>>;

    // Export a Python class for our parameters struct
    py::class_<Restraint<MDStringPotential>::input_param_type> mdstringParams(m, "MDStringRestraintParams");
    m.def("make_mdstring_params",
          &makeMDStringParams);

    // API object to build.
    py::class_<PyMDString, std::shared_ptr<PyMDString>> mdstring(m, "MDStringRestraint");
    // MDStringRestraint can only be created via builder for now.
    mdstring.def("bind",
                 &PyMDString::bind,
                 "Implement binding protocol");
    /*
     * To implement gmxapi_workspec_1_0, the module needs a function that a Context can import that
     * produces a builder that translates workspec elements for session launching. The object returned
     * by our function needs to have an add_subscriber(other_builder) method and a build(graph) method.
     * The build() method returns None or a launcher. A launcher has a signature like launch(rank) and
     * returns None or a runner.
     */

    // Generate the name operation that will be used to specify elements of Work in gmxapi workflows.
    // WorkElements will then have namespace: "myplugin" and operation: "mdstring_restraint"
    m.def("mdstring_restraint",
          [](const py::object element) { return createMDStringBuilder(element); });
    //
    // End MDStringRestraint
    ///////////////////////////////////////////////////////////////////////////
}

} // end namespace plugin
