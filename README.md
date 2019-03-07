MDString restraint plugin
==========================

Repository Contents
-------------------

This repository uses CMake to build and install a Python C++ extension
package.

-   `CMakeLists.txt`, `cmake/FindGROMACS.cmake`, and
    `src/CMakeLists.txt` provide necessary CMake infrastructure. You
    should not need to edit these.
-   `src/cpp` contains a header and `cpp` file for each restraint
    potential built with this module. When adding new potentials, you
    will update `CMakeLists.txt` to create build targets. Use the
    existing potentials as examples.
-   `src/pythonmodule/` contains `CMakeLists.txt`, `export_plugin.h`,
    and `export_plugin.cpp`. When you have written a new potential, you
    can add it to `CMakeLists.txt` and `export_plugin.cpp`. 
    `MDStringPotential` uses additional facilities provided by gmxapi.
-   `src/pybind11` is just a copy of the Python bindings framework from
    the Pybind project (ref <https://github.com/pybind/pybind11> ). It
    is used to wrap the C++ restraint code and give it a Python
    interface.
-   `tests/` contains C++ and Python tests for the provided code. Update
    `CMakeLists.txt` to add your own, based on these examples. C++ unit
    tests use [googletest](https://github.com/google/googletest). Python
    tests use the [pytest](https://docs.pytest.org/en/latest/). Refer to
    those respective projects for more about how they make test-writing
    easier.
-   `Dockerfile` is a recipe to build a Docker image from the root of
    the repository.

The basics
----------

### Build and install

To download, build, and install, you may need to first install `wget`,
`git`, and/or `cmake`.

The plugin requires libgmxapi to build. See
[gromacs-gmxapi](https://github.com/kassonlab/gromacs-gmxapi) :

    # install GROMACS. Instead of `master`, you can choose a specific release or the `devel` branch.
    wget https://github.com/kassonlab/gromacs-gmxapi/archive/master.zip
    unzip master.zip
    cd gromacs-gmxapi-master
    mkdir build
    cd mkdir build
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/gromacs -DGMX_THREAD_MPI=ON
    make install # use -j10 to build in parallel with 10 cores (or however many you have)
    cd ../..

We use CMake to configure and build a C++ library and a Python module
for interacting with it. After installing the modified GROMACS (see
above), either source the GMXRC file provided with the GROMACS
installation or provide the install location to CMake with the
`gmxapi_DIR` environment variable.

As with [gmxapi](https://github.com/kassonlab/gromacs-gmxapi), we
recommend installing and using this code in a Python virtual
environment. (See the documentation for your `gmxapi` distribution or
<http://gmxapi.readthedocs.io/en/latest/install.html> ) Accordingly, if
you choose to install the plugin rather than just to use it out of
its build directory, consider whether you want to have to set your
`PYTHONPATH` environment variable or where you can install it that
Python will find it. You can explicitly set the installation location by
setting `-DGMXPLUGIN_INSTALL_PATH=/path/to/install/directory` or you can
let CMake determine an appropriate location automatically for your
Python interpreter. If you have administrative privileges (such as when
running on a desktop computer) or if you are using a Python virtual
environment (recommended), you don\'t need to specify anything
additional. If you are an unprivileged user (such as on a shared
machine) and are not in a Python virtual environment, set
-DGMXPLUGIN\_USER\_INSTALL=ON to install into the \"user\" Python
packages directory in your home directory. (Equivalent to the `--user`
option to `pip`)

If you have multiple Python installations or just want to be
unambiguous, provide CMake with the Python interpreter you wish to use
(the same as you are using for `gmxapi`) with
`-DPYTHON_EXECUTABLE=/path/to/python`. For instance, if you have both
Python 3.x and Python 2.7, but you plan to use Python 2.7, use
`` -DPYTHON_EXECUTABLE=`which python2` `` or
`` -DPYTHON_EXECUTABLE=`which python` `` (if `python` points to the
Python 2 interpreter). :

    # build sample restraint
    git clone https://github.com/kassonlab/sample_restraint.git
    # optionally, check out the development branch
    # pushd sample_restraint ; git checkout devel ; popd
    # perform an out-of-source build
    mkdir build
    cd build
    # Get the GROMACS environment settings
    source $HOME/gromacs/bin/GMXRC
    # Configure the build environment with CMake
    cmake ../sample_restraint
    # or
    # cmake ../sample_restraint -DGMXPLUGIN_INSTALL_PATH=/path/to/install/directory
    # or
    # cmake ../sample_restraint -DGMXPLUGIN_USER_INSTALL=ON -DPYTHON_EXECUTABLE=`which python`
    make
    # run C++ tests
    make test
    # optionally, install
    make install

If you choose not to install the plugin module, you can tell Python
where to find it by setting your PYTHONPATH environment variable. For
instance, while still in the build directory:

    export PYTHONPATH=`pwd`/src/pythonmodule

The Python module `gmx` is required for testing. See
[gmxapi](https://github.com/kassonlab/gmxapi)

### What\'s going on

This sample project builds several C++ libraries with names such as
`harmonicpotential`. The actual filename will be something like
`libharmonicpotential.so` or `harmonicpotential.dll` or something
depending on your operating system. These libraries are used to build a
Python module named `myplugin`.

When setting up a workflow, a Python script provides gmxapi with
parameters and a factory function for a plugin restraint potential. This
Python interface is defined in `src/pythonmodule/export_plugin.cpp`.
When a Session is launched, an C++ object that performs restraint force
calculations is created and given to the GROMACS library. During each MD
step, part of the MD force evaluation includes a call to the
calculations performed by the restraint. For the pair restraints
demonstrated here, GROMACS provides relative coordinates of two atomic
sites to the calculation code in the plugin. If multiple restrained
pairs are needed, multiple restraints are attached to the simulation.
Coordination across an ensemble of simulations is possible using
resources provided by the Session.

Fundamentally, a new restraint potential is implemented by creating a
class that provides a `calculate()` method and using wrappers to give it
interfaces to GROMACS and to Python. C++ wrappers allow the basic class
implementing the potential to be presented to the GROMACS library in a
way that can be used to evaluate forces during a simulation. Other C++
template code wraps the potential in a portable way so that it can be
passed to GROMACS through a Python interface and to receive parameters
from the Python interpreter. Pybind11 syntax in `export_plugin.cpp`
provides the code to actually expose the plugin as a class in a Python
module that is compatible with the `gmx` package provided in the
`gmxapi` project.

Python tests
------------

For the Python-level testing, you will need `pytest` and `gmxapi`. We
recommend setting up a Python virtual environment as described at
[<https://github.com/kassonlab/gmxapi>](https://github.com/kassonlab/gmxapi)

You will also need a functioning MPI installation and the `mpi4py`
package.

Python tests can be run from the root directory of the repository after
building. Assuming you built in a subdirecory of the repository named
`build` (as above):

    PYTHONPATH=build/src/pythonmodule/ python -m pytest tests

This command causes the directory named `tests` to be explored for
Python files with names like `test_*.py` or `*_test.py`. Matching files
will be imported and any functions with similarly obvious names will be
run and errors reported. In particular, `assert` statements will be
evaluated to perform individual tests. See also
<https://docs.pytest.org/en/latest/goodpractices.html#test-discovery>

The tests assume that the package is already installed or is available
on the default Python path (such as by setting the `PYTHONPATH`
environment variable). If you just run `pytest` with no arguments, it
will discover and try to run tests from elsewhere in the repository that
were not intended, and they will fail.

To run the full set of tests for the ensemble workflow features, first
make sure that you have an MPI-capable environment and `mpi4py`
installed. Refer to <http://mpi4py.readthedocs.io/en/stable/> and
<https://github.com/kassonlab/gmxapi> for more information.

The ensemble tests assume that 2 ranks are available. After installing
the plugin, run (for example):

    mpiexec -n 2 python -m mpi4py -m pytest

If you do not have MPI set up for your system, you could build a docker
image using the Dockerfile in this repository.

    docker build -t samplerestraint . Dockerfile
    docker run --cpus 2 --rm -ti samplerestraint bash -c \
        "cd /home/jovyan/sample_restraint/tests && 
        mpiexec -n 2 python -m mpi4py -m pytest"

To test with a pre-built image from our docker hub repository, do

    docker run --cpus 2 --rm -ti gmxapi/sample_restraint bash -c \
            "cd /home/jovyan/sample_restraint/tests && 
            mpiexec -n 2 python -m mpi4py -m pytest"
