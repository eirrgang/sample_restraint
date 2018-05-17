#!/usr/bin/env python
"""Test BRMC restraint"""

import sys
# Of course, the location of the Python plugin module is user-specific and could be
# passed by PYTHONPATH instead of programatically here.
sys.path.append('/Users/jmh/dev/sample_restraint/cmake-build-release/src/pythonmodule')

import gmx

import logging
logging.getLogger().setLevel(logging.DEBUG)
# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handler
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s: %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logging.getLogger().addHandler(ch)
logger = logging.getLogger()

import myplugin

logger.info("myplugin is {}".format(myplugin.__file__))

# Work with single ensemble member for now.
tpr_list = ['/Users/jmh/dev/syx-2017/test_syx.tpr']

print("tpr list: {}".format(tpr_list))


# dt = 0.002
# First restraint applied between atoms 422 and 2955
# Second restraint applied between atom 1287 and 3059
# Restraint site coordinates relative to atom 2028

params1 = {
    'sites': [422, 2955, 2028],
    'A': 20,
    'tau': 0.5,
    'tolerance': 0.05,
    'target': 5.0,
    'nSamples': 3,
    'parameter_filename': "params52210.log"
    # 'samplePeriod': 0.1
}

params2 = {
    'sites': [1287, 3059, 2028],
    'A': 20,
    'tau': 0.5,
    'tolerance': 0.05,
    'target': 4.0,
    'nSamples': 3,
    'parameter_filename': "params105216.log"
    # 'samplePeriod': 0.1
}
#
potential1 = gmx.workflow.WorkElement(
    namespace="myplugin",
    operation="brmc_restraint",
    depends=[],
    params=params1
)
potential1.name = "brmc_restraint_1"

potential2 = gmx.workflow.WorkElement(
    namespace="myplugin",
    operation="brmc_restraint",
    depends=[],
    params=params2
)
potential2.name = "brmc_restraint_2"

md = gmx.workflow.from_tpr(tpr_list)
md.add_dependency(potential1)
md.add_dependency(potential2)

context = gmx.context.ParallelArrayContext(md)

with context as session:
    session.run()

