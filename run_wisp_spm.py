import os
import sys
import time
import argparse
import numpy as np
import textwrap
import pickle
from wisp.run import run_wisp
from wisp.contexts import ContextManager
from wisp.structure import Atom, Molecule
from wisp.traj import multi_threading_to_collect_data_from_frames, collect_data_from_frames
from wisp.paths import multi_threading_find_paths, GetPaths
from wisp.utils import GetCovarianceMatrix
from wisp.viz import Visualize
from wisp.io import output_dir_info


### SET DIRECTORIES ###
WORKING_DIR = os.path.dirname(os.path.abspath(__file__))
FILE_DIR = os.path.join(WORKING_DIR, "input/")
WRITING_DIR = os.path.join(WORKING_DIR,"tmp/")

if __name__ == "__main__":
    program_start_time = time.time()

    # get the commandline parameters
    parameters = UserInput()

    # compute the correlation matrix
    if parameters["load_wisp_saved_matrix"] == "TRUE":
        correlation_matrix_object = pickle.load(
            open(parameters["wisp_saved_matrix_filename"], "rb")
        )  # load the matrix instead of generating
    else:
        correlation_matrix_object = GetCovarianceMatrix(
            parameters
        )  # so generate the matrix instead of loading it
    correlation_matrix = correlation_matrix_object.correlations

    # always save a copy of the correlation matrix, regardless of how it was loaded/generated
    pickle.dump(
        correlation_matrix_object,
        open(
            parameters["output_directory"]
            + "functionalized_matrix_with_contact_map_applied.pickle",
            "wb",
        ),
    )

    # now get the source and sink locations from the parameters
    sources = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        parameters["source_residues"]
    )
    sinks = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        parameters["sink_residues"]
    )

    # compute the paths
    paths = GetPaths(
        correlation_matrix,
        sources,
        sinks,
        parameters,
        correlation_matrix_object.average_pdb.residue_identifiers_in_order,
    )

    # create the visualization
    vis = Visualize(parameters, correlation_matrix_object, paths)

    # provide the user information about the generated files
    output_dir_info(parameters)

    log(
        "\n# Program execution time: "
        + str(time.time() - program_start_time)
        + " seconds",
        parameters["logfile"],
    )

    parameters["logfile"].close()