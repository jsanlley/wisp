import os
import time
import argparse
import numpy as np
import pickle
import multiprocessing
from wisp.run import run_wisp
from wisp.contexts import ContextManager
from wisp.structure import Atom, Molecule
from wisp.traj import multi_threading_to_collect_data_from_frames, collect_data_from_frames
from wisp.paths import multi_threading_find_paths, GetPaths
from wisp.utils import GetCovarianceMatrix
from wisp.viz import Visualize
from wisp.io import output_dir_info, UserInput, log

### SET DIRECTORIES ###
WORKING_DIR = os.path.join(os.getcwd(),"wisp-spm-output")

if __name__ == "__main__":
    print(f'WORKING_DIR: {WORKING_DIR}')

    program_start_time = time.time()

    # get the commandline parameters
    argparser = argparse.ArgumentParser(description="WISP-SPM:Shortest Path Method using the Weighted Implementation of Suboptimal Paths")
    argparser.add_argument(
        "-prmtop",
        type=str,
        required=False,
        help="Topology file in AMBER format (prmtop)",)
    
    argparser.add_argument(
        "-nc",
        type=str,
        required=False,
        help="Trajectory file in AMBER format (nc)",
    )
    argparser.add_argument(
        "-pdb",
        type=str,
        required=True,
        help="The filename of the multi-frame PDB to analyze. Individual frames should be separated by 'END' or 'ENDMDL' lines.",
    )
    argparser.add_argument(
        "-distmat",
        type=str,
        required=True,
        help="The filename of the pre-computed distance matrix to use.",
    )
    argparser.add_argument(
        "-corrmat",
        type=str,
        required=True,
        help="The filename of the pre-computed correlation matrix to use.",
    )
    argparser.add_argument(
        "-config",
        type=str,
        required=False,
        help="The filename of the configuration file to use.",
    )

    args = argparser.parse_args()

    # Input settings from the command line
    context = ContextManager()
    
    # Default settings
    context.pdb_path=args.pdb
    context.contact_map_path=args.distmat
    context.functionalized_matrix_path=args.corrmat
    context.output_dir=WORKING_DIR    
    context.n_paths = 1
    context.contact_map_distance_limit = 4.5
    context.n_cores = 4
    context.frame_chunks = 96
    context.node_definition  = "RESIDUE COM" # Double check if CA in script

    if context.wisp_saved_matrix_path:
        correlation_matrix_object = pickle.load(
            open(context.wisp_saved_matrix_path, "rb")
        )
    else:
        correlation_matrix_object = GetCovarianceMatrix(context)

    correlation_matrix = correlation_matrix_object.correlations

    pickle.dump(
        correlation_matrix_object,
        open(
            os.path.join(
                context.output_dir,
                "functionalized_matrix_with_contact_map_applied.pickle",
            ),
            "wb",
        ),
    )
       # now get the source and sink locations from the parameters
    source_residues = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        context["source_residues"]
    )
    sink_residues = correlation_matrix_object.convert_list_of_residue_keys_to_residue_indices(
        context["sink_residues"]
    )

    # compute the paths
    paths = GetPaths(
        correlation_matrix_object.correlations,
        source_residues,
        sink_residues,
        context,
        correlation_matrix_object.average_pdb.residue_identifiers_in_order,
    )