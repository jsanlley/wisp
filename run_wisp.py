import os
import sys
import time
import argparse

from wisp.run import run_wisp
from wisp.contexts import ContextManager

### SET DIRECTORIES ###
WORKING_DIR = os.path.join(os.getcwd(), "wisp-spm-output")

if __name__ == "__main__":
    print(f'WORKING_DIR: {WORKING_DIR}')

    program_start_time = time.time()

    # get the commandline parameters
    argparser = argparse.ArgumentParser(description="WISP-SPM:Shortest Path Method using the Weighted Implementation of Suboptimal Paths")
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
        "-source",
        type=str,
        nargs="+",
        required=False,
        help="Source residues for path generation, e.g. '-source X_SER_1 X_LEU_4'. Overrides source_residues from -config if both are given.",
    )
    argparser.add_argument(
        "-sink",
        type=str,
        nargs="+",
        required=False,
        help="Sink residues for path generation, e.g. '-sink X_ARG_37'. Overrides sink_residues from -config if both are given.",
    )
    argparser.add_argument(
        "-config",
        type=str,
        required=False,
        help="The filename of the YAML configuration file to use.",
    )
    argparser.add_argument(
        "-output",
        type=str,
        required=False,
        help="The output directory to write results to. Defaults to 'wisp-spm-output' in the current directory.",
    )

    args = argparser.parse_args()

    # Input settings, optionally seeded from a YAML config file
    context = ContextManager(yaml_paths=args.config)

    context.pdb_path = args.pdb
    context.contact_map_path = args.distmat
    context.functionalized_matrix_path = args.corrmat
    context.output_dir = args.output or WORKING_DIR

    if args.source:
        context.source_residues = args.source
    if args.sink:
        context.sink_residues = args.sink

    context.n_paths = 1  # SPM: single shortest path only
    context.contact_map_distance_limit = 4.5
    context.n_cores = 4
    context.frame_chunks = 96
    context.node_definition = "RESIDUE_COM"

    if not context.source_residues or not context.sink_residues:
        print(
            "\nNo source and/or sink residues specified. Provide them via -source/-sink "
            "on the command line, or set source_residues/sink_residues in the -config YAML file."
        )
        sys.exit(1)

    os.makedirs(context.output_dir, exist_ok=True)

    paths = run_wisp(context)

    print(
        "\n# Program execution time: "
        + str(time.time() - program_start_time)
        + " seconds"
    )
