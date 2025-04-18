#!/usr/bin/env python3

import argparse
import yaml
import sys
import os

# local import
from src import preprocess


def  update_params(new_params, template_params):
    for proc, params in template_params['params'].items():
        if proc not in new_params['params']:
            continue
        for p, value in params.items():
            if p in new_params['params'][proc]:
                template_params['params'][proc][p] = new_params['params'][proc][p]
    return template_params


def main():
    parser = argparse.ArgumentParser()

    group1 = parser.add_argument_group('\nMandatory arguments')
    group1.add_argument(
        "-s", dest="spectra", type=str, required=True,
        help="Path to spectra files in mzML format"
    )
    group1.add_argument('-o', dest="output", type=str, required=True,
        help="Path to output file. Output format is a mzML file"
    )

    group2 = parser.add_argument_group("\nPreprocessing parameters")
    group2.add_argument('-p', '--paramfile', dest="param_file", type=str, required=False,
        help="Path yo YAML with parameters for preprocessing"
    )
    group2.add_argument('-r', '--repl_to_sample', dest="repl_to_sample", type=str, required=False,
        help="Table in csv mapping replicate spectra names/IDs to sample name for replicate merging."
             "Each row contains \"spectra_name,sample_name\", eg. \"abc_1,abc\""
             "If not given, it is inferred from the spectrum IDs (abc_1, abc_2, abc_3 -> abc)"
    )

    group3 = parser.add_argument_group("\nPlotting arguments")
    group3.add_argument("-q", "--plotnpeaks", dest="make_plot", required=False, action='store_true',
        help="Whether to plot number of detected peaks per spectrum. It is stored in the sme directory as the output."
             "It is useful to assess assess if the signal-to-noise threshold is reasonable."
    )

    args = parser.parse_args()

    default_file = os.path.dirname(sys.argv[0]) + '/preprocessing_params.yml'
    with open(default_file, 'r') as p:
        prep_params = yaml.safe_load(p)
    if args.param_file:
        with open(args.param_file, 'r') as p:
            user_params = yaml.safe_load(p)
        prep_params = update_params(user_params, prep_params)

    preprocess.main(args.spectra, args.output, prep_params, args.repl_to_sample, args.make_plot)

if __name__ == "__main__":
    main()