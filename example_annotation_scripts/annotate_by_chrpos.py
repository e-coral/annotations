"""
annotate regions that are in the format chr:start-end, held within multiple sheets of an excel file
"""
import os
from pathlib import Path
from src.annotate import annotate


outdir = (Path(__file__).parent.parent / 'example_output_files').resolve()
indir = (Path(__file__).parent.parent / 'example_input_files').resolve()


def main():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated output file
    """
    # for all relevant input files
    for filename in indir.glob(r'example_multisheet_chr-start-end_input.xlsx'):
        # generate the relevant output file name
        outname = f"annotated_{filename.stem}.xlsx"

        # generate the paths for inputs and outputs
        infile = os.path.join(indir, filename)
        outfile = os.path.join(outdir, outname)

        # annotate the input data and output as a file
        annotate.annotate_gene_location_data(infile, outfile, "Gene location")


if __name__ == '__main__':
    main()
