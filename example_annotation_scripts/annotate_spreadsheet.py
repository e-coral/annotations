import os
from pathlib import Path
from src.annotate import annotate


outdir = (Path(__file__).parent.parent / 'example_output_files').resolve()
indir = (Path(__file__).parent.parent / 'example_input_files').resolve()


def annotate_spreadsheet():
    """
    read in the data, annotate it, and return the annotated file
    :return: annotated version of the input file
    """
    # for all the input files
    for filename in indir.glob(r'*CaC.xlsx'):
        # generate the relevant output file names and paths
        outname = f"annotated_{filename.stem}"
        infile = os.path.join(indir, filename)
        outfile = os.path.join(outdir, outname)

        # annotate the input and write excel and csv output files, with exploded gene annotations (one row per gene)
        annotate.annotate_standard_excel_input_file(infile=infile, outfile=outfile, colname="chrom",
                                                    explode=False, make_csv=False)


if __name__ == '__main__':
    annotate_spreadsheet()
