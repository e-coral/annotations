import os
from pathlib import Path
from src.annotate import annotate
import pandas


outdir = (Path(__file__).parent.parent / 'example_output_files').resolve()
indir = (Path(__file__).parent.parent / 'example_input_files').resolve()


def annotate_standard_csv():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    # for whichever input csv files you want to annotate
    for filename in indir.glob('*standard.csv'):
        # generate the relevant output file names and paths
        outname = f"annotated_{filename.stem}"
        infile = os.path.join(indir, filename)
        outfile = os.path.join(outdir, outname)

        # annotate the input csv and output it to csv
        annotate.annotate_standard_csv_input_file(infile, outfile, colname="seqnames", explode=False)



if __name__ == '__main__':
    annotate_standard_csv()
