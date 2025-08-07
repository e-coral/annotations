# for annotating spreadsheets of Chr Start End Signal

import os
from pathlib import Path
from src.annotate import annotate
import pandas


indir = Path('/home/emmon/OneDrive/Julie/CapCruncher_results/')
outdir = Path('/home/emmon/OneDrive/Julie/CapCruncher_results/annotated_bedgraphs/')


def main():
    """
    For all relevant bedgraphs in the subdirectories, read in the data, annotate it, and return annotated spreadsheets
    to relevant locations
    :return:
    """
    for dir in indir.iterdir():
        if not str(dir).endswith('old') and not str(dir).endswith('raw_bedgraphs') and not str(dir).endswith('troubleshooting'):
            dirname = dir.name
            wholepath = Path(os.path.join(indir, dirname, f"comparisons/bedgraphs"))
            for filename in wholepath.glob(r'Late*'):
                outname = f"annotated_{dirname}_{filename.stem}"
                print(outname)
                outfile = os.path.join(outdir, outname)
                annotate_bedgraph_as_excel(filename, outfile)
            for filename in wholepath.glob(r'Early.*'):
                outname = f"annotated_{dirname}_{filename.stem}"
                print(outname)
                outfile = os.path.join(outdir, outname)
                annotate_bedgraph_as_excel(filename, outfile)
            for filename in wholepath.glob(r'UT*'):
                outname = f"annotated_{dirname}_{filename.stem}"
                print(outname)
                outfile = os.path.join(outdir, outname)
                annotate_bedgraph_as_excel(filename, outfile)


def annotate_bedgraph_as_excel(infile, outfile):
    """
    For all relevant bedgraphs in the subdirectories, read in the data, annotate it, and return annotated spreadsheets
    to relevant locations
    :return: annotated version of the input file
    """
    data = pandas.read_csv(infile, sep='\t')
    data.columns = ["chrom", "start", "end", "signal"]
    annotate.annotate_dataframe_output_excel(data, outfile)


if __name__ == '__main__':
    main()