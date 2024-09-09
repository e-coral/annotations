import os
from pathlib import Path
from src.annotate import annotate


outdir = (Path(__file__).parent.parent / 'example_output_files').resolve()
indir = (Path(__file__).parent.parent / 'example_input_files').resolve()


# def format_columns(df, orig_col_names):
#     """
#     Drop irrelevant columns from the dataframe, and rename to original names
#     :param df: dataframe of input values and annotations
#     :param orig_col_names: list of the original column names
#     :return: formatted dataframe of values and annotations
#     """
#     # add the names of the new columns to the names of the original columns
#     for col in ['tChr', 'tStart', 'tEnd', 'Centromere', 'region-telomere_distance', 'repeats', 'fragile_sites',
#                 'genes', 'gene_lengths', 'g-t_distance']:
#         orig_col_names.append(col)
#
#     # rename the columns accordingly
#     df.columns = orig_col_names
#
#     # drop the telomere and centromere columns that aren't needed in the final output
#     df = df.drop(columns=['tChr', 'tStart', 'tEnd', 'Centromere'])
#
#     # set values to strings to prevent trailing .0s
#     df = df.astype(str)
#
#     return df
#

def annotate_captureC_spreadsheet():
    """
    For all the captureC spreadsheets, read in the data, annotate it, and return the annotated file
    :return: annotated version of the input file
    """
    # for all the input files
    for filename in indir.glob(r'*CaC*.xlsx'):
        # generate the relevant output file names and paths
        outname = f"annotated_{filename.stem}"
        infile = os.path.join(indir, filename)
        outfile = os.path.join(outdir, outname)

        # annotate the input and write excel and csv output files, with exploded gene annotations (one row per gene)
        annotate.annotate_standard_excel_input_file(infile=infile, outfile=outfile, colname="chrom",
                                                    explode=False, make_csv=False)


if __name__ == '__main__':
    annotate_captureC_spreadsheet()
