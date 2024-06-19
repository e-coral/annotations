import os
from pathlib import Path
from annotate import annotate
import pandas


outdir = (Path(__file__).parent / 'other_capC').resolve()
refs_dir = (Path(__file__).parent / 'annotate/ref_files').resolve()


def format_columns(df, orig_col_names):
    """
    Drop irrelevant columns from the dataframe, and rename to original names
    :param df: dataframe of input values and annotations
    :param orig_col_names: list of the original column names
    :return: formatted dataframe of values and annotations
    """
    df.rename(columns={"seqid": 'chrom'}, inplace=True)
    df.fillna('', inplace=True)
    # print(df.columns)

    # drop the telomere and centromere columns that aren't needed in the final output
    df = df.drop(columns=['Chr', 'Start', 'End', 'Centromere'])

    # set values to strings to prevent trailing .0s
    df = df.astype(str)

    return df


def annotate_spreadsheet():
    """
    For all the spreadsheets, read in the data, annotate it, and return the annotated file
    :return: annotated version of the input file
    """
    # for all the input files
    for filename in refs_dir.glob(r'*100624.xlsx'):
        # generate the relevant output file name
        outname = f"annotated_{filename.stem}.xlsx"

        # create the output file
        with pandas.ExcelWriter(os.path.join(outdir, outname)) as writer:
            # read in the input file
            infile = pandas.read_excel(filename, sheet_name=None)
            # for each sheet in the input file (probably only one for each in this case)
            for sname, s in infile.items():
                # if the sheet is not empty
                if not s.empty:
                    # store the original column names to be able to add them back in later
                    orig_columns = [col.strip() for col in s.columns]

                    # rename the chrom column, because some of the methods may be based on 'seqid'
                    s.rename(columns={"chrom": "seqid"}, inplace=True)

                    # calculate the distance between the region and the nearest telomere
                    s = annotate.calculate_distances_to_telomeres(s)

                    # annotate the fragile sites
                    fs_anns = annotate.annotate_just_fragile_sites(s)

                    # add the gene lengths, distances, etc. to the genes in the input file
                    genes_df = annotate.annotate_existing_genes(fs_anns, "Gene")

                    genes_df.rename(columns={"orig_gene": "Gene"}, inplace=True)

                    fs_anns.to_csv(os.path.join(outdir, "fs_anns.csv"), index=False)
                    genes_df.to_csv(os.path.join(outdir, "genes_anns.csv"), index=False)

                    # merge the gene anns df with the orig + tel & fs df
                    full_df = pandas.merge(fs_anns, genes_df, how='left', on='Gene')
                    # print(fs_anns.columns)
                    full_df.to_csv(os.path.join(outdir, "merged.csv"), index=False)

                    # reformat the final df to match the input
                    final = format_columns(full_df, orig_columns)

                    # create output files
                    final.to_excel(writer, sheet_name=sname, index=False)
                    # final.to_csv(os.path.join(outdir, f"{sname}-exploded.csv"), index=False)


if __name__ == '__main__':
    annotate_spreadsheet()
