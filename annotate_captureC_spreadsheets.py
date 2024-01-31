import os
from pathlib import Path
from annotate import annotate
import pandas


outdir = (Path(__file__).parent / 'capC_output').resolve()
refs_dir = (Path(__file__).parent / 'annotate/ref_files').resolve()
print(refs_dir)


def format_columns(df, orig_col_names):
    """
    Drop irrelevant columns from the dataframe, and rename to original names
    :param df: dataframe of input values and annotations
    :param orig_col_names: list of the original column names
    :return: formatted dataframe of values and annotations
    """

    return df


def annotate_captureC_spreadsheet():
    for filename in refs_dir.glob(r'CaC*.xlsx'):
        outname = f"annotated_{filename.stem}.xlsx"
        # print(outname)

        with pandas.ExcelWriter(os.path.join(outdir, outname)) as writer:
            infile = pandas.read_excel(filename, sheet_name=None)
            print(f"infile length = {len(infile)}")
            for sname, s in infile.items():
                if not s.empty:
                    orig_columns = [col.strip() for col in s.columns]
                    s.rename(columns={"chrom": "seqid"}, inplace=True)
                    s = annotate.calculate_distances_to_telomeres(s)

                    regions_df = annotate.annotate_overlaps(s)

                    print(regions_df.head())
                    exit()

                    final = format_columns(s, orig_columns)
                    final.to_excel(writer, sheet_name=sname, index=False)



                    # trim away the additional columns

            #
            #         # add gene lengths, telomere distances and gene locations
            #         s_df = annotate.calculate_gene_lengths(s_df)
            #
            #         # print(s_df.columns)
            #         s_df = annotate.calculate_distances_to_telomeres(s_df)
            #         # print(s_df.columns)
            #         s_df = annotate.add_gene_location(s_df)
            #
            #         # remove the irrelevant columns
            #         s_df = s_df[["Orig_gene", "Alt_gene", "Gene", "Chr", "start", "end", "gene_length", "g-t_distance", "full_gene_loc"]]
            #
            #         # rename columns
            #         s_df.columns = ["Orig_Gene", "Alt_Gene", "Gene", "Chr", "Start", "End", "Gene_length", "Distance_to_telomere", "Gene_loc"]


if __name__ == '__main__':
    annotate_captureC_spreadsheet()