import os
from pathlib import Path
from annotate import annotate
import pandas
import glob


outdir = (Path(__file__).parent / 'other_capC').resolve()
refs_dir = (Path(__file__).parent / 'annotate/ref_files').resolve()


def format_columns(df, orig_col_names):
    """
    Drop irrelevant columns from the dataframe, and rename to original names
    :param df: dataframe of input values and annotations
    :param orig_col_names: list of the original column names
    :return: formatted dataframe of values and annotations
    """
    # add the names of the new columns to the names of the original columns
    for col in ['tChr', 'tStart', 'tEnd', 'Centromere', 'region-telomere_distance', 'repeats', 'fragile_sites',
                'genes', 'gene_lengths', 'g-t_distance', 'gene_positions']:
        orig_col_names.append(col)

    # rename the columns accordingly
    df.columns = orig_col_names

    # drop the telomere and centromere columns that aren't needed in the final output
    df = df.drop(columns=['tChr', 'tStart', 'tEnd', 'Centromere'])

    # set values to strings to prevent trailing .0s
    df = df.astype(str)

    return df


def annotate_results():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    # for all the input files
    for filename in refs_dir.glob(r'*100624.xlsx'):
        # generate the relevant output file name
        outname = f"basic_annotated_{filename.stem}.xlsx"

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

                    # rename the chr column, because some of the methods may be based on 'seqid'
                    s.rename(columns={"chrom": "seqid"}, inplace=True)

                    # calculate the distance between the region and the nearest telomere
                    res = annotate.calculate_distances_to_telomeres(s)

                    # annotate the genes, fragile sites, repeats and gene sizes
                    regions_df = annotate.annotate_overlaps(res)

                    # reformat the final df to match the input
                    final = format_columns(regions_df, orig_columns)

                    # write output file
                    final.to_excel(writer, sheet_name=sname, index=False)


if __name__ == '__main__':
    annotate_results()
