# for just annotating the boundaries, split the positions for each event into two, and then shove them back together

import os
from pathlib import Path
from annotate import annotate
import pandas


outdir = (Path(__file__).parent / '06082024').resolve()
refs_dir = (Path(__file__).parent / 'annotate/ref_files').resolve()


def format_columns(df, orig_col_names):
    """
    Drop irrelevant columns from the dataframe, and rename to original names
    :param df: dataframes of input values and annotations
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

    # remove nans and set values to strings to prevent trailing .0s
    df.fillna('', inplace=True)
    df = df.astype(str)

    return df


def annotate_results():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    # for all the input files
    for filename in refs_dir.glob(r'EC2*'):
        # generate the relevant output file name
        a_outname = f"posA_basic_annotated_{filename.stem}.xlsx"
        b_outname = f"posB_basic_annotated_{filename.stem}.xlsx"

        # read in the input file
        dysgu_df = pandas.read_csv(filename, sep="\t", header=0)

        if not dysgu_df.empty:
            # rename the columns to use the A boundary for analysis, and add an end column for compatibility
            dysgu_df.columns = ['seqid', 'start', 'chrB', 'posB']
            dysgu_df['end'] = dysgu_df['start'] + 1

            # annotate distances
            dysgu_data = annotate.calculate_distances_to_telomeres(dysgu_df)
            # print(dysgu_data.head())
            res = annotate.annotate_overlaps(dysgu_data)
            print(res.head())

            # reformat the final df to match the input
            a_final = format_columns(res, ['chr', 'start', 'chrB', 'posB', 'end'])

            # create the output file
            with (pandas.ExcelWriter(os.path.join(outdir, a_outname)) as writer):
                # write output file
                a_final.to_excel(writer, index=False)

            # rename to use B data
            dysgu_df.columns = ['chrA', 'posA', 'seqid', 'end', 'start']
            dysgu_df['start'] = dysgu_df['end'] - 1

            # annotate
            dysgu_data = annotate.calculate_distances_to_telomeres(dysgu_df)
            res = annotate.annotate_overlaps(dysgu_data)

            # reformat the final df to match the input
            b_final = format_columns(res, ['chrA', 'posA', 'chr', 'end', 'start'])

            # create the output file
            with (pandas.ExcelWriter(os.path.join(outdir, b_outname)) as writer):
                # write output file
                b_final.to_excel(writer, index=False)


if __name__ == '__main__':
    annotate_results()
