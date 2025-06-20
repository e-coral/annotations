# for region boundary files in the format:
# chr A, pos A, chr B, pos B
# annotate the positions at both boundaries

import os
from pathlib import Path
from src.annotate import annotate
import pandas


outdir = Path('/home/emmon/OneDrive/Kate/GI2/results/')
indir = Path('/home/emmon/OneDrive/Kate/GI2/from_kate/')


def add_distances(df, orig_cols):
    """
    dataframe including seqid and start for the position whose distance to telomere needs to be calculated
    @param df: pandas dataframe of positions
    @param orig_cols: list of original column names
    :return: df with tel distance added
    """
    # include 'end' as an original column for the formatting function
    orig_cols.append('end')
    # print(orig_cols)

    # make sure the start position is the right type
    df['start'] = df['start'].astype(float)

    # add an end column, required by some methods
    df['end'] = df['start'] + 1

    # annotate positions
    dist_data = annotate.calculate_distances_to_centromeres_and_telomeres(df)
    # print(dist_data.head())

    # reformat df to remove temp data, including the end column
    pos_distances = annotate.format_just_distance_output_columns(dist_data, orig_cols)
    pos_distances = pos_distances.drop(columns=['end'])

    return pos_distances


def annotate_results():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    # for all the input files
    for filename in indir.glob(r'*'):
        # generate the relevant output file name
        outname = f"just_distances-{filename.stem}.xlsx"
        # print(filename)

        # read in the input file
        # dysgu_df = pandas.read_csv(filename, sep="\t", header=0)
        infile = pandas.read_excel(filename, sheet_name=None)
        # for each sheet in the input file (probably only one for each in this case)
        for sname, sheet in infile.items():
            # if there are data in the file
            if not sheet.empty:
                # get the position info from the first columns
                s = sheet[sheet.columns[0:4]]
                # print(f"s = {s.dtypes}")

                # rename the columns to the required values to use the A boundary for analysis,
                s.columns = ['seqid', 'start', 'chrB', 'posB']

                # calculate and add the telomere distances
                dists = add_distances(s, orig_cols=['chr', 'start', 'chrB', 'posB'])
                # print(f"dists = {dists.dtypes}")

                # rename the columns to switch to B values
                dists.columns = ['chrA', 'posA', 'seqid', 'start', 'tel_distanceA']
                # print(dists.head())

                # calculate and add the telomere distances
                b_dists = add_distances(dists, orig_cols=['chrA', 'posA', 'chr', 'start', 'tel_distanceA'])
                b_dists.columns = ['chrA', 'posA', 'chrB', 'posB', 'tel_distanceA', 'tel_distanceB']
                # print(b_dists.head())

                # b_dists['posA'] = b_dists['posA'].astype(float)
                # b_dists['posB'] = b_dists['posB'].astype(float)

                # final = pandas.merge(sheet, b_dists, on=['chrA', 'posA', 'chrB', 'posB'])
                # print(final.head())
                # create the output file
                with (pandas.ExcelWriter(os.path.join(outdir, outname)) as writer):
                    # write output file
                    # final.to_excel(writer, index=False)
                    b_dists.to_excel(writer, index=False)


if __name__ == '__main__':
    annotate_results()
