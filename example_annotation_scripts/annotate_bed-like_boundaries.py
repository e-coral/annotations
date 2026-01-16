# for csv-format region boundary files, with the columns:
# chr A, pos A, chr B, pos B
# annotate the positions at both boundaries
# and output each as an excel file

import os
from pathlib import Path
from src.annotate import annotate
import pandas


outdir = (Path(__file__).parent.parent / 'example_output_files').resolve()
indir = (Path(__file__).parent.parent / 'example_input_files').resolve()


def annotate_results():
    """
    read in the file, add the annotations, and output annotated excel files
    :return: annotated file
    """
    # for all the input files
    for filename in indir.glob(r'example_bed-like*'):
        # generate the relevant output file name
        a_outname = f"posA_annotated_{filename.stem}.xlsx"
        b_outname = f"posB_annotated_{filename.stem}.xlsx"
        # print(filename)

        # read in the input file
        in_df = pandas.read_csv(filename, sep="\t", header=0)

        # if there are data in the file
        if not in_df.empty:
            # rename the columns to use the A boundary for analysis, and add an 'end' column for compatibility
            in_df.columns = ['seqid', 'start', 'chrB', 'posB']
            in_df['end'] = in_df['start'] + 1

            # annotate A positions
            dysgu_data = annotate.calculate_distances_to_centromeres_and_telomeres(in_df)
            res = annotate.annotate_overlaps(dysgu_data)

            # reformat the final df to match the input
            a_final = annotate.format_output_columns(res, ['chr', 'start', 'chrB', 'posB', 'end'])

            # create the output file
            with (pandas.ExcelWriter(os.path.join(outdir, a_outname)) as writer):
                # write output file
                a_final.to_excel(writer, index=False)

            # rename columns to use chr, start and end for the B positions
            in_df.columns = ['chrA', 'posA', 'seqid', 'end', 'start']
            # create start positions for compatibility with modules
            in_df['start'] = in_df['end'] - 1

            # annotate
            dysgu_data = annotate.calculate_distances_to_centromeres_and_telomeres(in_df)
            res = annotate.annotate_overlaps(dysgu_data)

            # reformat the final df to match the input
            b_final = annotate.format_output_columns(res, ['chrA', 'posA', 'chr', 'end', 'start'])

            # create the output file
            with (pandas.ExcelWriter(os.path.join(outdir, b_outname)) as writer):
                # write output file
                b_final.to_excel(writer, index=False)


if __name__ == '__main__':
    annotate_results()
