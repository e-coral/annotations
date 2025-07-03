# for annotating a csv file containing information about supplementary alignments

import pandas
import os
from src.annotate import annotate


mydir = os.path.dirname('/media/emmon/5918f239-6127-4c44-b6b3-7c8cd9ae538c/Sara/pacbio/to_annotate/')


def additional_filtering(ann_file, filename):
    """
    perform additional filtering based on annotations, mapq, etc.
    :param ann_file: annotated file
    :param filename: original input filename
    :return: filtered output file
    """
    # read the annotated file in as a df
    ann_df = pandas.read_csv(ann_file, sep=",")

    # filter to longer, more unique, higher apq reads not in repeat/centromere/telomere regions
    r_df = ann_df[ann_df['percentage_unique'] > 50]
    r_df = r_df[r_df['unique'] == True]
    r_df = r_df[r_df['centromeric'] == False]
    r_df = r_df[r_df['telomeric'] == False]
    r_df = r_df[r_df['SA_length'] >= 200]
    r_df = r_df[r_df['mapq'] >= 50]

    outname = os.path.join(mydir, f"filtered_annotated_{filename}")
    r_df.to_csv(outname, index=False)

    qnames_list_outname = os.path.join(mydir, f"qnames_filtered_for_{filename}.txt")
    qnames = set(r_df['qname'].tolist())
    with open(qnames_list_outname, 'w') as outfile:
        outfile.write("\n".join(qnames))


def annotate_standard_csv():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    # for whichever input csv files you want to annotate
    for filename in os.listdir(mydir):
        # print(filename)
        # generate the relevant output file names and paths
        outname = f"annotated_{filename}"
        infile = os.path.join(mydir, filename)
        outfile = os.path.join(mydir, outname)

        # annotate the input csv and output it to csv
        annotate.annotate_standard_csv_input_file(infile, outfile, colname="chr", explode=False)

        additional_filtering(outfile, filename)


if __name__ == '__main__':
    annotate_standard_csv()
