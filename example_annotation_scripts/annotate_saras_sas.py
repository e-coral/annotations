# for annotating a list of genes in a single column in an excel file

import pandas
import glob
from pathlib import Path
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
    print(len(r_df))
    r_df = r_df[r_df['unique'] == True]
    print(len(r_df))
    r_df = r_df[r_df['centromeric'] == False]
    r_df = r_df[r_df['telomeric'] == False]
    print(len(r_df))
    r_df = r_df[r_df['SA_length'] >= 200]
    r_df = r_df[r_df['mapq'] >= 50]
    print(len(r_df))
    # r_df = r_df[r_df['repeats'] == ]

    # for rname in r_df.groupby('qname'):
    #     print(rname)

    outname = os.path.join(mydir, f"filtered_annotated_{filename}")
    r_df.to_csv(outname, index=False)
    print(len(set(r_df['qname'])))

def annotate_standard_csv():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    # for whichever input csv files you want to annotate
    for filename in os.listdir(mydir):
        # print(filename)
        # generate the relevant output file names and paths
        outname = f"genes_annotated_{filename}"
        infile = os.path.join(mydir, filename)
        outfile = os.path.join(mydir, outname)

        # annotate the input csv and output it to csv
        annotate.annotate_standard_csv_input_file(infile, outfile, colname="chr", explode=False)

        additional_filtering(outfile, filename)


if __name__ == '__main__':
    annotate_standard_csv()
