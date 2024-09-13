# irrelevant following cleanup; use annotate_csv instead

import os
from pathlib import Path
from src.annotate import annotate
import pandas


outdir = (Path(__file__).parent / 'kas-seq_output').resolve()
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
                'genes', 'gene_lengths', 'g-t_distance']:
        orig_col_names.append(col)

    # rename the columns accordingly
    df.columns = orig_col_names

    # drop the telomere and centromere columns that aren't needed in the final output
    df = df.drop(columns=['tChr', 'tStart', 'tEnd', 'Centromere'])

    # set values to strings to prevent trailing .0s
    df = df.astype(str)

    return df


def annotate_kasseq_results():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    # for each of the ouput files
    for filename in refs_dir.glob(r'kas-seq_results/*'):
        print(filename)

        # create the output file name
        outname = f"pytest_annotated_{filename.stem}"

        # read in the data
        res = pandas.read_csv(filename, header=0)

        # store the original column names to be able to add them back in later
        orig_columns = [col.strip() for col in res.columns]

        # rename the seqnames column, because some of the methods may be based on 'seqid'
        res.rename(columns={"seqnames": "seqid"}, inplace=True)

        # calculate the distance between the region and the nearest telomere
        res = annotate.calculate_distances_to_telomeres(res)

        # annotate the genes, fragile sites, repeats and gene sizes
        regions_df = annotate.annotate_overlaps(res)

        # reformat the final df to match the input
        final = format_columns(regions_df, orig_columns)

        # split out the genes into separate columns, if preferred
        # final = annotate.split_multiple_genes(final)

        # write output file
        final.to_csv(os.path.join(outdir, f"{outname}-notexploded.csv"), index=False)


if __name__ == '__main__':
    annotate_kasseq_results()
