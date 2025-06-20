# keep the file collection and y/n sections of this - but move to the main annotations module


import os
from pathlib import Path
from src.annotate import annotate
import pandas
import glob
import numpy as np

outdir = (Path(__file__).parent / 'multiinter_output').resolve()
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

    df['any_fragile?'] = np.where(df['fragile_sites'] == "", "n", "y")
    df['any_repeats?'] = np.where(df['repeats'] == "", "n", "y")

    # drop the telomere and centromere columns that aren't needed in the final output
    new_df = df.drop(columns=['tChr', 'tStart', 'tEnd', 'Centromere', 'fragile_sites', 'repeats'])

    # set values to strings to prevent trailing .0s
    df = new_df.astype(str)

    return df


def annotate_results():
    """
    read in the file, add the annotations, and output the annotated file
    :return: annotated file
    """
    for infile in glob.glob('C:/Users/ec339/Downloads/capped_with_corrections/for_annotating/*'):
        outname = f"annotated_{os.path.basename(infile)}"
        winsize = outname.split("_")[2].strip("kb")
        if winsize.endswith("m"):
            winsize = 1000000
        else:
            winsize = int(winsize) * 1000

        # read in the data
        res = pandas.read_csv(infile, header=0, sep='\t')

        # replace the window start and end with the original junction position start and end
        # res.rename(columns={"start": "window_start"}, inplace=True)
        # res.rename(columns={"end": "window_end"}, inplace=True)

        # calculate the original position around which the window was constructed
        # #(do not use start, as this can be rationalised to 0 if start < window)
        # res['start'] = res['window_end'] - winsize
        # res['end'] = res['window_end'] - winsize

        # store the original column names to be able to add them back in later
        orig_columns = [col.strip() for col in res.columns]

        # rename the chr column, because some of the methods may be based on 'seqid'
        res.rename(columns={"chrom": "seqid"}, inplace=True)

        # calculate the distance between the region and the nearest telomere
        res = annotate.calculate_distances_to_centromeres_and_telomeres(res)

        # annotate the genes, fragile sites, repeats and gene sizes
        regions_df = annotate.annotate_overlaps(res)

        # reformat the final df to match the input
        final = format_columns(regions_df, orig_columns)

        # split out the genes into separate columns, if preferred
        # final = annotate.split_multiple_genes(final)  # not preferred for in_both!

        # write output file
        # final.to_csv(os.path.join(outdir, f"{outname}-exploded.csv"), index=False)
        final.to_csv(os.path.join(outdir, f"only_censat_{outname}.csv"), index=False)


if __name__ == '__main__':
    annotate_results()
