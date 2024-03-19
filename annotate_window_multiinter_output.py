import os
from pathlib import Path
from annotate import annotate
import pandas
import glob


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
                'genes', 'gene_lengths', 'g-t_distance']:
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
    for infile in glob.glob('C:/Users/ec339/Downloads/window_analysis/window_analysis/multiintersect_outputs/just_results/*'):
        outname = f"annotated_{os.path.basename(infile)}"
        winsize = outname.split("_")[1].strip("kb")
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
        res = annotate.calculate_distances_to_telomeres(res)

        # annotate the genes, fragile sites, repeats and gene sizes
        regions_df = annotate.annotate_overlaps(res)

        # reformat the final df to match the input
        final = format_columns(regions_df, orig_columns)

        # split out the genes into separate columns, if preferred
        # final = annotate.split_multiple_genes(final)  # not preferred for in_both!

        # write output file
        # final.to_csv(os.path.join(outdir, f"{outname}-exploded.csv"), index=False)
        final.to_csv(os.path.join(outdir, f"plain_{outname}.csv"), index=False)


if __name__ == '__main__':
    annotate_results()
