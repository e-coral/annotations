import os
import pandas
from pathlib import Path

# initialise relevant paths
refs_dir = (Path(__file__).parent / 'ref_files').resolve()
outdir = (Path(__file__).parent.parent / 'default_output').resolve()

input_file = 'GENE_LISTS_100124.xlsx'
genes_file = 'chm13v2.0_RefSeq_Liftoff_genes_only.csv'
extras = 'CAT_liftoff_from_table_browser.csv'
nontels_file = 'non-telomere_regions.bed'
centromeres_file = 'centromeres.csv'
ensembl_genes = 'ensembl_IDs_to_gene_names.csv'


def get_telomere_boundaries():
    """
    Read in the bed format (tsv) telomere boundaries file
    column headers = Chr, Start, End
    :return: df of values from telomere boundaries file
    """
    ann_tels = pandas.read_csv(os.path.join(refs_dir, nontels_file), sep='\t')

    # convert all positions to int type to allow relevant operations and removal of trailing .0
    ann_tels = ann_tels.astype({"Start": int, "End": int})

    return ann_tels


def get_centromere_positions():
    """
    Read in the centromere regions file (csv), generated from the cytobands ideogram file by get_centromeres.py
    Parse the file to get the centre of the centromere for each chromosome
    column headers = chrom,chromStart,chromEnd,name,gieStain,ChrArm
    :return: df of values from centromeres file
    """
    centromeres = pandas.read_csv(os.path.join(refs_dir, centromeres_file))
    # not yet sure what is useful, so TODO: parse file to get useful values

    # subset the df to only the p arm, because there is a p and q entry for each chromosome
    centromeres = centromeres[centromeres['ChrArm'] == "p"]

    # geta df of just 'centromere position' (as the end of the first centromere section) for each chromosome
    centromeres = centromeres[['chrom', 'chromEnd']]
    centromeres.columns = ["seqid", "Centromere"]

    # convert position to int type to allow relevant operations and removal of trailing .0
    centromeres = centromeres.astype({"Centromere": int})

    return centromeres


def get_gene_annotations():
    """
    Read in the files relevant to gene annotations:
        genes_df: just the genes from the chm13v2.0_RefSeq_Liftoff.csv, from UCSC
        column headers = '',seqid,source,type,start,end,score,strand,phase,attributes
        extras_df: fuller CAT/Liftoff file (csv), which contains some of the AC... AP... genes missing from genes_df
        column headers = seqid,start,end,name,Gene

    :return: concatenated data frame containing all the gene annotations
    """
    # read the files into dfs
    gene_anns_df = pandas.read_csv(os.path.join(refs_dir, genes_file))
    extra_anns_df = pandas.read_csv(os.path.join(refs_dir, extras))

    # extract the gene name from the attributes string
    gene_anns_df["gene_name"] = gene_anns_df["attributes"].str.extract(r'.*;gene_name=(.*?);.*').fillna('')
    gene_anns_df["gene"] = gene_anns_df["attributes"].str.extract(r'.*;gene=(.*?);.*').fillna('')

    # subset the dfs to the same, relevant columns
    gene_anns_df = gene_anns_df[["seqid", "start", "end", "gene_name", "gene"]]
    extra_anns_df = extra_anns_df[["seqid", "start", "end", "gene_name", "gene"]]

    # concatenate into a single df
    gene_anns = pandas.concat([gene_anns_df, extra_anns_df])

    return gene_anns


def convert_ensembl_ids(genes_df):
    """
    Read in the ensembl IDs and their corresponding gene names, and use this for ensgID > gene name conversion
    (positions cannot be collected from ensembl directly, because it doesn't yet have t2t positions)
    ensg_df: ensembl IDs and gene names, for converting ensembl IDs to gene names
    column headers = no headers
    :param pandas.DataFrame genes_df: initialised dataframe containing input genes
    :return: updated dataframe with new gene names for ensgIDs
    """
    # read in the data
    ensg_df = pandas.read_csv(os.path.join(refs_dir, ensembl_genes), header=None)

    # name the columns
    ensg_df.columns = ["orig_gene", "alt_gene"]

    # merge ID gene names into the input genes df
    genes_df = pandas.merge(genes_df, ensg_df, on="orig_gene", how="left")

    # return the df that now includes additional gene names
    return genes_df


def report_duplicated_genes(genes_df, outfile):
    """
    identify and report duplicated gene names and their locations
    :param genes_df: the df containing gene names and positions from the gene annotation resources
    :param outfile: path to output file for the duplicated genes
    :return: an output file with duplicated genes
    """
    duplicaterows = pandas.concat(g for _, g in genes_df.groupby("gene_name") if len(g) > 1)
    duplicaterows.to_csv(outfile, index=False)


def report_genes_with_no_annotations(genes_df, outfile):
    """
    identify and report genes in the input file that were not found in the annotations files
    :param genes_df: the df containing gene names and positions from the gene annotation resources
    :param outfile: path to output file for the genes with no available position information
    :return: output file
    """
    no_ann = genes_df[genes_df["seqid"].isna()]
    no_ann = no_ann.dropna(axis=1)
    no_ann.to_csv(outfile, index=False)


def find_unique_genes_in_column(genes):
    """
    find the unique gene names in a given column
    :param Pandas.Series genes: df column of gene names
    :return: list of unique gene names
    """
    # remove nan values
    no_nan = genes.dropna()

    # de-duplicate
    unique = list(set(list(no_nan)))

    # remove empty strings
    if "" in unique:
        unique.remove("")

    return unique


def add_all_positions(genes, annotations):
    """
    Merge the input genes and annotations to add the positions to the input genes, and ensure correct formatting of
    positions
    :param genes: input df of gene names
    :param annotations: df of gene annotations
    :return: df containing input genes and their positions, including all variants of gene names
    """
    # merge the dfs of gene names and gene annotations on each of the potential gene name columns
    # df1 will contain most entries
    pos_df1 = pandas.merge(genes, annotations, left_on=['orig_gene'], right_on=['gene'], how='left')

    # alt_gene only contains a few translated ensgIDs, so df2 and df3 will be very small
    pos_df2 = pandas.merge(genes, annotations, left_on=['alt_gene'], right_on=['gene'], how='left')

    # remove rows where gene name and gene are the same to prevent doubling up the merge work
    annotations = annotations[annotations.gene != annotations.gene_name]

    # alt_gene only contains a few translated ensgIDs, so df2 and df3 will be very small
    pos_df3 = pandas.merge(genes, annotations, left_on=['alt_gene'], right_on=['gene_name'], how='left')

    # df4 would have largely contained duplicates of df1, because gene and gene name are often the same
    pos_df4 = pandas.merge(genes, annotations, left_on=['orig_gene'], right_on=['gene_name'], how='left')

    # drop all na rows from the largely empty/duplicated dfs, but not from df1 to maintain input genes
    pos_df2.dropna(subset=['seqid', 'start'], inplace=True)
    pos_df3.dropna(subset=['seqid', 'start'], inplace=True)
    pos_df4.dropna(subset=['seqid', 'start'], inplace=True)

    pos_df4.to_csv(os.path.join(outdir, 'pos_df4.csv'), index=False)

    # concatenate the dfs
    pos_df = pandas.concat([pos_df1, pos_df2, pos_df3, pos_df4], ignore_index=True, sort=False)
    # print(pos_df.seqid.isna().sum())  # only 206 not found now

    return pos_df


def get_plain_gene(df, col_name):
    """
    Get a relevant row from the df of gene and position(s)
    :param Pandas.DataFrame df: df of positions for a single gene
    :param str col_name: "gene_name" or "gene"
    :return: single-row df of gene and position
    """
    # if an alternative name was found, use it
    if not df.alt_gene.isnull().all():
        plain_gene = df[df.alt_gene == df[col_name]]

    else:
        plain_gene = df[df.orig_gene == df[col_name]]

    if plain_gene.empty:
        plain_gene = df

    return plain_gene


def remove_additional_copies(pos_df):
    """
    Remove duplicates and additional copies of gene positions from the df
    :param pos_df:
    :return: clean df of genes and positions
    """
    # initialise the list of dfs of relevant entries
    relevant_entries = []

    # remove duplicates
    pos_df.drop_duplicates(inplace=True)

    # split the df by original gene name, to determine if they are associated with multiple entries
    grouped = pos_df.groupby('orig_gene')
    for group in grouped:  # group[0] = orig_gene value (str), group[1] = df containing all the entries for that gene
        # if there's only one entry for the gene, add it to the final df
        if len(group[1]) == 1:
            relevant_entries.append(group[1])
        # if there's more than one entry, check, de-duplicate, add relevant entries to the final df
        else:
            # count the number of unique gene/gene_name entries (excluding empty/NaN values)
            unique_gene_name = find_unique_genes_in_column(group[1].gene_name)
            unique_gene = find_unique_genes_in_column(group[1].gene)

            # if there's only one unique name in one of the columns, then test whether there's a single value in the
            # other column that exactly matches the gene name (orig or alt).
            # If so, add that single row. Else, add all rows.
            if len(unique_gene_name) == 1:
                plain_gene = get_plain_gene(group[1], "gene")

            elif len(unique_gene) == 1:
                plain_gene = get_plain_gene(group[1], "gene_name")

            else:
                plain_gene = pandas.DataFrame({})
                if len(unique_gene_name) > 1 and len(unique_gene) > 1:
                    print(f"non-na unique gene names in the group for {group[0]}")
                    group[1].to_csv(os.path.join(outdir, f"{group[0]}-non-unique.csv"))
                else:
                    print(f"The gene names column contained {len(unique_gene_name)} unique entries,"
                          f"and the gene column contained {len(unique_gene)} entries")
                    group[1].to_csv(os.path.join(outdir, f"{group[0]}-multiple_entries.csv"))

            if not plain_gene.empty:
                relevant_entries.append(plain_gene)

    rel_ents_df = pandas.concat(relevant_entries, ignore_index=True)

    rel_ents_df.to_csv(os.path.join(outdir, 'relevant_entries.csv'), index=False)

    return rel_ents_df


def add_position_annotations(genes, annotations):
    """
    Merge the input genes and annotations to add the positions to the input genes, and ensure correct formatting of
    positions
    :param genes: input df of gene names
    :param annotations: df of gene annotations
    :return: df containing input genes and their positions
    """
    pos_df = add_all_positions(genes, annotations)

    relevant_entries = remove_additional_copies(pos_df)

    return relevant_entries


def add_tel_positions(genes):
    """
    Merge the genes and positions with the telomere boundary positions
    :param genes: df containing input genes and their positions
    :param tels: df of telomere boundary positions
    :return: merged df, containing genes, gene positions, and telomere boundary positions
    """
    tels = get_telomere_boundaries()  # no trailing .0
    tels_added = pandas.merge(genes, tels, left_on=['seqid'], right_on=['Chr'], how='left')  # trailing .0 returns

    # print(tels_added.columns)
    return tels_added


def add_cent_positions(genes):
    """
    Merge the genes and positions with the centromere positions
    :param genes: df of genes and their positions
    :return: merged df, containing genes, gene positions, and centromere positions
    """
    cents = get_centromere_positions()
    cents_added = pandas.merge(genes, cents, on=['seqid'], how='left')
    # print(cents_added.columns)
    return cents_added


def calculate_gene_lengths(df):
    """
    Calculate the lengths of the genes by end position - start position
    :param df: the dataframe containing genes and their positions
    :return: the dataframe with added gene lengths
    """
    df['gene_length'] = df['end'] - df['start']
    return df


def calculate_distances_to_telomeres(df):
    """
    Calculate the distances between genes and the telomere on the same chromosome arm
    :param df: the dataframe containing the gene info
    :return: df containing telomere positions and distances between genes and telomere boundaries
    """
    # temporarily remove na values to allow successful conversion for calculations
    no_na_df = df.dropna(subset=['seqid'])

    # add the required telomere and centromere positions
    no_na_df = add_tel_positions(no_na_df)
    no_na_df = add_cent_positions(no_na_df)

    # temporarily remove na values to allow successful conversion for calculations
    na_cents = no_na_df[no_na_df['Centromere'].isna()]
    no_na_df = no_na_df.dropna(subset=['Centromere'])

    # generate a df of the genes for which no annotations were found
    not_annotated = pandas.concat([df[df['seqid'].isna()], na_cents])

    # convert all positions to int type to allow relevant operations and removal of trailing .0
    no_na_df = no_na_df.astype({"start": int, "end": int, "Start": int, "End": int,
                                "gene_length": int, "Centromere": int})

    # if the gene is on the p arm, calculate the distance to the p telomere
    no_na_df.loc[(no_na_df["start"] < no_na_df["Centromere"]), "g-t_distance"] = no_na_df["start"] - no_na_df["Start"]
    no_na_df.loc[(no_na_df["start"] > no_na_df["Centromere"]), "g-t_distance"] = no_na_df["End"] - no_na_df["end"]

    # convert the integers to strings to prevent trailing .0s reappearing after NaN values are added back in
    no_na_df = no_na_df.astype({"g-t_distance": int})
    no_na_df = no_na_df.astype({"start": str, "end": str, "Start": str, "End": str, "gene_length": str,
                                "Centromere": str, "g-t_distance": str})

    # add the removed rows back in
    df = pandas.concat([no_na_df, not_annotated])
    # df.to_csv(os.path.join(outdir, "mostly_annotated.csv"), index=False)

    return df


def add_gene_location(df):
    """
    Adds a column containing gene location in chrX:123-456 format
    :param df: df of position-annotated genes
    :return: the same df with an additional column containing the formatted location of the gene
    """
    df['full_gene_loc'] = df.agg(lambda x: f"{x['seqid']}:{x['start']}-{x['end']}", axis=1)
    return df


def reformat_for_output(df):
    """
    organise and rename the columns of the annotated df
    :param df: df of genes and annotations
    :return: reformatted df
    """
    # keep only relevant info
    df = df[['orig_gene', 'alt_gene', 'gene', 'gene_name', 'seqid', 'start', 'end',
             'gene_length', "g-t_distance", 'full_gene_loc']]

    # rename columns
    df.columns = ['Orig_gene', 'Alt_gene', 'Annotation_gene', 'Annotation_gene2', 'Chr', 'Start',
                  'End', 'Gene_length', "Distance_to_telomere", 'Gene_loc']

    return df


def output_full_csv(df):
    """
    output the df as a single csv file
    :param df: df of genes and annotations
    :return: output file
    """
    # output the df to csv
    df.to_csv(os.path.join(outdir, 'full_gene_distance_output.csv'), index=False)


def recreate_multisheet_excel_doc(df):
    """
    organise and rename the columns of the annotated df, and output as a single csv file
    :param df: df of genes and annotations
    :return: output file
    """
    # read the file in
    in_genes = pandas.read_excel(os.path.join(refs_dir, input_file))

    with pandas.ExcelWriter(os.path.join(outdir, "original_sheets_plus_distances.xlsx")) as writer:
        for _, s in in_genes.items():
            s = s.dropna()
            s_df = pandas.DataFrame(s.values, columns=['Orig_gene'])
            tmp_df = pandas.merge(s_df, df, how='left', on='Orig_gene')
            # print(tmp_df.columns)

            # reorder again
            tmp_df = tmp_df[['Orig_gene', 'Alt_gene', 'Annotation_gene', 'Annotation_gene2', 'Chr', 'Start', 'End', 'Gene_length', 'Distance_to_telomere', 'Gene_loc']]

            # print(tmp_df.head())
            tmp_df.to_excel(writer, sheet_name=s.name, index=False)


def input_file_to_dataframe():
    """
    Read in the input file and convert it into a df

    :return: df containing input file gene names
    """
    # read the file in, and add column names for ensembl genes
    in_genes = pandas.read_excel(os.path.join(refs_dir, input_file))

    # convert the excel file values into a single, unique list of gene names
    unique_genes = list(set([y.strip() for x in in_genes.values.tolist() for y in x if pandas.notna(y)]))

    # use the unique genes to initialise a dataframe for output, and strip spaces from the gene names just in case
    genes_df = pandas.DataFrame(unique_genes, columns=['orig_gene'])

    # no longer required: strip on y above does the same
    # genes_df = genes_df.map(lambda x: x.strip() if isinstance(x, str) else x)

    return genes_df








#
#
# # print(len(genes_df))  # should be 11210
# # print(genes_df.columns)
#
# # convert ensg IDs
# # malformed biomart xml warning sometimes - skip if so
# genes_df['alt_gene'] = genes_df['orig_gene'].apply(lambda x: ens_dict.get(x))
# genes_df['alt_gene'] = ""
#
# genes_df['Gene'] = numpy.where(genes_df['alt_gene'], genes_df['alt_gene'], genes_df['orig_gene'])
# # genes_df.to_csv(os.path.join(outdir, "intermediate_genesdf2.csv"), index=False)
# # print(genes_df.head())
# # print(len(genes_df))
#
# # in the gene annotations file, the gene name is hidden within the attributes column of the file
# # therefore, extract this value into its own df column
# # there are both gene and gene_name, which are almost identical - so, run with both
# ann_genes["Gene"] = ann_genes["attributes"].str.extract(r'.*;gene_name=(.*?);.*').fillna('')
# # ann_genes["Gene"] = ann_genes["attributes"].str.extract(r'.*;gene=(.*?);.*').fillna('')
# # ann_genes.to_csv(os.path.join(outdir, "intermediate_anngenes_full.csv"), index=False)
#
# # ann_genes.to_csv(os.path.join(outdir, "anns_by_gene.csv"), index=False)
# # ann_genes.to_csv(os.path.join(outdir, "anns_by_genename.csv"), index=False)
# df = pandas.merge(genes_df, ann_genes, on=['Gene'], how='left')
# # print(len(df))  # should be at least 11210 - though some have duplicates, and strip could/should assist the merge
#
#
# # # there are an additional 57 rows after the merge, suggesting multiple entries for some gene names in ann_genes
# # # find them and put them in an output file
# duplicaterows = pandas.concat(g for _, g in df.groupby("Gene") if len(g) > 1)
# # duplicaterows.to_csv(os.path.join(outdir, 'genes_with_multiple_entries_by_gene.csv'), index=False)
# duplicaterows.to_csv(os.path.join(outdir, 'genes_with_multiple_entries_by_genename.csv'), index=False)
# # # print(set(duplicaterows['Gene']))
#
# # ultimately ignore the duplicates for now; just get the distances for all the different positions
#
# # merge the dfs of gene annotations and telomere positions
# three_df = pandas.merge(df, ann_tels, left_on=['seqid'], right_on=['Chr'], how='left')
# # print(three_df.head())
#
# # remove the trailing .0s while preserving nan entries
# # split the data into genes with and without positional info
# no_na = three_df.dropna(subset=['seqid', 'Start'])
# # print(len(no_na))  # should be at least 10575
#
# nan_df = three_df[three_df['seqid'].isna()]
# nan_df = nan_df.dropna(axis=1)
# # nan_df.to_csv(os.path.join(outdir, 'genes_not_in_refseq_liftoff_by_gene.csv'), index=False)
# nan_df.to_csv(os.path.join(outdir, 'genes_not_in_refseq_liftoff_by_genename.csv'), index=False)
# # print(len(nan_df))  # should be at most 691
#
# # for the genes not in refseq/liftoff, try using the CAT table file to identify positions
# filled_gaps = pandas.merge(nan_df, extra_anns, how='left', on='Gene')
#
# # get the genes that match with name2 rather than with name2 (not many)
# filled_alt = pandas.merge(nan_df, extra_anns, how='left', left_on='Gene', right_on='name')
# filled_alt = filled_alt.dropna(subset=['seqid'])
# filled_alt = filled_alt[["orig_gene", "alt_gene", "Gene_x", "seqid", "start", "end", "name"]]
# filled_alt.columns = ['orig_gene', 'alt_gene', 'Gene', 'seqid', 'start', 'end', "name"]
# filled_alt.to_csv(os.path.join(outdir, 'filled_alt.csv'))
# # exit(filled_alt.columns)
#
# # assign the name as the alternative name
# filled_gaps['alt_gene'] = filled_gaps['name']
# print(filled_gaps.columns)
#
# # concatenate the dfs with different match columns
# filled_gaps = pandas.concat([filled_gaps, filled_alt])
# # merge the dfs of gene annotations and telomere positions
# filled_gaps = pandas.merge(filled_gaps, ann_tels, left_on=['seqid'], right_on=['Chr'], how='left')
#
# # split the file into annotated and unannotated for further analysis and output
# filled_gaps_no_na = filled_gaps.dropna(subset=['seqid', 'Start'])
# not_found_in_any_files = filled_gaps[filled_gaps['seqid'].isna()]
# # filled_gaps.to_csv(os.path.join(outdir, 'test_filled_gaps.csv'), index=False)
# not_found_in_any_files.to_csv(os.path.join(outdir, 'genes_not_found_in_any_ref_files.csv'), index=False)
#
# # subset to only the useful columns
# filled_gaps_no_na = filled_gaps_no_na[['orig_gene', 'alt_gene', 'Gene', 'seqid', 'start', 'end', "Start", "End"]]
# no_na = no_na[['orig_gene', 'alt_gene', 'Gene', 'seqid', 'start', 'end', "Start", "End"]]
# no_na = pandas.concat([filled_gaps_no_na, no_na])
# print(no_na.columns)
#
# # convert the values, so they don't have .0 at the end
# # print(no_na.dtypes)
# no_na = no_na.astype({"start": int, "end": int, "Start": int, "End": int})
#
# # calculate the gene length
# no_na['Gene_length'] = no_na['end'] - no_na['start']
#
# no_na['full_gene_loc'] = no_na.agg(lambda x: f"{x['seqid']}:{x['start']}-{x['end']}", axis=1)
# # calculate the differences between the gene start/end ('start' and 'end' in lowercase) and positions where the
# # telomeres begin ('Start' and 'End' in title case)
# no_na['gs-cs'] = no_na['start'] - no_na['Start']
# no_na['ge-cs'] = no_na['end'] - no_na['Start']
# no_na['gs-ce'] = no_na['End'] - no_na['start']
# no_na['ge-ce'] = no_na['End'] - no_na['end']
# no_na['gene-telomere_distance'] = no_na[['gs-cs', 'ge-cs', 'gs-ce', 'ge-ce']].min(axis=1)
# print(no_na.columns)
#
# # create the usual location string in a final column
# no_na['full_gene_loc'] = no_na.agg(lambda x: f"{x['seqid']}:{x['start']}-{x['end']}", axis=1)
#
# # subset to useful columns
# no_na = no_na[['orig_gene', 'alt_gene', 'Gene', 'seqid', 'start', 'end', 'Gene_length', "gene-telomere_distance", 'full_gene_loc']]
# # rename columns
# no_na.columns = ["Orig_Gene", "Alt_Gene", "Gene", "Chr", "Start", "End", 'Gene_length', "Distance_to_telomere", "Gene_loc"]
#
# # output to file (all genes)
# # no_na.to_csv(os.path.join(outdir, 'genes_and_distances_by_gene.csv'), index=False)
# no_na.to_csv(os.path.join(outdir, 'genes_and_distances_by_genename.csv'), index=False)
#
# # # output to file (for each list)
# with pandas.ExcelWriter(os.path.join(outdir, "per-list_allgenes_and_distances_by_genename.xlsx")) as writer:
# # with pandas.ExcelWriter(os.path.join(outdir, "per-list_allgenes_and_distances_by_gene")) as writer:
#     for _, s in in_genes.items():
#         s = s.dropna()
#         # tmp_df = no_na[no_na['Gene'].isin(s)]  # this drops kate's genes if they aren't in the df
#         s_df = pandas.DataFrame(s.values, columns=['Orig_Gene'])
#         tmp_df = pandas.merge(s_df, no_na, how='left', on='Orig_Gene')
#         # print(tmp_df.columns)
#
#         # reorder again
#         tmp_df = tmp_df[['Orig_Gene', 'Alt_Gene', 'Gene', 'Chr', 'Start', 'End', 'Gene_length', 'Distance_to_telomere', 'Gene_loc']]
#
#         # print(tmp_df.head())
#         tmp_df.to_excel(writer, sheet_name=s.name, index=False)
#
#
#     # print(len(s))
#     # list_genes = pandas.merge(s, no_na, how='left')
#     # print(list_genes)

# start = ann_tels.loc[ann_tels['Chr'] == 'chr1', 'Start'].values[0]
# end = ann_tels

# test1 = df['start'][2] - ann_tels.loc[ann_tels['Chr'] == df['seqid'][2], ['Start']]
# df['gs-cs'] = df['start'] - ann_tels.loc[ann_tels['Chr'] == df['seqid'].values, 'Start']
# print(test1)
#
# print(df.columns)

# # make the gene annotations file usable, based on make_annotation_ncls in eccDNA pipelines
# # initialise
# regions = {}
# # split by chromosomes
# ann_genes_by_chrom = {k: v for k, v in ann_genes.groupby("seqid")}
# print(type(ann_genes_by_chrom['chr1'].iloc[2]))
# print(ann_genes_by_chrom['chr1'].iloc[2])
# # for k, v in ann_genes_by_chrom.items():
# #     regions[k] = NCLS(v.start.values, v.end.values, v.index.values)
# # gene_regions = get_annotation_regions(self.gene_df)
#
#
#
# # print(ann_genes.columns)
# # print(ann_tels)
#

