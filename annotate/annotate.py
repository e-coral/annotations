import os
import pandas
import numpy
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
    # no editing required, because the file is very minimal
    return ann_tels


def get_centromere_positions():
    """
    Read in the centromere regions file (csv), generated from the cytobands ideogram file by _____
    column headers = chrom,chromStart,chromEnd,name,gieStain,ChrArm
    :return: df of values from centromeres file
    """
    centromeres = pandas.read_csv(os.path.join(refs_dir, centromeres_file))
    # not yet sure what is useful, so TODO: parse file to get useful values
    return centromeres


def get_gene_annotations():
    """
    Read in the files relevant to gene annotations:
        genes_df: just the genes from the chm13v2.0_RefSeq_Liftoff.csv, from UCSC
        column headers = '',seqid,source,type,start,end,score,strand,phase,attributes
        extras_df: fuller CAT/Liftoff file (csv), which contains some of the AC... AP... genes missing from genes_df
        column headers = seqid,start,end,name,Gene

    :return:
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


def add_position_annotations(genes, annotations):
    """
    Merge the input genes and annotations to add the positions to the input genes, and ensure correct formatting of
    positions
    :param genes: input df of gene names
    :param annotations: df of gene annotations
    :return: df containing input genes and their positions
    """
    # merge the dfs of gene names and gene annotations on each of the potential gene name columns
    # print(genes.columns)
    # print(annotations.columns)

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

    # remove the duplicate rows
    pos_df.drop_duplicates(inplace=True)
    # print(pos_df.seqid.isna().sum())  # only 206 not found now

    # no_alt_transcripts = pos_df[pos_df.gene_name == pos_df.gene]
    # no_alt_transcripts = no_alt_transcripts[no_alt_transcripts.gene_name == no_alt_transcripts.gene]

    # no_alt_transcripts.to_csv(os.path.join(outdir, 'no_alt_transcripts.csv'), index=False)
    # pos_df.to_csv(os.path.join(outdir, 'final_pos_df.csv'))

    relevant_entries = []

    grouped = pos_df.groupby('orig_gene')
    for group in grouped:
        # group[0] is orig_gene value (str), group[1] is a df containing all the entries for that gene
        # print(type(group[0]))
        if len(group[1]) == 1:
            # print(group[1])
            relevant_entries.append(group[1])
        else:
            # count the number of unique gene/gene_name entries
            same_gene_name = group[1].gene_name.nunique()
            same_gene = group[1].gene.nunique()

            if same_gene_name == 1:
                # print(f"unique g n")
                plain_gene = group[1][group[1].orig_gene == group[1].gene]
                # print(f"unique gene name, length: {len(plain_gene)}")

            elif same_gene == 1:
                # print(f"unique g")
                plain_gene = group[1][group[1].orig_gene == group[1].gene_name]
                # print(f"unique g, len(plain_gene) = {len(plain_gene)}")

            else:
                # if it's just one gene and NaN in a column, then trim off the names and add them as for above
                same_gene_name = group[1].gene_name.dropna().nunique()
                same_gene = group[1].gene.dropna().nunique()
                if not same_gene == 1 or same_gene_name == 1:
                    print(f"non-na unique gene names in the group for {group[0]}")
                    group[1].to_csv(os.path.join(outdir, f"{group[0]}un-unique.csv"))
                    plain_gene = pandas.DataFrame({})
                else:
                    pass

            if not plain_gene.empty:
                relevant_entries.append(plain_gene)

    rel_ents_df = pandas.concat(relevant_entries, ignore_index=True)

    rel_ents_df.to_csv(os.path.join(outdir, 'relevant_entries.csv'), index=False)
        # print(len(group[1]))
        # print(group[2])
        # print(len(group[1]))
        # exit()

    return pos_df


def add_tel_lengths(genes, tels):
    """
    Merge the genes and positions with the telomere boundary positions
    :param genes: df containing input genes and their positions
    :param tels: df of telomere boundary positions
    :return:
    """


def set_gene_name_for_annotation():
    pass


def input_file_to_dataframe():
    """
    Read in the input file and convert it into a df

    :return:
    """
    # read the files in, and add column names for ensembl genes
    in_genes = pandas.read_excel(os.path.join(refs_dir, input_file))

    # convert the excel file values into a single, unique list of gene names
    unique_genes = list(set([y for x in in_genes.values.tolist() for y in x if pandas.notna(y)]))
    # print(len(unique_genes))

    # use the unique genes to initialise a dataframe for output, and strip spaces from the gene names just in case
    genes_df = pandas.DataFrame(unique_genes, columns=['orig_gene'])
    genes_df = genes_df.map(lambda x: x.strip() if isinstance(x, str) else x)

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

