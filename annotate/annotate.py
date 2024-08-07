import os
import pandas
from pathlib import Path
from ncls import NCLS
import re

# initialise relevant paths
refs_dir = (Path(__file__).parent / 'ref_files').resolve()
outdir = (Path(__file__).parent.parent / 'default_output').resolve()

input_file = 'GENE_LISTS_100124.xlsx'
genes_file = 'chm13v2.0_RefSeq_Liftoff_genes_only.csv'
extras = 'CAT_liftoff_from_table_browser.csv'
nontels_file = 'non-telomere_regions.bed'
centromeres_file = 'centromeres.csv'
ensembl_genes = 'ensembl_IDs_to_gene_names.csv'

# for region annotations, in addition to genes_file above
# reps_file = 'chm13v2.0_rmsk.bed'
reps_file = 'censat.bed'
fs_file = 'all_fragile_sites_positions.tsv'

# set the df printing width to allow viewing all columns
pandas.set_option('display.max_columns', None)


def get_separate_chr_pos(df, colname):
    """
    if positions are supplied in chr:start-end format, separate these values into separate columns
    :param df: dataframe containing chr:start-pos format positions
    :param str colname: name of the relevant column containing chr:start-pos data
    :return: dataframe containing separate columns for chr, start and end
    """
    # create the pattern on which to split the data
    pat = re.compile(r':(\d+)-')

    # remove commas from the numbers, and any trailing spaces
    df[colname] = df[colname].str.replace(',', '')
    df[colname] = df[colname].str.strip()

    # split into seqid, start and end based on : and -
    df[["seqid", "start", "end"]] = df[colname].str.split(pat, expand=True)

    # return the edited df
    return df


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
    Read in the files relevant to gene position annotations:
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


def get_genes_df():
    """
    read the genes file into a relevant df
    :return: df of gene annotations
    """
    genes_df = pandas.read_csv(os.path.join(refs_dir, genes_file))
    genes_df = calculate_gene_lengths(genes_df)
    genes_df = calculate_distances_to_telomeres(genes_df, keep_int=True, isgene=True)

    return genes_df.astype({'gene_length': str, 'g-t_distance': str})


def get_fs_df():
    """
    read the fragile sites file into a relevant df
    :return: df of fragile site annotations
    """
    fs_df = pandas.read_csv(os.path.join(refs_dir, fs_file), sep='\t', header=0)
    return fs_df


def get_reps_df():
    """
    read the reps file into a relevant df
    :return: df of reps annotations
    """
    reps_df = pandas.read_csv(os.path.join(refs_dir, reps_file), sep='\t',
                              names=["chrom", "start", "end", "repeats", "score", "strand", "thickStart", "thickEnd",
                                     "reserved", "blockCount", "blockSizes", "blockStarts", "id", "description"])

    return reps_df


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

    # pos_df4.to_csv(os.path.join(outdir, 'pos_df4.csv'), index=False)

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


def remove_additional_copies(pos_df, only_longest_transcript = False):
    """
    Remove duplicates and additional copies of gene positions from the df
    :param pos_df: df containing gene positions with possible duplicates
    :param bool only_longest_transcript: whether to only keep rows for the longest transcript
    :return: clean df of genes and positions
    """
    # initialise the list of dfs of relevant entries
    relevant_entries = []

    # remove duplicates
    pos_df.drop_duplicates(inplace=True)

    if only_longest_transcript:
        pos_df = calculate_gene_lengths(pos_df)

    # split the df by original gene name, to determine if they are associated with multiple entries
    grouped = pos_df.groupby('orig_gene')
    for group in grouped:  # group[0] = orig_gene value (str), group[1] = df containing all the entries for that gene
        # if there's only one entry for the gene, add it to the final df
        if len(group[1]) == 1:
            relevant_entries.append(group[1])
        # if there's more than one entry, check, de-duplicate, add relevant entries to the final df
        else:
            if only_longest_transcript:
                plain_gene = group[1][group[1]['gene_length'] == group[1]['gene_length'].max()]
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

    relevant_entries = remove_additional_copies(pos_df, only_longest_transcript=True)

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


def get_regions_per_chr(df):
    """
    to find genes etc. that overlap with the regions of interest, the regions must be split by chromosome
    to allow for number overlap searches
    :param df: df of the input file
    :return: dict with chromosomes as keys and the rest of the df values as values
    """
    if "seqid" in df.columns:
        regions_by_chr = {k: v for k, v in df.groupby('seqid')}
    elif "chr" in df.columns:
        regions_by_chr = {k: v for k, v in df.groupby('chr')}
    elif "chrom" in df.columns:
        regions_by_chr = {k: v for k, v in df.groupby('chrom')}
    else:
        exit(f"unable to find chromosome column in the columns:\n{df.columns}")

    return regions_by_chr


def find_reps_overlaps(chrom, pos, reps, reps_regions, reps_df, posend=None):
    """
    find positions of variant junctions that overlap with repeats regions
    (similar to methods for censat and fsites, but each data structure requires its own specific df access command)
    :param chrom: chromosome number of a variant
    :param pos: position of the junction of a variant
    :param list reps: list of annotated variants to be added to
    :param reps_regions: repeats regions dict
    :param reps_df: df of rep regions
    :param posend: position of the end of the variant junction, if specified (otherwise, set to be pos + 1)
    :return: list of variant junction-overlapping repeats regions for annotation
    """
    # initialise list
    rep = []
    if posend is None:
        posend = pos + 1

    try:
        # use the ncls find_overlap method to check whether the var positions overlap with repeat regions
        reps_overlaps = list(reps_regions[chrom].find_overlap(pos, posend))
        # if there are any overlaps,
        # then for each overlap in the list, use the index of the region to identify the relevant annotation
        if reps_overlaps:
            for item in reps_overlaps:  # usually only 1 pos-annotation overlap, but it's possible for more
                # use the index to find the relevant annotation(s), and add them to the list for the function
                rep.append(reps_df.iloc[item[2]].repeats)

    except KeyError:
        print(f"No repeats annotations available for {chrom}.")

    # join the function list (even if empty) into a comma-separated string and add it to the annotations list
    reps.append(", ".join(rep))

    return reps


def find_fs_overlaps(chrom, pos, f_sites, fs_regions, fs_df, posend=None):
    """
    find positions of variant junctions that overlap with fragile sites regions
    (similar to methods for censat and reps, but each data structure requires its own specific df access command)
    :param chrom: chromosome number of a variant
    :param pos: position of the junction of a variant
    :param list f_sites: list of annotated variants to be added to
    :param dict fs_regions: the dict of NCLSs for the of fragile sites regions
    :param pandas.DataFrame fs_df: fragile sites dataframe
    :param posend: if using an end position, add here. else, None, so that it can be set to start + 1
    :return: list of variant junction-overlapping fragile sites regions for annotation
    """
    # initialise list
    fs = []
    if posend is None:
        posend = pos + 1

    # it won't work for chromosomes Y or M, because there are no data for those
    if chrom not in ["chrM", "chrY"]:
        try:
            # use the ncls find_overlap method to check whether the var positions overlap with repeat regions
            fs_overlaps = list(fs_regions[chrom].find_overlap(pos, posend))
            # if there are any overlaps,
            # then for each overlap in the list, use the index of the region to identify the relevant annotation
            if fs_overlaps:
                for item in fs_overlaps:  # usually only 1 pos-annotation overlap, but it's possible for more
                    # use the index to find the relevant annotation(s), and add them to the list for the function
                    fs.append(fs_df.iloc[item[2]].fragile_site)

        except KeyError:
            print(f"No fragile sites annotations available for {chrom}:{pos}-{posend}.")

    # join the function list (even if empty) into a comma-separated string and add it to the annotations list
    f_sites.append(", ".join(fs))

    return f_sites


def find_gene_overlaps(chrom, pos, gene_names, gene_regions, gene_df, gene_sizes, gtds, gene_positions, posend=None):
    """
    find positions of variant junctions that overlap with genes regions
    :param chrom: chromosome number of a variant
    :param pos: position of the junction of a variant
    :param list gene_names: list of gene name annotations to be added to
    :param gene_regions: NCLS of the gene regions
    :param gene_df: dataframe containing gene data for the chromosome
    :param list gene_sizes: list of gene size annotations to be added to
    :param list gtds: list of gene-telomere distances to be added to
    :param list gene_positions: list of positions of genes to be added to
    :param posend: position of the end of the variant, if supplied (otherwise, set to be start pos + 1)
    :return: list of variant junction-overlapping gene sites regions for annotation
    """
    gtd = []
    gene_length = []
    gene_name = []
    gene_pos = []

    if posend is None:
        posend = pos + 1

    try:
        genes_overlaps = list(gene_regions[chrom].find_overlap(pos, posend))
        if genes_overlaps:
            for item in genes_overlaps:
                relevant_data = gene_df.iloc[item[2]]

                # if the type of the annotation is a gene, then get the gene name and size
                if relevant_data.type == "gene":
                    gene_length.append(relevant_data.gene_length)
                    gtd.append(relevant_data["g-t_distance"])
                    gene_pos.append(f"{relevant_data['seqid']}:{relevant_data['start']}-{relevant_data['end']}")
                    # extract the gene name from the attributes field
                    atts = relevant_data["attributes"]
                    # print(atts)
                    if match := re.search(r'.*;gene_name=(.*?);.*', atts):
                        # print(atts)
                        # print(match.group(1))
                        gene_name.append(match.group(1))
                    else:
                        gene_name.append("")

    except KeyError as err:
        print(f"No gene annotations available for {chrom}:{pos}-{posend}.")

    gene_sizes.append(", ".join(gene_length))
    gene_names.append(", ".join(gene_name))
    gtds.append(", ".join(gtd))
    gene_positions.append(", ".join(gene_pos))

    return gene_names, gene_sizes, gtds, gene_positions


def add_annotation_column(df, data_list, column_name):
    """
    validate the length of the data to be added, and add it
    :param pandas.DataFrame df: the by-chromosome subset of the dysgu results df to be modified
    :param list data_list: list of new information to add
    :param str column_name: name for the new column
    :return pandas.DataFrame: edited df
    """
    # if the list isn't empty, and is the length of the df, then add it to the df as a column
    if data_list:
        if len(df) == len(data_list):
            df[column_name] = data_list
        else:
            # if the lengths differ, something has gone wrong
            print(f"The length of the annotations list for '{column_name}' does not equal the length "
                  f"of the variants csv.")
            print(f"Continuing without '{column_name}' annotations")
            df[column_name] = ""
    # if the list is empty (e.g. chrM), just add an empty column
    else:
        df[column_name] = ""

    return df


# def find_overlapping_features_for_dysgu_output(dysgu_data):
#     """
#     annotate each of the regions/boundaries in the dysgu output with features that overlap them
#     # the dysgu output has two chrom columns and two position columns representing variant boundaries
#
#     :return:
#     """
#
#     # create the dfs of annotations and NCLSs of regions covered by features to annotate
#     reps_df = get_reps_df()
#     genes_df = get_genes_df()
#     fs_df = get_fs_df()
#
#     # create the NCLSs for the repeats, genes and fragile sites
#     rep_regions, genes_regions, fs_regions = make_annotation_ncls(reps_df, genes_df, fs_df)
#
#     # split the data by chromosome
#     vars_by_chrom = {k: v for k, v in dysgu_data.groupby("chrA")}
#
#     # for each chromosome, annotate the positions of the variants
#     for k, v in vars_by_chrom.items():
#         # initialise the lists of annotations
#         a_repeats = []
#         b_repeats = []
#         a_f_sites = []
#         b_f_sites = []
#         a_genes = []
#         b_genes = []
#         a_gene_lengths = []
#         b_gene_lengths = []
#         a_gtd = []
#         b_gtd = []
#         a_gene_positions = []
#         b_gene_positions = []
#
#         for i, r in v.iterrows():  # i = df index (int), r is a df row
#             # get the start position for chrA (chrA was used to split the df, so = k)
#             pos_a = r.posA
#
#             # get the chr and start position for chrB
#             chr_b = r.chrB
#             pos_b = r.posB - 1  # shift it down 1 for functions to get the correct range later on
#
#             # if the boundaries around the variant are on the same chromosome, then use the provided positions
#             # to identify overlaps
#             if pos_b == k:
#                 print(f"same chr: {pos_a} = {pos_b}")
#                 end_pos = chr_b
#                 start_pos = pos_a
#
#                 # find the annotations for the region
#                 a_repeats = find_reps_overlaps(k, start_pos, a_repeats, rep_regions, reps_df, posend=end_pos)
#                 a_f_sites = find_fs_overlaps(k, start_pos, a_f_sites, fs_regions, fs_df, posend=end_pos)
#                 a_genes, a_gene_lengths, a_gtd, a_gene_positions = find_gene_overlaps(k, start_pos, a_genes,
#                                                                                       genes_regions, genes_df,
#                                                                                       a_gene_lengths, a_gtd,
#                                                                                       a_gene_positions, posend=end_pos)
#                 # add empty entries to the b lists
#                 b_repeats.append("")
#                 b_f_sites.append("")
#                 b_genes.append("")
#                 b_gene_lengths.append("")
#                 b_gtd.append("")
#                 b_gene_positions.append("")
#
#             # else, use the boundaries only, by using the given position +/-1
#             else:
#                 start_pos = pos_b - 1
#                 end_pos = pos_a + 1
#
#                 a_repeats = find_reps_overlaps(k, pos_a, a_repeats, rep_regions, reps_df, posend=end_pos)
#                 a_f_sites = find_fs_overlaps(k, pos_a, a_f_sites, fs_regions, fs_df, posend=end_pos)
#                 a_genes, a_gene_lengths, a_gtd, a_gene_positions = find_gene_overlaps(k, start_pos, a_genes,
#                                                                                       genes_regions, genes_df,
#                                                                                       a_gene_lengths, a_gtd,
#                                                                                       a_gene_positions, posend=end_pos)
#
#                 b_repeats = find_reps_overlaps(chr_b, start_pos, b_repeats, rep_regions, reps_df, posend=pos_b)
#                 b_f_sites = find_fs_overlaps(chr_b, start_pos, b_f_sites, fs_regions, fs_df, posend=pos_b)
#                 b_genes, b_gene_lengths, b_gtd, b_gene_positions = find_gene_overlaps(chr_b, start_pos, b_genes,
#                                                                                       genes_regions,
#                                                                                       genes_df, b_gene_lengths, b_gtd,
#                                                                                       b_gene_positions, posend=pos_b)
#         v = add_annotation_column(v, a_repeats, "posA_repeats")
#         v = add_annotation_column(v, b_repeats, "posB_repeats")
#         v = add_annotation_column(v, a_f_sites, "posA_fragileSites")
#         v = add_annotation_column(v, b_f_sites, "posB_fragileSites")
#
#         v = add_annotation_column(v, a_genes, "posA_geneName")
#         v = add_annotation_column(v, b_genes, "posB_geneName")
#         v = add_annotation_column(v, a_gene_lengths, "posA_geneSize(bp)")
#         v = add_annotation_column(v, b_gene_lengths, "posB_geneSize(bp)")
#         v = add_annotation_column(v, a_gtd, "posA_gene-telomere_distances")
#         v = add_annotation_column(v, b_gtd, "posB_gene-telomere_distances")
#         v = add_annotation_column(v, a_gene_positions, "gene_positions")
#         v = add_annotation_column(v, b_gene_positions, "gene_positions")
#
#     final_df = pandas.concat(vars_by_chrom, ignore_index=True)
#
#     return final_df

def find_overlapping_features(regions, rep_regions, rep_df, gene_regions, gene_df, fs_regions, fs_df):
    """"
    annotate each region with the features that overlap with it
    # based on eccDNA pipeline.pipelines.annotate_dysgu_ouptut
    :param regions: dict of the df split by chromosome
    :param rep_regions: repeats regions NCLS
    :param rep_df: df of rep regions
    :param gene_regions: gene regions NCLS
    :param gene_df: df of gene regions
    :param fs_regions: fs regions NCLS
    :param fs_df: df of fs regions
    :return: annotations within the 'regions' dict
    """
    # for each chromosome in the input data
    for chrom, v in regions.items():
        # initialise the lists in which to record annotations
        repeats = []
        f_sites = []
        genes = []
        gene_lengths = []
        gtd = []
        gene_positions = []

        for i, r in v.iterrows():  # i = df index (int), r = df row
            # get the start and end positions for the region
            start_pos = int(r['start'])
            end_pos = int(r['end'])
            # print(start_pos, end_pos)
            repeats = find_reps_overlaps(chrom, start_pos, repeats, rep_regions, rep_df, posend=end_pos)
            f_sites = find_fs_overlaps(chrom, start_pos, f_sites, fs_regions, fs_df, posend=end_pos)
            genes, gene_lengths, gtd, gene_positions = find_gene_overlaps(chrom, start_pos, genes, gene_regions, gene_df, gene_lengths, gtd, gene_positions, posend=end_pos)

        v = add_annotation_column(v, repeats, "repeats")
        v = add_annotation_column(v, f_sites, "fragile_sites")
        v = add_annotation_column(v, genes, "genes")
        v = add_annotation_column(v, gene_lengths, "gene_lengths")
        v = add_annotation_column(v, gtd, "gene-telomere_distances")
        v = add_annotation_column(v, gene_positions, "gene_positions")

    final_df = pandas.concat(regions, ignore_index=True)

    return final_df


def annotate_just_fragile_sites(df):
    """
    annotate the regions of interest with overlapping fragile sites
    :param df: the input df
    :return: df annotated with fragile sites
    """
    # read the fragile sites reference file into a relevant df
    fs_df = get_fs_df()

    # create the NCLS for the fragile sites regions
    fs_regions = get_annotation_regions(fs_df)

    # create the NCLS for the regions of interest
    regions = get_regions_per_chr(df)

    for chrom, v in regions.items():
        # initialise the list of annotations
        f_sites = []
        for i, r in v.iterrows():  # i = df index (int), r = df row
            # get the start and end positions for the region
            start_pos = int(r['start'])
            end_pos = int(r['end'])
            f_sites = find_fs_overlaps(chrom, start_pos, f_sites, fs_regions, fs_df, posend=end_pos)

        v = add_annotation_column(v, f_sites, "fragile_sites")

    final_df = pandas.concat(regions, ignore_index=True)

    return final_df


def annotate_existing_genes(df, column_name):
    """
    if the input file already contains genes, then match them to the relevant reference genes
    in order to obtain positions and lengths
    :param df: input dataframe of positions etc.
    :param str column_name: name of the column containing the genes
    :return: df with added gene lengths, positions and distance to tels
    """
    # get the genes from the ref files
    genes = get_gene_annotations()

    # get the gene names from the input file
    df_genes = pandas.DataFrame(df, columns=[column_name])
    df_genes.rename(columns={column_name: "orig_gene"}, inplace=True)

    # convert any ensembl genes/create the alternative gene name column for downstream compatibility
    df_genes = convert_ensembl_ids(df_genes)

    # add the positions to the input genes
    pos_anns = add_position_annotations(df_genes, genes)

    # calculate the lengths and gene-telomere distances from the positions
    len_annotated = calculate_gene_lengths(pos_anns)
    # print(len_annotated.head())
    tel_annotated = calculate_distances_to_telomeres(len_annotated, isgene=True)

    # add the location format (dropped later if not needed, but adding it here makes other steps easier/more universal)
    full_df = add_gene_location(tel_annotated)

    trimmed_df = full_df.drop(columns=["alt_gene", "seqid", "gene_name", "gene", "Chr", "Start", "End", "Centromere"])
    trimmed_df.rename(columns={"start": "gene_start", "end": "gene_end"}, inplace=True)

    return trimmed_df


def get_annotation_regions(df):
    """
    get the NCLS for the regions provided in the annotations resources
    :param pandas.DataFrame df: annotation data df
    :return: dict of NCLS for the annotation type, for each chromosome
    """
    # initialise dict for region NCLS
    regions = {}

    # group the df by chromosome
    grouped_df = get_regions_per_chr(df)

    # add the regions to the regions dict
    for k, v in grouped_df.items():
        regions[k] = NCLS(v.start.values, v.end.values, v.index.values)

    return regions


def make_annotation_ncls(reps_df, genes_df, fs_df):
    """
    create the NCLS objects for annotations
    :param pandas.DataFrame reps_df: annotation data df
    :param pandas.DataFrame genes_df: annotation data df
    :param pandas.DataFrame fs_df: annotation data df
    :return: the NCLS objects for annotations
    """
    reps_regions = get_annotation_regions(reps_df)
    genes_regions = get_annotation_regions(genes_df)
    fs_regions = get_annotation_regions(fs_df)

    return reps_regions, genes_regions, fs_regions


def annotate_overlaps(df):
    """
    find all the overlaps between the regions of interest and features of interest
    :param df: df derived from input
    :return: same df containing annotations
    """
    reps_df = get_reps_df()
    genes_df = get_genes_df()
    fs_df = get_fs_df()

    # create the NCLSs for the repeats, genes and fragile sites
    reps_regions, genes_regions, fs_regions = make_annotation_ncls(reps_df, genes_df, fs_df)

    # split the regions to be annotated by chr
    regions = get_regions_per_chr(df)

    # find the overlaps between the features (NCLSs) and regions (list of dfs) to perform annotations
    anns = find_overlapping_features(regions, reps_regions, reps_df, genes_regions, genes_df, fs_regions, fs_df)

    return anns


def calculate_gene_lengths(df):
    """
    Calculate the lengths of the genes by end position - start position
    :param df: the dataframe containing genes and their positions
    :return: the dataframe with added gene lengths
    """
    df['gene_length'] = df['end'] - df['start']
    return df


def calculate_distances_to_telomeres(df, keep_int=False, isgene=False):
    """
    Calculate the distances between specified locations and the telomere on the same chromosome arm
    :param df: the dataframe containing the location info
    :param keep_int: whether to convert values to strings for printing, or maintain integers for further operations
    :param bool isgene: whether the type of region of interest is a gene
    :return: df containing telomere positions and distances between specified locations and telomere boundaries
    """
    if isgene:
        colname = "g-t_distance"
    else:
        colname = "reg-tel_distance"

    # temporarily remove na values to allow successful conversion for calculations
    no_na_df = df.dropna(subset=['seqid'])

    # add the required telomere and centromere positions
    no_na_df = add_tel_positions(no_na_df)
    no_na_df = add_cent_positions(no_na_df)

    # temporarily remove na values to allow successful conversion for calculations
    na_cents = no_na_df[no_na_df['Centromere'].isna()]
    no_na_df = no_na_df.dropna(subset=['Centromere'])

    # generate a df of the regions for which no annotations were found
    not_annotated = pandas.concat([df[df['seqid'].isna()], na_cents])

    # convert all positions to int type to allow relevant operations and removal of trailing .0
    if "gene_length" in no_na_df.columns:
        no_na_df = no_na_df.astype({"start": int, "end": int, "Start": int, "End": int,
                                    "gene_length": int, "Centromere": int})
    else:
        no_na_df = no_na_df.astype({"start": int, "end": int, "Start": int, "End": int, "Centromere": int})

    # if the region is on the p arm, calculate the distance to the p telomere. else, calculate to the q
    no_na_df.loc[(no_na_df["start"] < no_na_df["Centromere"]), colname] = no_na_df["start"] - no_na_df["Start"]
    no_na_df.loc[(no_na_df["start"] > no_na_df["Centromere"]), colname] = no_na_df["End"] - no_na_df["end"]

    # to get rid of the trailing .0s in this column, convert to int AND to string, to prevent the .0s returning
    # after re-addition of the temporarily removed entries
    no_na_df = no_na_df.astype({colname: int})
    no_na_df = no_na_df.astype({colname: str})
    # print(type(no_na_df["g-t_distance"][5]), type(no_na_df["start"][5]), type(no_na_df["end"][5]))
    # print(no_na_df.head())

    if not keep_int:
        # convert the integers to strings to prevent trailing .0s reappearing after NaN values are added back in

        if "gene_length" in no_na_df.columns:
            no_na_df = no_na_df.astype({"start": str, "end": str, "Start": str, "End": str, "gene_length": str,
                                        "Centromere": str, colname: str})
        else:
            no_na_df = no_na_df.astype({"start": str, "end": str, "Start": str, "End": str, "Centromere": str,
                                        colname: str})

        not_annotated = not_annotated.astype(str)

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
    if "gene_positions" in df.columns:
        df = df[['orig_gene', 'alt_gene', 'gene', 'gene_name', 'seqid', 'start', 'end',
                 'gene_length', "g-t_distance", "gene_positions", 'full_gene_loc']]
    else:
        df = df[['orig_gene', 'alt_gene', 'gene', 'gene_name', 'seqid', 'start', 'end',
                 'gene_length', "g-t_distance", "full_gene_loc", 'full_gene_loc']]

    # rename columns
    df.columns = ['Orig_gene', 'Alt_gene', 'Annotation_gene', 'Annotation_gene2', 'Chr', 'Start',
                  'End', 'Gene_length', "Distance_to_telomere", "Gene_positions", 'Gene_loc']

    return df


def output_full_csv(df, outname):
    """
    output the df as a single csv file
    :param df: df of genes and annotations
    :return: output file
    """
    # output the df to csv
    df.to_csv(os.path.join(outdir, outname), index=False)


def recreate_multisheet_excel_doc(df, outname):
    """
    organise and rename the columns of the annotated df, and output as a single csv file
    :param df: df of genes and annotations
    :return: output file
    """
    # read the file in
    in_genes = pandas.read_excel(os.path.join(refs_dir, input_file))

    with pandas.ExcelWriter(os.path.join(outdir, outname)) as writer:
        for _, s in in_genes.items():
            s = s.dropna()
            s_df = pandas.DataFrame(s.values, columns=['Orig_gene'])
            tmp_df = pandas.merge(s_df, df, how='left', on='Orig_gene')
            # print(tmp_df.columns)

            # reorder again
            tmp_df = tmp_df[['Orig_gene', 'Alt_gene', 'Annotation_gene', 'Annotation_gene2', 'Chr', 'Start', 'End',
                             'Gene_length', 'Distance_to_telomere', "Gene_positions", 'Gene_loc']]

            # print(tmp_df.head())
            tmp_df.to_excel(writer, sheet_name=s.name, index=False)


def recreate_manually_updated_excel_doc(df, outname):
    """
    organise and rename the columns of the annotated df, and output as a single csv file
    :param df: df of genes and annotations
    :return: output file
    """
    # read the file in
    in_genes = pandas.read_excel(os.path.join(refs_dir, input_file))
    # in_genes = pandas.read_excel('C:/Users/ec339/Downloads/ECKL_ALLGENES_SIZES_REPLETE_190124.xlsx', sheet_name=None)

    with pandas.ExcelWriter(os.path.join(outdir, outname)) as writer:
        for sname, s in in_genes.items():
            if not s.empty and not sname == "Mean gene lengths":
                # print(sname, s.columns)
                # if "Gene" in s.columns:
                #     s_df = s[["Gene"]]
                # else:
                #     s_df = s[["Orig_Gene"]]
                s_df = s[["Gene"]]
                s_df = s_df.dropna()
                s_df.columns = ["Orig_gene"]
                tmp_df = pandas.merge(s_df, df, how='left', on='Orig_gene')
                # print(tmp_df.columns)

                # reorder again
                tmp_df = tmp_df[['Orig_gene', 'Alt_gene', 'Annotation_gene', 'Annotation_gene2', 'Chr', 'Start', 'End', 'Gene_length', 'Distance_to_telomere', "Gene_positions", 'Gene_loc']]

                # print(tmp_df.head())
                tmp_df.to_excel(writer, sheet_name=sname, index=False)


def input_file_to_dataframe(infilepath=None):
    """
    Read in the input file and convert it into a df
    :param infilepath: path to input file
    :return: df containing input file gene names
    """
    if infilepath:
        in_genes = pandas.read_excel(infilepath)
    else:
        # read the file in, and add column names for ensembl genes
        in_genes = pandas.read_excel(os.path.join(refs_dir, input_file))

    # print(in_genes)

    # convert the excel file values into a single, unique list of gene names
    unique_genes = list(set([y.strip() for x in in_genes.values.tolist() for y in x if pandas.notna(y)]))

    # use the unique genes to initialise a dataframe for output, and strip spaces from the gene names just in case
    genes_df = pandas.DataFrame(unique_genes, columns=['orig_gene'])

    # no longer required: strip on y above does the same
    # genes_df = genes_df.map(lambda x: x.strip() if isinstance(x, str) else x)

    return genes_df


def split_multiple_genes(ann_df):
    """
    if there are multiple genes for any of the regions, split them into single rows
    :param ann_df: the dataframe containing all relevant annotations
    :return: edited df
    """
    to_explode = ['genes', 'gene_lengths', 'g-t_distance', 'gene_positions']
    for heading in to_explode:
        ann_df[heading] = ann_df[heading].str.split(', ')
    ann_df = ann_df.explode(to_explode)

    return ann_df





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

