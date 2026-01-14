import os
import pandas
from pathlib import Path
from ncls import NCLS
import re

# initialise relevant paths
refs_dir = (Path(__file__).parent / 'ref_files').resolve()
outdir = (Path(__file__).parent.parent.parent / 'default_output').resolve()

genes_file = 'chm13v2.0_RefSeq_Liftoff_genes_only.csv'
extras = 'CAT_liftoff_from_table_browser.csv'
nontels_file = 'non-telomere_regions.bed'
centromeres_file = 'centromeres.csv'
ensembl_genes = 'ensembl_IDs_to_gene_names.csv'

# for region annotations, in addition to genes_file above
detailed_reps_file = 'chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed'
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

    # because there is a p and q entry for each chromosome, keep only the start for p and only the end of q
    pside = centromeres[centromeres['ChrArm'] == "p"]
    pside = pside[['chrom', 'chromStart', 'chromEnd']]

    qside = centromeres[centromeres['ChrArm'] == "q"]
    qside = qside[['chrom', 'chromEnd']]

    # merge to get chr, start and end of each centromere
    cen_df = pside.merge(qside, on='chrom')
    cen_df.columns = ["seqid", "censtart", "Centromere", "cenend"]

    # convert position to int type to allow relevant operations and removal of trailing .0
    cen_df = cen_df.astype({"censtart": int, "Centromere":int, "cenend": int})

    return cen_df


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
    genes_df = calculate_distances_to_centromeres_and_telomeres(genes_df, keep_int=True, isgene=True)

    return genes_df.astype({'gene_length': str, 'g-t_distance': str})


def get_fs_df():
    """
    read the fragile sites file into a relevant df
    :return: df of fragile site annotations
    """
    fs_df = pandas.read_csv(os.path.join(refs_dir, fs_file), sep='\t', header=0)
    return fs_df


def get_reps_df(repfile):
    """
    read the reps file into a relevant df
    :param repfile: source of the repeats to collect (either censat or rmsk)
    :return: df of reps annotations
    """
    reps_df = pandas.read_csv(os.path.join(refs_dir, repfile), sep='\t',
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
    Merge the current df of genes and positions with the centromere positions
    :param genes: df of genes and their positions
    :return: merged df, containing genes, gene positions, and centromere positions
    """
    cents = get_centromere_positions()
    cents_added = pandas.merge(genes, cents, on=['seqid'], how='left')

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

    # the repeats use title case Chr, so ensure the source matches this format
    chromkey = chrom.replace("Chr", "chr")
    # print(chromkey)

    try:
        # use the ncls find_overlap method to check whether the var positions overlap with repeat regions
        reps_overlaps = list(reps_regions[chromkey].find_overlap(pos, posend))
        # print(f"repeats found for {chrom}")
        # if there are any overlaps,
        # then for each overlap in the list, use the index of the region to identify the relevant annotation
        if reps_overlaps:
            for item in reps_overlaps:  # usually only 1 pos-annotation overlap, but it's possible for more
                # use the index to find the relevant annotation(s), and add them to the list for the function
                rep.append(reps_df.iloc[item[2]].repeats)

    except KeyError as err:
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
    if chrom.lower() not in ["chrm", "chry"]:
        # ensure lowercase chr is used, as this is how the fragile sites are recorded
        chromkey = chrom.replace("Chr", "chr")
        try:
            # use the ncls find_overlap method to check whether the var positions overlap with repeat regions
            fs_overlaps = list(fs_regions[chromkey].find_overlap(pos, posend))
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


def find_gene_overlaps(chrom, pos, gene_names, gene_regions, gene_df, gene_sizes, gtds, gcds, gene_positions, posend=None):
    """
    find positions of variant junctions that overlap with genes regions
    :param chrom: chromosome number of a variant
    :param pos: position of the junction of a variant
    :param list gene_names: list of gene name annotations to be added to
    :param gene_regions: NCLS of the gene regions
    :param gene_df: dataframe containing gene data for the chromosome
    :param list gene_sizes: list of gene size annotations to be added to
    :param list gtds: list of gene-telomere distances to be added to
    :param list gcds: list of gene-centromere distances to be added to
    :param list gene_positions: list of positions of genes to be added to
    :param posend: position of the end of the variant, if supplied (otherwise, set to be start pos + 1)
    :return: list of variant junction-overlapping gene sites regions for annotation
    """
    gtd = []
    gcd = []
    gene_length = []
    gene_name = []
    gene_pos = []

    if posend is None:
        posend = pos + 1

    chromkey = chrom.replace("Chr", "chr")

    try:
        genes_overlaps = list(gene_regions[chromkey].find_overlap(pos, posend))
        if genes_overlaps:
            for item in genes_overlaps:
                relevant_data = gene_df.iloc[item[2]]

                # if the type of the annotation is a gene, then get the gene name and size
                if relevant_data.type == "gene":
                    gene_length.append(relevant_data.gene_length)
                    gtd.append(relevant_data["g-t_distance"])
                    gcd.append(relevant_data["g-c_distance"])
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
        print(f"{err}\nNo gene annotations available for {chrom}:{pos}-{posend}.")

    gene_sizes.append(", ".join(gene_length))
    gene_names.append(", ".join(gene_name))
    gtds.append(", ".join(gtd))
    try:
        gcds.append(", ".join(gcd))
    except TypeError:
        # hack to deal with chromosomes that don't have centromeres
        gcd_strings = [str(i) for i in gcd]
        gcds.append(", ".join(gcd_strings))

    gene_positions.append(", ".join(gene_pos))

    return gene_names, gene_sizes, gtds, gcds, gene_positions


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


def format_output_columns(df, orig_col_names, rmsk=False):
    """
    For the eventual output, drop irrelevant columns from the dataframe, and rename others to original names
    :param df: dataframe of input values and annotations
    :param orig_col_names: list of the original column names
    :return: formatted dataframe of values and annotations
    """
    if orig_col_names:
        if rmsk:
        # add the names of the new columns to the names of the original columns
            for col in ['tChr', 'tStart', 'tEnd', 'censtart', 'Centromere', 'cenend', 'reg-tel_distance', 'reg-cen_distance',
                        'censat_repeats', 'rmsk_repeats', 'fragile_sites', 'genes', 'gene_lengths', 'gene-telomere_distances',
                        'gene_positions']:
                orig_col_names.append(col)
        else:
            for col in ['tChr', 'tStart', 'tEnd', 'censtart', 'Centromere', 'cenend', 'reg-tel_distance', 'reg-cen_distance',
                        'censat_repeats', 'fragile_sites', 'genes', 'gene_lengths', 'gene-telomere_distances',
                        'gene_positions']:
                orig_col_names.append(col)
        # rename the columns accordingly
        df.columns = orig_col_names
        # drop the telomere and centromere columns that aren't needed in the final output
        df = df.drop(columns=['tChr', 'tStart', 'tEnd', 'censtart', 'Centromere', 'cenend'])

    else:
        df = df.drop(columns=['Chr', 'Start', 'End', 'censtart', 'Centromere', 'cenend'])

    # remove nans and set values to strings to prevent trailing .0s
    df.fillna('', inplace=True)
    df = df.astype(str)

    return df


def format_just_distance_output_columns(df, orig_col_names):
    """
    For the eventual output, drop irrelevant columns from the dataframe, and rename others to original names
    :param df: dataframe of input values and annotations
    :param orig_col_names: list of the original column names
    :return: formatted dataframe of values and annotations
    """
    # add the names of the new columns to the names of the original columns
    # for col in ['gene chrom', 'gene start', 'gene end', 'tChr', 'tStart', 'tEnd', 'Centromere',
    #             'region-telomere_distance', 'repeats', 'fragile_sites',
    #             'genes', 'gene_lengths', 'g-t_distance', 'gene_positions']:
    for col in ['tChr', 'tStart', 'tEnd', 'censtart', 'Centromere', 'cenend', 'reg-tel_distance', 'reg-cen_distance']:
        orig_col_names.append(col)

    # rename the columns accordingly
    df.columns = orig_col_names

    # drop the telomere and centromere columns that aren't needed in the final output
    df = df.drop(columns=['tChr', 'tStart', 'tEnd', 'censtart', 'Centromere', 'cenend'])

    # remove nans and set values to strings to prevent trailing .0s
    df.fillna('', inplace=True)
    df = df.astype(str)

    return df


def find_overlapping_features(regions, rep_regions, rep_df, d_rep_regions, d_rep_df, gene_regions, gene_df, fs_regions, fs_df):
    """"
    annotate each region with the features that overlap with it
    # based on eccDNA pipeline.pipelines.annotate_dysgu_ouptut
    :param regions: dict of the df split by chromosome
    :param rep_regions: repeats regions NCLS
    :param rep_df: df of rep regions
    :param d_rep_regions: repeatmasker rep regions NCLS
    :param d_rep_df: df of repeatmasker rep regions
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
        d_repeats = []
        f_sites = []
        genes = []
        gene_lengths = []
        gtd = []
        gcd = []
        gene_positions = []

        for i, r in v.iterrows():  # i = df index (int), r = df row
            # get the start and end positions for the region
            start_pos = int(r['start'])
            end_pos = int(r['end'])
            repeats = find_reps_overlaps(chrom, start_pos, repeats, rep_regions, rep_df, posend=end_pos)
            if d_rep_regions:
                d_repeats = find_reps_overlaps(chrom, start_pos, d_repeats, d_rep_regions, d_rep_df, posend=end_pos)
            f_sites = find_fs_overlaps(chrom, start_pos, f_sites, fs_regions, fs_df, posend=end_pos)
            genes, gene_lengths, gtd, gcd, gene_positions = find_gene_overlaps(chrom, start_pos, genes, gene_regions,
                                                                               gene_df, gene_lengths, gtd, gcd,
                                                                               gene_positions, posend=end_pos)

        v = add_annotation_column(v, repeats, "censat_repeats")
        if d_rep_regions:
            v = add_annotation_column(v, d_repeats, "rmsk_repeats")
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
    tel_annotated = calculate_distances_to_centromeres_and_telomeres(len_annotated, isgene=True)

    # add the location format (dropped later if not needed, but adding it here makes other steps easier/more universal)
    full_df = add_gene_location(tel_annotated)

    trimmed_df = full_df.drop(columns=["alt_gene", "seqid", "gene_name", "gene", "Chr", "Start", "End", "censtart", "Centromere", "cenend"])
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


def make_annotation_ncls(reps_df, detailed_reps_df, genes_df, fs_df):
    """
    create the NCLS objects for annotations
    :param pandas.DataFrame reps_df: censat annotation data df
    :param pandas.DataFrame detailed_reps_df: repeatmasker annotation data df
    :param pandas.DataFrame genes_df: genes annotation data df
    :param pandas.DataFrame fs_df: fragile site annotation data df
    :return: the NCLS objects for annotations
    """
    reps_regions = get_annotation_regions(reps_df)
    detailed_reps_regions = get_annotation_regions(detailed_reps_df)
    genes_regions = get_annotation_regions(genes_df)
    fs_regions = get_annotation_regions(fs_df)

    return reps_regions, detailed_reps_regions, genes_regions, fs_regions


def annotate_overlaps(df, rmsk=False):
    """
    find all the overlaps between the regions of interest and features of interest
    :param df: df derived from input
    :param bool rmsk: whether to use rmsk or not
    :return: same df containing annotations
    """
    reps_df = get_reps_df(reps_file)
    genes_df = get_genes_df()
    fs_df = get_fs_df()

    # create the NCLSs for the repeats, genes and fragile sites
    reps_regions = get_annotation_regions(reps_df)
    genes_regions = get_annotation_regions(genes_df)
    fs_regions = get_annotation_regions(fs_df)

    # split the regions to be annotated by chr
    regions = get_regions_per_chr(df)

    # create the NCLS for rmsk repeats, if requested
    if rmsk:
        d_reps_df = get_reps_df(detailed_reps_file)
        d_reps_regions = get_annotation_regions(d_reps_df)
    else:
        d_reps_df = None
        d_reps_regions = None

    # find the overlaps between the features (NCLSs) and regions (list of dfs) to perform annotations
    anns = find_overlapping_features(regions, reps_regions, reps_df, d_reps_regions, d_reps_df, genes_regions,
                                     genes_df, fs_regions, fs_df)

    return anns


def calculate_gene_lengths(df):
    """
    Calculate the lengths of the genes by end position - start position
    :param df: the dataframe containing genes and their positions
    :return: the dataframe with added gene lengths
    """
    df['gene_length'] = df['end'] - df['start']
    return df

def annotate_dataframe_output_excel(df, outfile, colname='chrom', explode=False):
    """
    for a supplied dataframe, annotate and return as excel
    :param df: a dataframe with minimum columns chrom, start, end
    :param outfile: full path to output file
    :param colname: name of the chromosome column, if not chrom
    :param explode: whether to explode data out to multiple rows (default False)
    :return: excel file with annotated data
    """
    # store the original column names to be able to add them back in later
    orig_columns = [col.strip() for col in df.columns]

    # rename the seqnames column, because some of the methods may be based on 'seqid'
    df.rename(columns={colname: "seqid"}, inplace=True)

    # calculate the distance between the region and the nearest telomere
    res = calculate_distances_to_centromeres_and_telomeres(df)

    # annotate the genes, fragile sites, repeats and gene sizes
    regions_df = annotate_overlaps(res)

    # reformat the final df to match the input
    final = format_output_columns(regions_df, orig_columns)

    # split out the genes into separate columns, if preferred
    if explode:
        final = split_multiple_genes(final)

    # create the output file
    with pandas.ExcelWriter(f"{outfile}.xlsx") as writer:
        final.to_excel(writer, sheet_name='sheet1', index=False)


def annotate_standard_csv_input_file(infile, outfile, colname="chrom", explode=False):
    """
    read in an input file which contains at least chr, start and end columns for regions of interest,
    annotate it, and write an output file containing annotated input regions
    :param infile: path to input file
    :param outfile: path to output file
    :param colname: name of the column containing chromosomes, if not chrom
    :param explode: whether to 'explode' the gene annotations to one per row
    :return: annotated file
    """
    # read in the data
    res = pandas.read_csv(infile, header=0)

    # store the original column names to be able to add them back in later
    orig_columns = [col.strip() for col in res.columns]

    # rename the seqnames column, because some of the methods may be based on 'seqid'
    res.rename(columns={colname: "seqid"}, inplace=True)

    # calculate the distance between the region and the nearest telomere
    res = calculate_distances_to_centromeres_and_telomeres(res)

    # annotate the genes, fragile sites, repeats and gene sizes
    regions_df = annotate_overlaps(res)

    # reformat the final df to match the input
    final = format_output_columns(regions_df, orig_columns)

    # ensure not '.csv.csv'
    if outfile.endswith(".csv"):
        extension = ""
    elif outfile.endswith(".tsv"):
        extension = ""
    else:
        extension = ".csv"

    # split out the genes into separate columns, if preferred
    if explode:
        final = split_multiple_genes(final)
        # write output file
        final.to_csv(os.path.join(outdir, f"{outfile}-exploded{extension}"), index=False)
    else:
        final.to_csv(os.path.join(outdir, f"{outfile}{extension}"), index=False)


def annotate_standard_excel_input_file(infile, outfile, colname="chrom", startname="start", endname="end",
                                       explode=False, make_csv=False, headerin=True, headerout=True):
    """
    read in a standard input file, which contains chr, start and end of regions of interest at minimum,
    annotate it, and write the annotations to an output file

    :param infile: path to input file
    :param outfile: path to output file
    :param colname: name of the column containing the chromosomes, if not chrom
    :param startname: name of the start column
    :param endname: name of the end column
    :param explode: whether to 'explode' the gene annotations to one per row
                    (if false, all gene annotations will be written to one cell of the output file)
    :param make_csv: whether to also write a csv file
    :param headerin: whether the input file has a header (default True)
    :param headerout: whether the output file needs a header (default True)
    :return: annotated file
    """
    if explode:
        ext = "-exploded"
    else:
        ext = ""

    header_row = 0 if headerin else None

    # create the output file
    with pandas.ExcelWriter(f"{outfile}{ext}.xlsx") as writer:
        # read in the input file
        infile = pandas.read_excel(infile, sheet_name=None, header=header_row)
        # for each sheet in the input file (probably only one for each in this case)
        for sname, s in infile.items():
            # if the sheet is not empty
            if not s.empty:
                if headerin:
                    # store the original column names to be able to add them back in later
                    orig_columns = [col.strip() for col in s.columns]
                    # standardise the names of the columns
                    s.rename(columns={colname: "seqid"}, inplace=True)
                    s.rename(columns={startname: "start"}, inplace=True)
                    s.rename(columns={endname: "end"}, inplace=True)
                else:
                    # rename the default index values to be the standard
                    s.rename(columns={0: "seqid"}, inplace=True)
                    s.rename(columns={1: "start"}, inplace=True)
                    s.rename(columns={2: "end"}, inplace=True)
                    s.rename(columns={3: "signal"}, inplace=True) # quite specific - may need commenting out

                    if headerout:
                        # flag no orig columns
                        orig_columns = s.columns.tolist()
                    else:
                        orig_columns = None

                # calculate the distance between the region and the nearest telomere
                s = calculate_distances_to_centromeres_and_telomeres(s)

                # annotate the genes, fragile sites, repeats and gene sizes
                regions_df = annotate_overlaps(s)

                # reformat the final df to match the input
                final = format_output_columns(regions_df, orig_columns)

                # split out the genes into separate columns, if preferred
                if explode:
                    final = split_multiple_genes(final)

                # create output files
                final.to_excel(writer, sheet_name=sname, index=False)

                if make_csv:
                    final.to_csv(f"{outfile}_{sname}{ext}.csv", index=False)


def annotate_gene_location_data(infile, outfile, colname, rmsk=False):
    """
    read in a file containing chr:start-end formatted data, annotate it, and write the annotations to an output file
    :param infile: full path to input file
    :param outfile: full path to output file
    :param str colname: name of column containing the relevant data
    :param rmsk: whether to add rmsk annotations
    :return:
    """
    # create the output file
    with pandas.ExcelWriter(outfile) as writer:
        # read in the input file
        infile = pandas.read_excel(infile, sheet_name=None)
        # for each sheet in the input file (probably only one for each in this case)
        for sname, s in infile.items():
            # if the sheet is not empty
            if not s.empty:
                # convert the positions into separate columns for compatibility with the rest of the pipeline
                res = get_separate_chr_pos(s, colname)
                # res.to_csv(os.path.join(outdir, csvoutname), sep='\t', index=False)

                # store the original column names, plus the column names for new chr pos columns above
                orig_columns = [col.strip() for col in res.columns]

                # calculate the distance between the region and the nearest telomere
                res = calculate_distances_to_centromeres_and_telomeres(res)

                # annotate the genes, fragile sites, repeats and gene sizes
                regions_df = annotate_overlaps(res, rmsk)

                # reformat the final df to match the input
                final = format_output_columns(regions_df, orig_columns)

                # write output file
                final.to_excel(writer, sheet_name=sname, index=False)


def calculate_distances_to_centromeres_and_telomeres(df, keep_int=False, isgene=False):
    """
    Calculate the distances between specified locations and the telomere on the same chromosome arm
    :param df: the dataframe containing the location info
    :param keep_int: whether to convert values to strings for printing, or maintain integers for further operations
    :param bool isgene: whether the type of region of interest is a gene
    :return: df containing telomere positions and distances between specified locations and telomere boundaries
    """
    if isgene:
        colname = "g-t_distance"
        cendist_name = 'g-c_distance'
    else:
        colname = "reg-tel_distance"
        cendist_name = 'reg-cen_distance'

    # temporarily remove na values to allow successful conversion for calculations
    no_na_df = df.dropna(subset=['seqid'])

    # add the required telomere and centromere positions
    no_na_df = add_tel_positions(no_na_df)
    no_na_df = add_cent_positions(no_na_df)

    # temporarily remove na values to allow successful conversion for calculations
    na_centromeres = no_na_df[no_na_df['Centromere'].isna()]
    no_na_df = no_na_df.dropna(subset=['Centromere'])

    # generate a df of the regions for which no annotations were found
    not_annotated = pandas.concat([df[df['seqid'].isna()], na_centromeres])

    # convert all positions to int type to allow relevant operations and removal of trailing .0
    if "gene_length" in no_na_df.columns:
        no_na_df = no_na_df.astype({"start": int, "end": int, "Start": int, "End": int,
                                    "gene_length": int, "censtart": int, "Centromere": int, "cenend": int})
    else:
        no_na_df = no_na_df.astype({"start": int, "end": int, "Start": int, "End": int,"censtart": int,
                                    "Centromere": int, "cenend": int})

    # if the region is on the p arm, calculate the distance to the p telomere. else, calculate to the q
    no_na_df.loc[(no_na_df["start"] < no_na_df["Centromere"]), colname] = no_na_df["start"] - no_na_df["Start"]
    no_na_df.loc[(no_na_df["start"] < no_na_df["Centromere"]), cendist_name] = no_na_df["censtart"] - no_na_df["start"]
    no_na_df.loc[(no_na_df["start"] > no_na_df["Centromere"]), colname] = no_na_df["End"] - no_na_df["end"]
    no_na_df.loc[(no_na_df["start"] > no_na_df["Centromere"]), cendist_name] = no_na_df["end"] - no_na_df["cenend"]

    # to get rid of the trailing .0s in the new columns, convert to int AND to string, to prevent the .0s returning
    # after re-addition of the temporarily removed entries
    no_na_df = no_na_df.astype({colname: int})
    no_na_df = no_na_df.astype({colname: str})
    no_na_df = no_na_df.astype({cendist_name: int})
    no_na_df = no_na_df.astype({cendist_name: str})
    # print(type(no_na_df["g-t_distance"][5]), type(no_na_df["start"][5]), type(no_na_df["end"][5]))
    # print(no_na_df.head())

    if not keep_int:
        # convert the integers to strings to prevent trailing .0s reappearing after NaN values are added back in

        if "gene_length" in no_na_df.columns:
            no_na_df = no_na_df.astype({"start": str, "end": str, "Start": str, "End": str, "gene_length": str,
                                        "censtart": str, "Centromere": str, "cenend": str, colname: str, cendist_name: str})
        else:
            no_na_df = no_na_df.astype({"start": str, "end": str, "Start": str, "End": str, "censtart": str,
                                        "Centromere": str, "cenend": str, colname: str, cendist_name: str})

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
    # print(df.columns)
    # keep only relevant info
    if "gene_positions" in df.columns:
        df = df[['orig_gene', 'alt_gene', 'gene', 'gene_name', 'seqid', 'start', 'end',
                 'gene_length', "reg-tel_distance", "reg-cen_distance", "gene_positions", 'full_gene_loc']]
    else:
        df = df[['orig_gene', 'alt_gene', 'gene', 'gene_name', 'seqid', 'start', 'end',
                 'gene_length', "reg-tel_distance", "reg-cen_distance", "full_gene_loc", 'full_gene_loc']]

    # rename columns
    df.columns = ['Orig_gene', 'Alt_gene', 'Annotation_gene', 'Annotation_gene2', 'Chr', 'Start',
                  'End', 'Gene_length', "Distance_to_telomere", "Distance_to_centromere", "Gene_positions", 'Gene_loc']

    return df


def output_full_csv(df, infile, outname=None, output_dir=outdir):
    """
    output the df as a single csv file
    :param infile: path to input file
    :param df: df of genes and annotations
    :param output_dir: output directory (default, or overwritten by function call)
    :param outname: output file name
    :return: output file
    """
    # make the full outfile path
    if outname:
        outfile = os.path.join(outdir, outname)
    else:
        outfile = os.path.join(output_dir, f'annotated-{Path(infile).stem}.csv')

    # output the df to csv
    df.to_csv(outfile, index=False)


def annotate_original_excel_doc(df, infile, outname=None, output_dir=outdir):
    """
    merge the new annotation data with the existing data in the input file, and output
    :param output_dir: path to the output directory (default, or overwritten by function call)
    :param df: df of genes and annotations
    :param infile: input filepath
    :param outname: output file name
    :return: output file
    """
    # read the file in, and check it's an excel doc
    try:
        # read in the file as excel
        in_data = pandas.read_excel(infile, sheet_name=None)

    except ValueError:
        print("Cannot output an annotated version of original excel file, because the original file was in csv format")
        exit()

    # make the full outfile path
    if outname:
        outfile = os.path.join(output_dir, outname)
    else:
        outfile = os.path.join(output_dir, f'annotated-{Path(infile).stem}.xlsx')

    # create and write to the output file
    with pandas.ExcelWriter(outfile) as writer:
        # for each sheet
        for sname, s in in_data.items():
            if not s.empty:
                # drop blank cells to allow for merging
                try:
                    s_df = s[['Gene']]

                except ValueError:
                    s_df = s[['Orig_Gene']]

                # drop blank cells and rename the gene column, for easier merging
                s_df = s_df.dropna()
                s_df.columns = ["Orig_gene"]
                # merge the new (annotation) data with the original data
                tmp_df = pandas.merge(s_df, df, how='left', on='Orig_gene')

                # reorder again
                tmp_df = tmp_df[['Orig_gene', 'Alt_gene', 'Annotation_gene', 'Annotation_gene2', 'Chr', 'Start', 'End',
                                 'Gene_length', 'Distance_to_telomere', "Distance_to_centromere", "Gene_positions", 'Gene_loc']]

                # write the sheet to the workbook
                tmp_df.to_excel(writer, sheet_name=sname, index=False)
            else:
                blank = pandas.DataFrame()
                blank.to_excel(writer, sheet_name=sname, index=False)


def find_genes_in_df(df, genes, sheetname="input file"):
    """
    extract the genes from the input file
    :param df: automatic dataframe from reading in csv/excel
    :param genes: list of gene names (empty at first)
    :param sheetname: string name of sheet, if input has multiple sheets
    :return: list of genes from the data
    """
    if "Gene" in df.columns:
        genes.append(df['Gene'])
    elif "gene" in df.columns:
        genes.append(df['gene'])
    elif "Orig_Gene" in df.columns:
        genes.append(df['Orig_Gene'])
    elif "Orig_gene" in df.columns:
        genes.append(df['Orig_gene'])
    elif "orig_gene" in df.columns:
        genes.append(df['orig_gene'])
    elif "Alt_Gene" in df.columns:
        genes.append(df['Alt_Gene'])
    elif "Alt_gene" in df.columns:
        genes.append(df['Alt_gene'])
    elif "alt_gene" in df.columns:
        genes.append(df['alt_gene'])
    else:
        print(f"Warning: no relevant genes column identified in {sheetname}. "
              f"The sheet will be ignored. Please check whether this is correct behaviour.")

    return genes


def input_file_to_genes_df(infilepath):
    """
    Read in the input file and convert it into a df containing only unique input genes
    :param infilepath: path to input file
    :return: df containing gene names from the input file
    """
    genes = []

    try:
        # read in the file as excel
        in_data = pandas.read_excel(infilepath, sheet_name=None)
        # for each sheet, extract the genes from a relevant genes column, in order of preference
        for sheetname, df in in_data.items():
            if not df.empty:
                genes = find_genes_in_df(df, genes, sheetname)

    except ValueError:
        in_data = pandas.read_csv(infilepath)
        genes = find_genes_in_df(in_data, genes)

    # if genes is a list of lists, then concatenate it. Else, skip.
    try:
        genes = pandas.concat(genes)
    except ValueError:
        pass

    # make the single list of genes unique
    unique_genes = list(set([y.strip() for y in genes if pandas.notna(y)]))

    # use the unique genes to initialise a dataframe for output, and strip spaces from the gene names just in case
    genes_df = pandas.DataFrame(unique_genes, columns=['orig_gene'])

    return genes_df


def split_multiple_genes(ann_df):
    """
    if there are multiple genes for any of the regions, split them into single rows
    :param ann_df: the dataframe containing all relevant annotations
    :return: edited df
    """
    to_explode = ['genes', 'gene_lengths', 'gene-telomere_distances', 'gene_positions']
    for heading in to_explode:
        ann_df[heading] = ann_df[heading].str.split(', ')
    ann_df = ann_df.explode(to_explode)

    return ann_df

