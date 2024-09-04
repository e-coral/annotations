import os
from pathlib import Path
from src.annotate import annotate
import pandas


outdir = (Path(__file__).parent / 'default_output').resolve()
outname = "fresh_annotate_tel_distances.xlsx"

with pandas.ExcelWriter(os.path.join(outdir, outname)) as writer:
    infile = pandas.read_excel('C:/Users/ec339/Downloads/col_edited_ECKL_ALLGENES_SIZES_REPLETE_190124.xlsx', sheet_name=None)
    for sname, s in infile.items():
        if not s.empty and not sname == "Mean gene lengths":
            # get them in the right order/drop the potentially incorrect lengths
            s_df = s[["Orig_Gene", "Alt_Gene", "Gene", "Chr", "Start", "End"]]
            # rename for annotation methods to find relevant values
            s_df.columns = ["Orig_gene", "Alt_gene", "Gene", "seqid", "start", "end"]

            # add gene lengths, telomere distances and gene locations
            s_df = annotate.calculate_gene_lengths(s_df)

            # print(s_df.columns)
            s_df = annotate.calculate_distances_to_telomeres(s_df)
            # print(s_df.columns)
            s_df = annotate.add_gene_location(s_df)

            # remove the irrelevant columns
            s_df = s_df[["Orig_gene", "Alt_gene", "Gene", "Chr", "start", "end", "gene_length", "g-t_distance", "full_gene_loc"]]

            # rename columns
            s_df.columns = ["Orig_Gene", "Alt_Gene", "Gene", "Chr", "Start", "End", "Gene_length", "Distance_to_telomere", "Gene_loc"]
            # print(s_df.head())

            # output to file
            s_df.to_excel(writer, sheet_name=sname, index=False)




