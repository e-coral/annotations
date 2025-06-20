# for annotating a list of genes in a single column in an excel file

import pandas
import glob
from pathlib import Path
import os
from src.annotate import annotate

outdir = (Path(__file__).parent / 'default_output').resolve()


for eachfile in glob.glob('C:/Users/ec339/Downloads/genes_to_annotate/*'):
    outname = f"annotated_{os.path.basename(eachfile)}"
    with pandas.ExcelWriter(os.path.join(outdir, outname)) as writer:
        indata = pandas.read_excel(eachfile, sheet_name=None, header=None)
        for sname, s in indata.items():
            if not s.empty:
                # print(type())
                unique_genes = list(set([y.strip() for y in s[0] if pandas.notna(y)]))
                genes_df = pandas.DataFrame(unique_genes, columns=['orig_gene'])
                converted_ensgs = annotate.convert_ensembl_ids(genes_df)
                genes = annotate.get_gene_annotations()
                pos_anns = annotate.add_position_annotations(converted_ensgs, genes)
                len_annotated = annotate.calculate_gene_lengths(pos_anns)
                tel_annotated = annotate.calculate_distances_to_centromeres_and_telomeres(len_annotated)
                df = annotate.add_gene_location(tel_annotated)
                df = annotate.reformat_for_output(df)
                print(df.head())

                df.to_excel(writer, sheet_name=sname, index=False)



# def get_all_anns():
#     """
#     Get all annotation types
#     :return: dfs containing relevant annotations and associated information
#     """
#     genes = annotate.get_gene_annotations()
#     tel_boundaries = annotate.get_telomere_boundaries()
#     centromeres = annotate.get_centromere_positions()
#
#     return genes, tel_boundaries, centromeres
#
#
#
# def main():
#     initial_df = annotate.input_file_to_dataframe()
#     # initial_df = get_genes_from_multiple_sheets()
#     # print(f"initial df: {len(initial_df)}")
#     converted_ensgs = annotate.convert_ensembl_ids(initial_df)
#     # print(converted_ensgs.head(10))
#     # converted_ensgs.to_csv(os.path.join(outdir, "converted_ensgs.csv"), index=False)
#     # genes, tels, cents = get_all_anns()
#     genes = annotate.get_gene_annotations()
#     # annotate.report_duplicated_genes(genes, os.path.join(outdir, "duplicated_genenames.csv"))
#     pos_anns = annotate.add_position_annotations(converted_ensgs, genes)
#     # print(f"final anns: {len(pos_anns)}")
#     # annotate.report_genes_with_no_annotations(pos_anns, os.path.join(outdir, "annotations_unavailable.csv"))
#     len_annotated = annotate.calculate_gene_lengths(pos_anns)
#     tel_annotated = annotate.calculate_distances_to_telomeres(len_annotated)
#     # print(tel_annotated.columns)
#     df = annotate.add_gene_location(tel_annotated)
#
#     df = annotate.reformat_for_output(df)
#
#     annotate.output_full_csv(df, "new_full_output.csv")
#
#     # annotate.recreate_multisheet_excel_doc(df, "new_sheet_output.xlsx")
#     annotate.recreate_manually_updated_excel_doc(df, "new_sheet_output.xlsx")
#     # cent_annotated = annotate.add_cent_positions(tel_annotated, cents)
#
#
# if __name__ == '__main__':
#     main()





