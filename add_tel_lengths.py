import os
from pathlib import Path
from annotate import annotate
import pandas


# initialise the output directory
outdir = (Path(__file__).parent / 'default_output').resolve()


def get_all_anns():
    """
    Get all annotation types
    :return: dfs containing relevant annotations and associated information
    """
    genes = annotate.get_gene_annotations()
    tel_boundaries = annotate.get_telomere_boundaries()
    centromeres = annotate.get_centromere_positions()

    return genes, tel_boundaries, centromeres


def get_genes_from_multiple_sheets():
    """
    For the edited workbook, get the names from the first column in every sheet
    :return: single-column df of all gene names
    """
    genes = []
    sheets = pandas.read_excel('C:/Users/ec339/Downloads/ECKL_ALLGENES_SIZES_REPLETE_190124.xlsx', sheet_name=None)
    for sheetname, df in sheets.items():
        if "Orig_Gene" in df.columns:
            genes.append(df['Orig_Gene'])
        elif "Gene" in df.columns:
            genes.append(df.Gene)
        else:
            print(f"no gene or orig gene column in {sheetname}")

    genes_list = pandas.concat(genes)
    # print(type(genes_list))
    # print(genes_list.head())

    # print(genes_list)

    unique_genes = list(set([y.strip() for y in genes_list if pandas.notna(y)]))
    # print(unique_genes)

    genes_df = pandas.DataFrame(unique_genes, columns=['orig_gene'])

    return genes_df


def main():
    initial_df = annotate.input_file_to_dataframe()
    # initial_df = get_genes_from_multiple_sheets()
    # print(f"initial df: {len(initial_df)}")
    converted_ensgs = annotate.convert_ensembl_ids(initial_df)
    # print(converted_ensgs.head(10))
    # converted_ensgs.to_csv(os.path.join(outdir, "converted_ensgs.csv"), index=False)
    # genes, tels, cents = get_all_anns()
    genes = annotate.get_gene_annotations()
    # annotate.report_duplicated_genes(genes, os.path.join(outdir, "duplicated_genenames.csv"))
    pos_anns = annotate.add_position_annotations(converted_ensgs, genes)
    # print(f"final anns: {len(pos_anns)}")
    # annotate.report_genes_with_no_annotations(pos_anns, os.path.join(outdir, "annotations_unavailable.csv"))
    len_annotated = annotate.calculate_gene_lengths(pos_anns)
    tel_annotated = annotate.calculate_distances_to_telomeres(len_annotated)
    # print(tel_annotated.columns)
    df = annotate.add_gene_location(tel_annotated)

    df = annotate.reformat_for_output(df)

    annotate.output_full_csv(df, "new_full_output.csv")

    # annotate.recreate_multisheet_excel_doc(df, "new_sheet_output.xlsx")
    annotate.recreate_manually_updated_excel_doc(df, "new_sheet_output.xlsx")
    # cent_annotated = annotate.add_cent_positions(tel_annotated, cents)


if __name__ == '__main__':
    main()
