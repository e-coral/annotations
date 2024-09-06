# add other annotations to a file that already contains gene names

from pathlib import Path
from src.annotate import annotate


# initialise the output directory
out_dir = (Path(__file__).parent.parent / 'example_output_files').resolve()
in_dir = (Path(__file__).parent.parent / 'example_input_files').resolve()
# infile = f'{in_dir}/example_multisheet_genes_input.xlsx'
# outfile = f'{out_dir}/example_multisheet_genes_output.xlsx'
#
infile = f'{in_dir}/example_genes_input.csv'
outfile = f'{out_dir}/example_genes_output.csv'


def main():
    # extract all genes from the original file into a df
    initial_df = annotate.input_file_to_genes_df(infile)

    # convert any ensemble genes to gene names
    converted_ensgs = annotate.convert_ensembl_ids(initial_df)

    # annotate the genes with their positions, lengths and distances to telomeres
    genes = annotate.get_gene_annotations()
    pos_anns = annotate.add_position_annotations(converted_ensgs, genes)
    len_annotated = annotate.calculate_gene_lengths(pos_anns)
    tel_annotated = annotate.calculate_distances_to_telomeres(len_annotated)
    df = annotate.add_gene_location(tel_annotated)

    # reformat to the standard output
    df = annotate.reformat_for_output(df)

    # output standard format to csv
    annotate.output_full_csv(df, infile, output_dir=out_dir)

    # output original format to excel
    annotate.annotate_original_excel_doc(df, infile, output_dir=out_dir)


if __name__ == '__main__':
    main()
