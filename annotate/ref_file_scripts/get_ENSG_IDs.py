import csv
import os
import biomart
from pathlib import Path

# initialise relevant path
outdir = (Path(__file__).parent.parent / 'ref_files').resolve()


def get_ensembl_mappings():
    """
    use biomart to get gene names from ensembl IDs, adapted from
    https://gist.github.com/ben-heil/cffbebf8865795fe2efbbfec041da969
    :return: dict of ensembl ID and gene symbol
    """
    server = biomart.BiomartServer('http://ensembl.org/biomart')
    mart = server.datasets['hsapiens_gene_ensembl']

    # Get the mapping between the attributes
    # biomart doesn't yet have t2t, so positions would be for hg38. Therefore, only get names
    # response = mart.search({'attributes': ['hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 'start_position',
    #                                        'end_position']})
    response = mart.search({'attributes': ['hgnc_symbol', 'ensembl_gene_id']})
    data = response.raw.data.decode('ascii')

    # initialise dict
    ensembl_to_genesymbol = {}

    # collect data and add to the dict
    for line in data.splitlines():
        line = line.split('\t')
        gene_symbol = line[0]
        ensembl_gene = line[1]
        ensembl_to_genesymbol[ensembl_gene] = gene_symbol

    return ensembl_to_genesymbol


# get the dict of ensembl IDs and gene names
ens_dict = get_ensembl_mappings()

# print(ens_dict)

# print the dict to a csv to circumvent issues related to connection errors
# ens_df = pandas.DataFrame.from_dict(ens_dict)
#
# print(ens_df.head())
#
# exit()
# ens_df.to_csv(os.path.join(outdir, "ensembl_IDs.csv"), index=False)

with open(os.path.join(outdir, "test_ensembl_IDs_to_gene_names.csv"), 'w') as f:
    writer = csv.writer(f)
    for g in ens_dict.items():
        writer.writerow(g)
