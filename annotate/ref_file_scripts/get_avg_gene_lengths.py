import os
import pandas
from pathlib import Path


refs_dir = (Path(__file__).parent.parent / 'ref_files').resolve()
outdir = (Path(__file__).parent.parent.parent / 'default_output').resolve()

genes_file = 'chm13v2.0_RefSeq_Liftoff_genes_only.csv'
extras = 'CAT_liftoff_from_table_browser.csv'

ann_genes = pandas.read_csv(os.path.join(refs_dir, genes_file))
extra_anns = pandas.read_csv(os.path.join(refs_dir, extras))
# print(len(extra_anns))

ann_genes["Gene"] = ann_genes["attributes"].str.extract(r'.*;gene_name=(.*?);.*').fillna('')
# print(ann_genes.head(10))
# print(ann_genes.columns)
# print(set(ann_genes["type"].to_list()))
#
# ann_genes = ann_genes[["Gene", "seqid", "start", "end"]]
# ann_genes = ann_genes[["Gene", "start", "end"]]
# ann_genes["Length"] = ann_genes["end"] - ann_genes["start"]
# ann_genes.sort_values(by="Length", ascending=True, inplace=True)
# print(ann_genes.head(10))
# ann_genes.to_csv(os.path.join(outdir, "input1_length_ascending.csv"), index=False)
# ann_genes.sort_values(by="Length", ascending=False, inplace=True)
# print(ann_genes.head(10))

# print(ann_genes["Length"].mean())
# print(ann_genes["Length"].std())
# print(ann_genes["Length"].median())
# print(ann_genes["Length"].quantile([0.25, 0.75]))
# print(len(ann_genes["Gene"].unique()))
# print(len(ann_genes))


# extra_anns = extra_anns[["Gene", "start", "end"]]

# print(f"ann_genes is: {len(ann_genes)}")
# print(f"extra_anns is: {len(extra_anns)}")
# dups = extra_anns[extra_anns['name'].str.contains('AC016629.2-1')]
# dups["length"] = dups["end"] - dups["start"]
# print(f"AC016629.2-1:\n{dups}")

df = pandas.concat([ann_genes, extra_anns])
# print(len(df))
df.drop_duplicates(inplace=True)
# print(len(df))
# 6001 duplicates

df["length"] = df["end"] - df["start"]
print(df["length"].mean())
print(df["length"].std())
print(df["length"].median())
print(df["length"].quantile([0.25, 0.75]))
print(len(df))
print(len(df["Gene"].unique()))

# print(df.head(10))
