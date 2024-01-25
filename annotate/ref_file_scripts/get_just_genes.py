import os
import pandas
from pathlib import Path


annotations_dir = (Path(__file__).parent.parent / 'ref_files').resolve()

# import the full data
full_genes_df = pandas.read_csv(os.path.join(annotations_dir, 'originals/chm13v2.0_RefSeq_Liftoff_v4.gff3'), sep='\t', comment='#', names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

# subset to only gene entries
genes_df = full_genes_df[full_genes_df.type.isin(["gene"])]

# write out the new file
genes_df.to_csv(os.path.join(annotations_dir, 'chm13v2.0_RefSeq_Liftoff_genes_only.csv'))
# full_genes_df.to_csv(os.path.join(output_dir, 'chm13v2.0_RefSeq_Liftoff_full.csv'))
