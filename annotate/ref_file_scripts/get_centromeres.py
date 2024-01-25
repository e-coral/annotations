import os
import pandas
from pathlib import Path


refs_dir = (Path(__file__).parent.parent / 'ref_files').resolve()

cytobands = pandas.read_csv(os.path.join(refs_dir, 'originals/cytobands.tsv'), sep='\t', index_col=0)

centromeres = cytobands[cytobands['gieStain'] == 'acen']
centromeres['ChrArm'] = centromeres.name.str[:1]
print(centromeres.head(10))

centromeres.to_csv(os.path.join(refs_dir, 'centromeres.csv'))
