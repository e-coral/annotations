import os
import pandas
from pathlib import Path




# the telomere start and end positions aren't in a useful format (two rows for each chr)
# therefore, reformat to get the positions for the non-telomere section (i.e. first entry end and second start)

# initialise relevant paths
refs_dir = (Path(__file__).parent.parent / 'ref_files').resolve()
tels_file = 'originals/chm13v2.0_telomere.bed'

# read the telomere positions in as a df
ann_tels = pandas.read_csv(os.path.join(refs_dir, tels_file), sep='\t', header=None)


# add column names for easier processing
ann_tels.columns = ['Chr', 'Start', 'End']

# split the df into p and q telomeres, by whether they start at 0
p_tels = ann_tels[ann_tels['Start'] == 0]
q_tels = ann_tels[ann_tels['Start'] != 0]

# rename the columns for successful merging
p_tels.columns = ['Chr', 'telStart', 'Start']
q_tels.columns = ['Chr', 'End', 'telEnd']

# merge back into a single df with one row per chromosome
ann_tels = p_tels.merge(q_tels, on='Chr')

ann_tels = ann_tels[['Chr', 'Start', 'End']]
ann_tels.to_csv(os.path.join(refs_dir, 'non-telomere_regions.bed'), index=None, sep='\t')
