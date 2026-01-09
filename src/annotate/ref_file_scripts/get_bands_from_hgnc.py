import pandas

# read the bed file in as a df
liftover_df = pandas.read_csv('/home/emmon/OneDrive/eccDNA/dev/annotations/fragile_sites/hgnc_fragile_sites.txt', sep='\t', header=0)

# create the search to extract the chromosome number from the fragile site - the only place in the results file that contains chromosome info
bandsearch = '.*,.s*fra\s*\((?P<chr>.*?)\)\((?P<band>.*?)\)'

# create the band column, which extracts chr number and band position from the fragile site name
liftover_df['band'] = (liftover_df.Name.str.extract(bandsearch, expand=True))['chr'].astype(str) + (liftover_df.Name.str.extract(bandsearch, expand=True))['band'].astype(str)

# view the df
liftover_df

# write to file
liftover_df.to_csv('/home/emmon/OneDrive/eccDNA/dev/annotations/fragile_sites/hgnc_bandnames.txt', sep='\t', header=True)
