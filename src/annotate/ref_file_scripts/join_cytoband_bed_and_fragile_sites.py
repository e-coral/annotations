import pandas

# read the fragile sites in as a df
fsites_df = pandas.read_csv('/home/emmon/OneDrive/eccDNA/dev/annotations/fragile_sites/fragile_site_and_cytoband.csv', sep='\t', header=0)

### don't split the values in this file: merge the values in the other file instead! ###
## split the chromosome and arm+position from 1p23 format cytobands
#chrsearch = '(?P<chr>.*?)[p-q].*'
#cytosearch = '[1-9X][0-9]?(?P<cyt>.*)'

## add the chromosome and cytoband columns
#fsites_df['chr'] = 'chr' + fsites_df.cytoband.str.extract(chrsearch, expand=True)
#fsites_df['name'] = fsites_df.cytoband.str.extract(cytosearch, expand=True)

# read the cytobands bed file in as a df
cytobands_df = pandas.read_csv('/home/emmon/OneDrive/eccDNA/dev/annotations/fragile_sites/chm13v2.0_cytobands_allchrs.bed', sep='\t', header=0)

# construct the cytoband name from the chromosome and the name column
cytobands_df['cytoband'] = cytobands_df.chrom.str.strip('chr') + cytobands_df.name


# merge the dfs on the cytoband
full_df = pandas.merge(fsites_df, cytobands_df, on='cytoband', how='left')

# view the df
full_df

# write to file
full_df.to_csv('/home/emmon/OneDrive/eccDNA/dev/annotations/fragile_sites/sites_and_positions.tsv', sep='\t', header=True)
