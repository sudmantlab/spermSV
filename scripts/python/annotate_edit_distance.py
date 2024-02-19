from parsing_utils import *
from analysis_utils import *

# pull variants per specimen
# for right now, restrict to SINE/Alu, LINE/L1, Retroposon/SVA + only single motifs found
table = pd.read_table(snakemake.input['tsv']).query(f"origin == '{snakemake.wildcards['specimen']}' & family in ['SINE/Alu', 'LINE/L1', 'Retroposon/SVA']")
single = table[table['ID'].isin(table.value_counts('ID').loc[lambda x: x == 1].index)].reset_index()
bamfile = pysam.AlignmentFile(snakemake.input['bam'], 'rb')

variants = single[['ID', '#CHROM', 'POS', 'RNAMES']].to_records(index = False)
calculated = dict()
for record in variants:
    id, contig, position, alt_rnames = record
    calculated[id] = calc.edit_distance(contig, position, alt_rnames.split(','), bamfile)
distances = pd.DataFrame.from_dict(calculated).T

# distances['origin'] = distances.index.to_series().astype(str).str.split('_', expand = True)[0]
distances.index.name = "ID" # set index name, can be read in as a normal column later

write.table(distances, snakemake.output['tsv'], index = True) # retain index because it contains variant IDs