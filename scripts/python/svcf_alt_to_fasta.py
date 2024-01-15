from parsing_utils import *

vcf = parse.sniffles_vcf(snakemake.input['vcf'], header_len = 249)
selected = vcf.query("SVTYPE not in ['BND', 'DEL', 'DUP'] & ALT not in ['N', '<INS>']")
write.alt_to_fasta(selected, snakemake.output['fa'])