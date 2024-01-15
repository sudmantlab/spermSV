from parsing_utils import *

vcf = parse.sniffles_vcf(snakemake.input['vcf'], header_len = 249)
transform.gff3(vcf, save = True, filepath = snakemake.output['gff'])