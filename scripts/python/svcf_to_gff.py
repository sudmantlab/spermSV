import parsing_utils as parse
from parsing_utils import SnifflesVCF as sVCF

vcf= sVCF().create(snakemake.input['vcf'], header_len = 247)
parse.as_gff3(vcf, save = True, filepath = snakemake.output['gff'])