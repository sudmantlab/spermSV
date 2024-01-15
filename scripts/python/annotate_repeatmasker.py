from parsing_utils import *

annotate.repeatmasker(snakemake.input['vcf'], snakemake.input['out'], header_len = 249, 
                      save = True, filepath = snakemake.output['tsv'])