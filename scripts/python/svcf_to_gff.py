import pandas as pd
import numpy as np
import parsing_utils as parse
from parsing_utils import SnifflesVCF as sVCF

sample = sVCF().create(snakemake.input['vcf'], header_len = 245)
parse.as_gff3(sample, save = True, filepath = snakemake.output['gff'])