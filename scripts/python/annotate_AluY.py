import pandas as pd
import numpy as np
import parsing_utils as parse
from parsing_utils import AnnotatedGFF as aGFF

AluY = aGFF().create(snakemake.input['vcf'], snakemake.input['gff'], header_len = 247)
AluY = parse.filter_motif(AluY, 'AluY')
AluY.to_csv(snakemake.output['tsv'], index = False, sep = '\t')