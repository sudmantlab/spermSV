import pandas as pd
import numpy as np
import parsing_utils as parse
from parsing_utils import AnnotatedGFF as aGFF

AluY = aGFF().create(snakemake.input['vcf'], snakemake.input['gff'])
AluY = parse.filter_motif(AluY, 'AluY')
AluY.to_hdf(snakemake.output['h5'], index = False)