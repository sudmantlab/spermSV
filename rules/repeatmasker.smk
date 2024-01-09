# rule alt_to_fasta:
#     # Takes a (sniffles) vcf and writes INS alleles to a fasta file for RepeatMasker. 
#     input: 
#         vcf='output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz'
#     output: 
#         fa = 'analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa'
#     wildcard_constraints:
#         specimen = '[A-Za-z0-9]+'
#     threads:
#         10
#     run:
#         def write_fasta(vcf, outfile):
#             # Takes the alt output of Sniffles2 and writes a fasta file for input through RepeatMasker.
#             table = vcf[['ID', 'ALT']]
#             with open(outfile, 'w') as f:
#                 for header, seq in table.to_records(index=False):
#                     f.write(f'>{header}\n{seq}\n')
        
#         vcf = pd.read_table(input.vcf, skiprows = 241) # skip vcf header
#         print("Writing ALT sequences to ", output.fa)
#         write_fasta(vcf, output.fa)
        
# rule annotate_alts:
#     # Takes a fasta file containing ALT allele sequences from Sniffles2 output, then analyzes the sequences using RepeatMasker.
#     # Relevant outputs: a .tbl (difficult to parse) and a .gff (easier to parse) denoting identified repeat motifs.
#     input:
#         "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa"
#     output:
#         "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa.tbl",
#         "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa.out.gff",
#         "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa.masked",
#         "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa.cat"
#     wildcard_constraints:
#         specimen = '[A-Za-z0-9]+'
#     params:
#         outdir = "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker",
#         species = config['repeatmasker']['species'],
#         engine = config['repeatmasker']['engine']
#     log:
#         "logs/analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.log"
#     threads: 10
#     conda:
#         '../envs/RepeatMasker.yml'
#     shell:
#         """
#         RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff \
#         -species {params.species} -dir {params.outdir} {input} &> {log}
#         """

# use rule alt_to_fasta as alt_to_fasta_duplomap with:
#     input:
#         vcf = 'output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.filt.vcf.gz'
#     output:
#         fa = 'analysis/{refalias}/{analysis}/files/{mapper}/duplomap/repeatmasker/{specimen}.filt.alt.fa'

# use rule annotate_alts as annotate_alts_duplomap with:
#     input:
#         'analysis/{refalias}/{analysis}/files/{mapper}_duplomap/repeatmasker/{specimen}.filt.alt.fa'
#     output:
#         "analysis/{refalias}/{analysis}/files/{mapper}_duplomap/repeatmasker/{specimen}.filt.alt.fa.tbl",
#         "analysis/{refalias}/{analysis}/files/{mapper}_duplomap/repeatmasker/{specimen}.filt.alt.fa.out.gff",
#         "analysis/{refalias}/{analysis}/files/{mapper}_duplomap/repeatmasker/{specimen}.filt.alt.fa.masked",
#         "analysis/{refalias}/{analysis}/files/{mapper}_duplomap/repeatmasker/{specimen}.filt.alt.fa.cat"
#     log:
#         "logs/analysis/{refalias}/{analysis}/files/{mapper}_duplomap/repeatmasker/{specimen}.filt.log"
#     params:
#         outdir = "analysis/{refalias}/{analysis}/files/{mapper}_duplomap/repeatmasker/",
#         species = config['repeatmasker']['species'],
#         engine = config['repeatmasker']['engine']