import parsing_utils as parse
from parsing_utils import SnifflesVCF as sVCF

def write_fasta(vcf, outfile):
    # Takes the selected ALT sequences writes a fasta file for input through RepeatMasker.
    table = vcf[['ID', 'ALT']]
    with open(outfile, 'w') as f:
        for header, seq in table.to_records(index=False):
            f.write(f'>{header}\n{seq}\n')

vcf = sVCF().create(snakemake.input['vcf'], header_len = 247)
selected = vcf.query("SVTYPE not in ['BND', 'DEL', 'DUP'] & ALT not in ['N', '<INS>']")
print("Writing selected ALT sequences to ", snakemake.output['fa'])
write_fasta(selected, snakemake.output['fa'])