#Map reads to assembly and assembly to self
./purge_dups/src/split_fa input.fasta > split.fasta

./minimap2 -I 200G -t 24 -xasm5 -DP split.fasta split.fasta > split.genome.paf

./minimap2 -I 200G -x map-pb -t 24 split.fasta reads.ccs.fastq.gz > reads.paf

#Calculate haploid/diploid coverage threshold and remove haplotype duplicates from assembly
./purge_dups/src/pbcstat -O coverage reads.paf

./purge_dups/src/calcuts PB.stat > cutoffs

./purge_dups/src/purge_dups -2 -c PB.base.cov -T cutoffs split.genome.paf > dups.bed

./purge_dups/src/get_seqs -e -p asm_mTadBra.purged dups.bed input.fasta

