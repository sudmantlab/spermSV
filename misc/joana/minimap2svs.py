import numpy as np
import pandas as pd
import os
import itertools
ref_paths = {
    'hg38': "/global/scratch/users/joana_rocha/PANPAN/reference/human_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fa",
    'clint': "/global/scratch/users/joana_rocha/PANPAN/reference/GCF_002880755.1_Clint_PTRv2_genomic.fa",
    'ht2t': "/global/scratch/users/joana_rocha/PANPAN/reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
}

sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', sep='\t')
sample_table_hprc=pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_names.txt', names=['names']) ### removed Y chromosome
references = ['hg38', 'clint']
samples = sample_table['Specimen'].unique()
human_samples =sample_table_hprc['names'].unique()
haps = ["hap1", "hap2"]
sample_haps = expand("{sample}.{hap}", sample=samples, hap=haps)
print(sample_haps)
print(human_samples)
all_samples = list(sample_haps) + list(human_samples)
print(all_samples)


rule all:
    input:
        #expand('svs/{ref}/{anything}/{anything}.vcf', ref=references, anything=all_samples),
        #'svs/survivor/panpanhprc_mapped2hg38_merged.vcf',
        #'svs/survivor/panpanhprc_mapped2clint_merged.vcf',
        #'svs/survivor/panpanhprc_mapped2hg38_merged_multisample_multicallers.vcf',
        #'svs/survivor/panpanhprc_mapped2clint_merged_multisample_multicallers.vcf',
        #'svs/survivor/panpanhprc_mapped2hg38_merged_multisample_multicallers_overlapp.txt',
        #'svs/survivor/panpanhprc_mapped2clint_merged_multisample_multicallers_overlapp.txt',
        #'svs/survivor/panpanhprc_mapped2hg38_merged_multisample_multicallers_overlapp.tiff',
        #'svs/survivor/panpanhprc_mapped2clint_merged_multisample_multicallers_overlapp.tiff',
        'svs/survivor/panpan_mapped2clint_merged_modified.vcf',
        'svs/survivor/hprc_mapped2hg38_merged_modified.vcf',
        'plink/panpan_mapped2clint_svs.eigenvec',
        'plink/panpan_mapped2clint_svs.eigenval',
        'plink/hprc_mapped2hg38_svs.eigenvec',
        'plink/hprc_mapped2hg38_svs.eigenval',
        #'plink/panpanhprc_mapped2hg38_svs.eigenvec',
        #'plink/panpanhprc_mapped2hg38_svs.eigenval'
   
        
rule panpan2ref:
    input:
        "/global/scratch/users/joana_rocha/PANPAN/Hifiasm-fasta_shortcut/{sample_hap}.p_ctg.fa",
    output: 
        'asm/{ref}/{sample_hap}.sam',
        'asm/{ref}/{sample_hap}.srt.bam',
        'asm/{ref}/{sample_hap}.srt.rg.bam',
    params:
        ref = lambda wildcards, output: ref_paths[wildcards.ref]
    shell: """
    minimap2 -a -x asm5 --cs -r2k -t 20 {params.ref} {input[0]} > {output[0]} &&  
    samtools sort -m4G -@4 -o {output[1]} {output[0]} &&
    samtools addreplacerg -r ID:{wildcards.sample_hap} -r PL:PACBIO -o {output[2]} {output[1]} 
    samtools index {output[2]}
    """

rule hprc2ref:
    input:
        "/global/scratch/users/joana_rocha/PANPAN/HPRC-yr1/HPRC_assemblies/{human_sample}.fa",
    output: 
        'asm/{ref}/{human_sample}.sam',
        'asm/{ref}/{human_sample}.srt.bam',
        'asm/{ref}/{human_sample}.srt.rg.bam',
    params:
        ref = lambda wildcards, output: ref_paths[wildcards.ref]
    shell: """
    minimap2 -a -x asm5 --cs -r2k -t 20 {params.ref} {input[0]} > {output[0]} &&  
    samtools sort -m4G -@4 -o {output[1]} {output[0]} &&
    samtools addreplacerg -r ID:{wildcards.human_sample} -r PL:PACBIO -o {output[2]} {output[1]} 
    samtools index {output[2]}
    """
    
rule asm2svs:
    input: 
        'asm/{ref}/{anything}.srt.rg.bam',
    output: 
        'svs/{ref}/{anything}/{anything}.vcf',
    params:
        ref = lambda wildcards, output: ref_paths[wildcards.ref],
        anything = lambda wildcards, output: wildcards.anything,
    shell: """
        svim-asm haploid 'svs/{wildcards.ref}/{wildcards.anything}' {input} {params.ref} &&
        mv 'svs/{wildcards.ref}/{wildcards.anything}/variants.vcf' 'svs/{wildcards.ref}/{wildcards.anything}/{wildcards.anything}.vcf' &&
        sed -i 's|Sample|{wildcards.anything}|g' 'svs/{wildcards.ref}/{wildcards.anything}/{wildcards.anything}.vcf'
    """


#rules to merge svim-asm sv-called individual vcfs: 
rule list_vcfs_panpan:
    output: 'svs/survivor/panpan_mapped2clint_vcf_files_raw_calls.txt'
    params:
        sample_haps = sample_haps  # Use the previously defined sample_haps
    run:
        shell("rm -f {output}")  # Remove existing output file if exists
        for sample in params.sample_haps:
            shell("ls svs/clint/{sample}/*.vcf >> {output}")

rule list_vcfs_hprc:
    output: 'svs/survivor/hprc_mapped2hg38_vcf_files_raw_calls.txt'
    params:
        human_samples = human_samples  # Use the previously defined
    run:
        shell("rm -f {output}")  # Remove existing output file if exists
        for sample in params.human_samples:
            shell("ls svs/hg38/{sample}/*.vcf >> {output}")

rule list_vcfs_panpanhprc2hg38:
    output: 'svs/survivor/panpanhprc_mapped2hg38_vcf_files_raw_calls.txt'
    params:
        all_samples = all_samples  # Use the previously defined
    run:
        shell("rm -f {output}")  # Remove existing output file if exists
        for sample in params.all_samples:
            shell("ls svs/hg38/{sample}/*.vcf >> {output}")

rule list_vcfs_panpanhprc2clint:
    output: 'svs/survivor/panpanhprc_mapped2clint_vcf_files_raw_calls.txt'
    params:
        all_samples = all_samples  # Use the previously defined
    run:
        shell("rm -f {output}")  # Remove existing output file if exists
        for sample in params.all_samples:
            shell("ls svs/clint/{sample}/*.vcf >> {output}")
           

rule merge_vcfs_panpan:
    input: 'svs/survivor/panpan_mapped2clint_vcf_files_raw_calls.txt'
    output: 'svs/survivor/panpan_mapped2clint_merged.vcf'
    shell: """
    SURVIVOR merge {input} 1000 1 1 0 0 30 {output}
    """

rule merge_vcfs_hprc:
    input: 'svs/survivor/hprc_mapped2hg38_vcf_files_raw_calls.txt'
    output: 'svs/survivor/hprc_mapped2hg38_merged.vcf'
    shell: """
    SURVIVOR merge {input} 1000 1 1 0 0 30 {output}
    """

rule merge_vcfs_panpanhprc2hg38:
    input: 'svs/survivor/panpanhprc_mapped2hg38_vcf_files_raw_calls.txt'
    output: 'svs/survivor/panpanhprc_mapped2hg38_merged.vcf'
    shell: """
    SURVIVOR merge {input} 1000 1 1 0 0 30 {output}
    """

rule merge_vcfs_panpanhprc2clint:
    input: 'svs/survivor/panpanhprc_mapped2clint_vcf_files_raw_calls.txt'
    output: 'svs/survivor/panpanhprc_mapped2clint_merged.vcf'
    shell: """
    SURVIVOR merge {input} 1000 1 1 0 0 30 {output}
    """

#rules to merge svim-asm+pbsv+sniffles sv-called multisample vcfs: 
rule list_vcfs_panpanhprc_survivor_pbsv_sniffles:
    output: 'svs/survivor/panpanhprc_{anything}_survivor_pbsv_sniffles_raw_calls.txt'
    shell: """
    ls svs/survivor/panpanhprc_{wildcards.anything}_merged* >> {output}
    """

rule merge_vcfs_panpanhprc_survivor_pbsv_sniffles:
    input: 'svs/survivor/panpanhprc_{anything}_survivor_pbsv_sniffles_raw_calls.txt'
    output: 'svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers.vcf'
    shell: """
    SURVIVOR merge {input} 1000 1 1 1 0 50 {output}
    """
rule make_VennDiagram_input:
    input: 'svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers.vcf'
    output: 'svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers_overlapp.txt'
    shell: """
    perl perl_script.pl {input} {output}
    """

rule venn_diagram:
    input: 
        t = "svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers_overlapp.txt"
    output: 
        "svs/survivor/panpanhprc_{anything}_merged_multisample_multicallers_overlapp.tiff"
    params:
        r_script = "venn_diagram.R"  
    shell: """
    Rscript {params.r_script} {input.t} {output}
    """


### recode vcfs for pca plots
rule recode_genotypes_hprc_mapped2hg38:
    input:
        "svs/survivor/hprc_mapped2hg38_merged.vcf"
    output:
        "svs/survivor/hprc_mapped2hg38_merged_modified.vcf"
    shell: "python recode_vcf.py {input} {output}"


rule recode_genotypes_panpan_mapped2clint:
    input:
        "svs/survivor/panpan_mapped2clint_merged.vcf"
    output:
        "svs/survivor/panpan_mapped2clint_merged_modified.vcf"
    shell: "python recode_vcf.py {input} {output}"


rule pcas_panpan:
    input:
        'svs/survivor/panpan_mapped2clint_merged_modified.vcf'
    output: 
        'plink/panpan_mapped2clint_svs.eigenvec',
        'plink/panpan_mapped2clint_svs.eigenval'
    params: 'plink/panpan_mapped2clint_svs'
    shell: """
    plink --vcf {input} --pca  --double-id --allow-extra-chr --out {params}
    """

rule pcas_hprc:
    input: 
        'svs/survivor/hprc_mapped2hg38_merged_modified.vcf'
    output: 
        'plink/hprc_mapped2hg38_svs.eigenvec',
        'plink/hprc_mapped2hg38_svs.eigenval',
    params: 'plink/hprc_mapped2hg38_svs'
    shell: """
    plink --vcf {input} --pca --double-id --allow-extra-chr --out {params}
    """


#rule pcas_panpan_hprc:
#    input:
#        'svs/survivor/panpanhprc_mapped2hg38_merged_modified.vcf'
#    output: 
#        'plink/panpanhprc_mapped2hg38_svs.eigenvec',
#        'plink/panpanhprc_mapped2hg38_svs.eigenval'
#    params: 'plink/panpanhprc_mapped2hg38_svs'
#    shell: """
#    plink --vcf {input} --pca  --double-id --allow-extra-chr --out {params}
#    """

#cat hprc_mapped2hg38_merged.vcf | sed 's|./.|0/0|g' | sed 's|./1|0/1|g' | sed 's|1/.|1/0|g' > hprc_mapped2hg38_merged_modified.vcf
#cat panpan_mapped2clint_merged.vcf | sed 's|./.|0/0|g' | sed 's|./1|0/1|g' | sed 's|1/.|1/0|g' > panpan_mapped2clint_merged_modified.vcf
#cat panpanhprc_mapped2hg38_merged.vcf | sed 's|./.|0/0|g' | sed 's|./1|0/1|g' | sed 's|1/.|1/0|g' > panpanhprc_mapped2hg38_merged_modified.vcf


#Max distance between breakpoints: 1000
#Minimum number of supporting caller: 1
#Take the type into account (1==yes, else no): 1
#Take the strands of SVs into account (1==yes, else no): 0
#Estimate distance based on the size of SV (1==yes, else no): 0
#Minimum size of SVs to be taken into account: 30


#After the merge, you can use the merged vcf file to obtain a graphical representation of the overlaps.
#If you have many VCF files that you want to compare: The first step is to obtain a pairwise comparison matrix like this:


#rule bgzip_vcfs:
#    input:
#        'svs/{ref}/{anything}/{anything}_variants.vcf',
#    output:
#        vcf_gz='svs/{ref}/{anything}/{anything}.vcf.gz',
#        vcf_gz_tbi='svs/{ref}/{anything}/{anything}.vcf.gz.tbi',
#    shell: """
#    bgzip -c {input} > {output.vcf_gz} &&
#    tabix -p vcf {output.vcf_gz}
#    """