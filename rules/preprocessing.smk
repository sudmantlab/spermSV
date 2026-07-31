# Usage notes:
# HiFiAdapterFilt is not managed under bioconda – it needs to be installed by cloning from the repo directly.
# I've made changes under a new fork (stacy-l/HiFiAdapterFilt) to fix a couple of directory path issues w/ temp files.
# that seemed to be specific to the bam input code. nicolas submitted a pull request for this, which hasn't
# been accepted yet. This is likely fine because our samples all appear to be adapter-free (Revio pre-filtering).

# below: adds HiFiAdapterFilt code directory to path if not already present
repo_path = config['paths']['HiFiAdapterFilt']['repo']
db_path = config['paths']['HiFiAdapterFilt']['db']

if os.system("echo $PATH | grep /HiFiAdapterFilt") == 256:
    os.environ["PATH"] += os.pathsep + os.pathsep.join([repo_path, db_path])


rule uBAMtoFastq:
    # Takes the unaligned BAM file from the sequencing facility and transforms it into a fastq.
    # This transformation retains the read quality (rq) + MM & ML (methylation information) tags for each read, 
    # such that they can be carried forward through the mapping step & recorded in the final BAM.
    input:
        "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.ccs.bam"
    output:
        "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
    conda: "../envs/preprocessing.yml"
    threads: 10
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T MM,ML {input} -0 {output}
        """

# rule HiFiAdapterFilt:
#     # Notes:
#     # This script must be executed from within the input file directory.
#     # It *can* write outputs to a specified directory
#     # Thus, values for below params are determined "dynamically" from input/output paths,
#     # in order to guard against issues arising from input/output path restructuring.
#     # For example, the below params evaluate to:
#     # outDir = 'output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}'
#     # inDir = 'data/PacBio-HiFi/{specimen}/{lane}'
#     # inPref ='{smrtcell}.ccs'
#     input: 
#         "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
#     output:
#         "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
#     log: "logs/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.log"
#     params:
#         # outDir = lambda wildcards, output: os.path.dirname(output[0]),
#         # inDir = lambda wildcards, output: os.path.dirname(input[0]),
#         # inPref = lambda wildcards, output: re.sub('(?<=ccs).*', '', os.path.basename(input[0])),
#         outDir = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}",
#         inDir = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}",
#         inPref = "{smrtcell}.ccs"
#     conda: "../envs/preprocessing.yml"
#     threads: 10
#     shell: 
#         """
#         export PATH=`echo $PATH | tr ":" "\n" | grep -v "sl-7.x86_64" | tr "\n" ":"`

#         ROOTPROJDIR=$(pwd -P)
#         cd {params.inDir}
#         pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
#         cd -
#         """