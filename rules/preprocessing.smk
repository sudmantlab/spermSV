# Usage notes:
# HiFiAdapterFilt is not managed under bioconda – it needs to be installed by cloning from the repo directly
# note that I've made changes under a new fork to fix a couple of directory path issues w/ temp files
# that seemed to be specific to the bam input code. nicolas submitted a pull request for this, which hasn't
# been accepted yet...
# the fix is under a fork (stacy-l/HiFiAdapterFilt) – to be contributed to main later if this fixes the issue?

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
        temp("output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz")
    conda: "../envs/preprocessing.yml"
    threads: 10
    shell:
        # extremely janky workaround right now because the savio module of samtools is OLD but also unremovable from my path rn
        # without forcing it like below
        # https://unix.stackexchange.com/questions/108873/removing-a-directory-from-path
        # there's some WEIRD SHIT going on right now where there's conda (???) default on savio in this path
        # what is HAPPENING
        # /global/software/sl-7.x86_64/modules/langs/python/3.9/envs
        """
        export PATH=`echo $PATH | tr ":" "\n" | grep -v "sl-7.x86_64" | tr "\n" ":"`
        samtools fastq -@ {threads} -c 6 -T MM,ML {input} -0 {output}
        """


rule HiFiAdapterFilt:
    # Notes:
    # This script must be executed from within the input file directory.
    # It *can* write outputs to a specified directory
    # Thus, values for below params are determined "dynamically" from input/output paths,
    # in order to guard against issues arising from input/output path restructuring.
    # For example, the below params evaluate to:
    # outDir = 'output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}'
    # inDir = 'data/PacBio-HiFi/{specimen}/{lane}'
    # inPref ='{smrtcell}.ccs'
    input: 
        "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
    output:
        "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz",
        "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.stats",
    log: "logs/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.log"
    params:
        # outDir = lambda wildcards, output: os.path.dirname(output[0]),
        # inDir = lambda wildcards, output: os.path.dirname(input[0]),
        # inPref = lambda wildcards, output: re.sub('(?<=ccs).*', '', os.path.basename(input[0])),
        outDir = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}",
        inDir = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}",
        inPref = "{smrtcell}.ccs"
    conda: "../envs/preprocessing.yml"
    threads: 10
    shell: 
        """
        export PATH=`echo $PATH | tr ":" "\n" | grep -v "sl-7.x86_64" | tr "\n" ":"`

        ROOTPROJDIR=$(pwd -P)
        cd {params.inDir}
        pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
        cd -
        """