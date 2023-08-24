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

# output/preprocessing/HiFiAdapterFilt/sudmant/894/PBmixRevio1412_1_D01_PEWA_24hours_19kbExpressCCSv3190pM2hrPE_200pM_HumanSudmant1894_CCSExpressIndex/m84053_230601_202536_s4.hifi_reads.bc2050.ccs.filt.fastq.gz

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
    input: "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.ccs.bam"
    output:
        filtered = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz",
        stats = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.stats",
    log: "logs/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.bamfilt.log"
    params:
        outDir = lambda wildcards, output: os.path.dirname(output[0]),
        inDir = lambda wildcards, input: os.path.dirname(input[0]),
        inPref = lambda wildcards, input: os.path.splitext(os.path.basename(input[0]))[0],        
    conda: "../envs/HiFiAssembly.yml"
    threads: 10
    shell: 
        """
        ROOTPROJDIR=$(pwd -P)
        cd {params.inDir}
        pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
        cd -
        """