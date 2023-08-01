import os
import subprocess
import pandas as pd

# Usage notes:
# HiFiAdapterFilt is not managed under bioconda – it needs to be installed by cloning from the repo directly
# note that I've made changes under a new fork to fix a couple of directory path issues w/ temp files
# that seemed to be specific to the bam input code. nicolas submitted a pull request for this, which hasn't
# been accepted yet...
# the fix is under a fork (stacy-l/HiFiAdapterFilt) – to be contributed to main later if this fixes the issue?

# below: adds HiFiAdapterFilt code directory to path if not already present
root_path = subprocess.run(["pwd", "-P"], stdout=subprocess.PIPE).stdout.decode('utf-8').strip('\n')
repo_path = "{root_dir}/code/HiFiAdapterFilt".format(root_dir = root_path)
db_path = "{root_dir}/code/HiFiAdapterFilt/DB".format(root_dir = root_path)

if os.system("echo $PATH | grep /HiFiAdapterFilt") == 256:
    os.environ["PATH"] += os.pathsep + os.pathsep.join([repo_path, db_path])

rule HiFiAdapterFilt:
    version: subprocess.check_output("pbadapterfilt.sh --version", shell=True)
    input: "data/PacBio-HiFi/homo_sapiens/{specimen}/{lane}/{smrtcell}.ccs.bam"
    output: 
        filtered = "output/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz",
        stats = "output/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.stats",
    log: "logs/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.bamfilt.log"
    params:
        outDir = lambda wildcards, output: os.path.dirname(output[0]),
        inBaseName = lambda wildcards, input: os.path.basename(input[0]),
        inDir = lambda wildcards, input: os.path.dirname(input[0]),
        inPref = lambda wildcards, input: os.path.splitext(os.path.basename(input[0]))[0],        
    conda: "../envs/HiFiAssembly.yml"
    threads: 10
    shell: """
        ROOTPROJDIR=$(pwd -P)
        cd {params.inDir}
        pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
        cd -
    """