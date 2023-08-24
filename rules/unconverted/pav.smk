

# only necessarily if using library:// remote protocol
# # checks for SylabsCloud remote endpoint access
# # enables library:// protocol on savio apptainer
# # removal: apptainer remote remove SylabsCloud
# import os
# status = subprocess.run("cat /global/home/users/stacy-l/.apptainer/remote.yaml", shell = True, stdout = subprocess.PIPE).stdout.decode('utf-8').split('\n')
# if "Active: SylabsCloud" in status:
#     pass
# else:
#     subprocess.run("apptainer remote add --no-login SylabsCloud cloud.sycloud.io", shell = True)
#     subprocess.run("apptainer remote use SylabsCloud", shell = True)

rule pav_setup:
    # Creates symlinks for pav run.
    input:
        config = "/global/scratch/users/stacy-l/spermSV/config/pav/config.json",
        asm_table = "/global/scratch/users/stacy-l/spermSV/config/pav/assemblies.tsv",
    output:
        config = "output/hg38_no_alts/pav/config.json",
        asm_table = "output/hg38_no_alts/pav/assemblies.tsv",
    params:
        out_dir = "output/hg38_no_alts/pav",
    shell:
        """
        mkdir -p {params.out_dir}
        ln -s {input.config} {output.config}
        ln -s {input.asm_table} {output.asm_table}
        """

rule pav:
    input:
        config = "output/hg38_no_alts/pav/config.json",
        asm_table = "output/hg38_no_alts/pav/assemblies.tsv",
    output:
        # placeholder for run finish right now
        flag = "output/hg38_no_alts/pav/run.success"
    threads: 40
    params:
        mount_dir = "/global/scratch/users/stacy-l",
        out_dir = "output/hg38_no_alts/pav",
        cache_dir = "/global/scratch/users/stacy-l/software/singularity"
    shell:
        """
        cd {params.out_dir}
        APPTAINER_CACHEDIR={params.cache_dir}

        # DO NOT USE -c OPTION. singularity run interprets as "contain"
        singularity run --writable-tmpfs \
        --bind {params.mount_dir}:/{params.mount_dir} \
        docker://becklab/pav:latest \
        --cores {threads}

        cd -
        touch {output.flag}
        """