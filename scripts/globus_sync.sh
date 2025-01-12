
# globus login should already be done
# run this script to recursively download new data without having to deal with symlinks via the browser system

ep1='86147ca7-94db-45c3-bffc-e80d3df16ee5' # hudsonalpha
ep2='8dfcaab6-fe4e-4376-9430-785631512d4a' # scratch storage

# Get the list of directories
directories=$(globus ls 86147ca7-94db-45c3-bffc-e80d3df16ee5:/psudmant@berkeley.edu/Human_20_062723 | grep /)

# Loop through each directory
for dir in $directories; do
    # Get the list of files matching the pattern
    files=$(globus ls "${ep1}:/psudmant@berkeley.edu/Human_20_062723/${dir}" | grep '\.hifi_reads\.[^/]*\.bam$')
    
    # Loop through each matching file and initiate a transfer
    for file in $files; do
        echo "Transferring ${file} from ${dir}"
        globus transfer \
        "${ep1}:/psudmant@berkeley.edu/Human_20_062723/${dir}${file}" \
        "${ep2}:/global/scratch/users/stacy-l/spermSV/globus_data/${dir}${file}" \
        --label "HudsonAlpha Data Sync - ${file}" \
        --sync-level mtime --verbose
    done
done