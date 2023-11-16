
# globus login should already be done
# run this script to recursively download new data without having to deal with symlinks via the browser system

ep1='86147ca7-94db-45c3-bffc-e80d3df16ee5' # hudsonalpha
ep2='8dfcaab6-fe4e-4376-9430-785631512d4a' # scratch storage

globus transfer 86147ca7-94db-45c3-bffc-e80d3df16ee5:/psudmant@berkeley.edu/Human_20_062723 \
8dfcaab6-fe4e-4376-9430-785631512d4a:/global/scratch/users/stacy-l/spermSV/globus_data \
--label "HudsonAlpha Data Sync" \
--recursive --sync-level exists --verbose