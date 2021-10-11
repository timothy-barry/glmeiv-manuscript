# sets up the offsite directory structure
source ~/.research_config

# 1. Create a nextflow config file, which puts the "work"
# directory in the top-level glmeiv data directory.
if [ -f nextflow.config ]; then
  rm nextflow.config
fi
touch nextflow.config
echo workDir = \"$LOCAL_GLMEIV_DATA_DIR\work\" >> nextflow.config

# 2. Initialize several directories offsite:
# - simulations/{spec_objects, results}
# - gasperini/{data, results}
# - xie/{data, results}
public_dir=$LOCAL_GLMEIV_DATA_DIR"public/"
mkdir -p $public_dir"simulations/spec_objects" \
$public_dir"simulations/results" \
$public_dir"gasperini/data" \
$public_dir"gasperini/results" \
$public_dir"xie/data" \
$public_dir"xie/results" \
$public_dir"aux"
