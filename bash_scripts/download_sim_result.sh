sim_no=$1
file_name="sim_res_$sim_no.rds"
echo "downloading simulation result $file_name" 

scp timbar@hpcc.wharton.upenn.edu:~/data/projects/glmeiv/public/simulations/results/"$file_name" /Users/timbarry/research_offsite/projects/glmeiv/public/simulations/results