#Script for processing MBPS data through MethPanel https://github.com/thinhong/MethPanel

# Please replace "/user/define/path/" by the path in your computer

# Name the project directory
project="ProstateLethal"

# prior to running the following commands download fastq files into new folder called 'raw' in "/user/define/path/${project}/MBPS/raw/"

### 1. Create sample config file
path="/user/define/path/${project}/MBPS/raw/"  # this directory contains fastq files
mkdir -p "/user/define/path/${project}/MBPS/config"
sample_config="/user/define/path/${project}/MBPS/config/sample.${project}.pre.config"

echo -e "[ $project ]\n" > $sample_config
for i in `ls $path`; do
echo -e "\t[[ $i ]]" >> $sample_config;
echo >> $sample_config;
done

### 2. Create system config file
system_config="/user/define/path/${project}/MBPS/config/system.${project}.pre.config"


### 3. Run MethPanel: Please replace "/path/to/MethPanel" by the path which MethPanel was installed
module load python/2.7.11
module load java/jdk1.8.0_60
module load bpipe/0.9.9.7

python "/path/to/MethPanel/pipe/run_Bpipe.py" $sample_config $system_config
