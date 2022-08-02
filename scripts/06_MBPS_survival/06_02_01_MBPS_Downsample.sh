#Script for downsampling raw MBPS data from bam files to  1K, 10K, 50K and 100K with for loop

#bam files too large for github, available at users request

#sequencing run 1
module load samtools/1.12
# path

p= *user_defined_file_path*
project=*user_defined_project name*
source="$p/$project/aligned"

# Number of reads to extract
for N in 1000 10000 50000 100000; do
destination="$p/${project}_${N}/aligned"
mkdir -p $destination
chrom="lethal_CACNA2D4::chr12:1906318-1906785"
for s in `ls $source`; do
echo $s;
# in bam path
inbam="$source/$s/${s}.bismark.bam";
# create destination path
mkdir -p "$destination/$s";
# out bam path
outbam="$destination/$s/${s}.bismark.bam";
# Calculate fraction to extract the N reads
frac=$( samtools idxstats $inbam | grep "lethal_CACNA2D4::chr12:1906318-1906785"| cut -f3 | awk -v N=$N '{frac=N/$0; if(frac>1){print 1}else{print frac}}' );
# Downsampling
if [[ ${frac%.*} -eq 1 ]]; then 
samtools view -b -F 4 -q 40 $inbam "$chrom:0-467" > $outbam;
else 
samtools view -b -F 4 -q 40 -s $frac $inbam "$chrom:0-467" > $outbam;
fi
samtools index $outbam;
fastq="$p/${project}_${N}/raw/$s";
fastq_trimmed="$p/${project}_${N}/raw_trimmed/$s";
mkdir -p $fastq;
mkdir -p $fastq_trimmed;
# samtools fastq -@ 8 -0 $fastq/${s}_R1.fastq -2 /dev/null -1 /dev/null -s /dev/null -n $outbam;
samtools fastq -@ 8 $outbam > $fastq/${s}_R1.fastq;
samtools fastq -@ 8 $outbam > $fastq_trimmed/${s}_trimmed_R1.fastq; 
gzip $fastq/${s}_R1.fastq;
gzip $fastq_trimmed/${s}_trimmed_R1.fastq;
done
done

### Use these fastq files as input for MethPanel to get the DNA methylation beta values


#Sequencing run 2

module load samtools/1.12
# path
p= *user_defined_file_path*
project=*user_defined_project name*
source="$p/$project/aligned"

# Number of reads to extract
for N in 1000 10000 50000 100000; do
destination="$p/${project}_${N}/aligned"
mkdir -p $destination
chrom="lethal_CACNA2D4::chr12:1906318-1906785"
for s in `ls $source`; do
echo $s;
# in bam path
inbam="$source/$s/${s}.bismark.bam";
# create destination path
mkdir -p "$destination/$s";
# out bam path
outbam="$destination/$s/${s}.bismark.bam";
# Calculate fraction to extract the N reads
frac=$( samtools idxstats $inbam | grep "lethal_CACNA2D4::chr12:1906318-1906785"| cut -f3 | awk -v N=$N '{frac=N/$0; if(frac>1){print 1}else{print frac}}' );
# Downsampling
if [[ ${frac%.*} -eq 1 ]]; then 
samtools view -b -F 4 -q 40 $inbam "$chrom:0-467" > $outbam;
else 
samtools view -b -F 4 -q 40 -s $frac $inbam "$chrom:0-467" > $outbam;
fi
samtools index $outbam;
fastq="$p/${project}_${N}/raw/$s";
fastq_trimmed="$p/${project}_${N}/raw_trimmed/$s";
mkdir -p $fastq;
mkdir -p $fastq_trimmed;
# samtools fastq -@ 8 -0 $fastq/${s}_R1.fastq -2 /dev/null -1 /dev/null -s /dev/null -n $outbam;
samtools fastq -@ 8 $outbam > $fastq/${s}_R1.fastq;
samtools fastq -@ 8 $outbam > $fastq_trimmed/${s}_trimmed_R1.fastq; 
gzip $fastq/${s}_R1.fastq;
gzip $fastq_trimmed/${s}_trimmed_R1.fastq;
done
done

### Use these fastq files as input for MethPanel to get the DNA methylation beta values
