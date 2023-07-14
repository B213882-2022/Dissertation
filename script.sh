# 从GEO获取SRA数据并转换为fastq格式
target_dir="./raw_data"
mkdir -p ${target_dir}/sra
parallel --verbose -j 20 prefetch -O ${target_dir}/sra ::: $(awk 'BEGIN{FS=","}{if(NR>1){print $1}}' annotation.csv)
mkdir -p ${target_dir}/fastq
parallel --verbose -j 20 fastq-dump --split-3 ${target_dir}/sra/{} -O ${target_dir}/fastq ::: $(ls ${target_dir}/sra | sort)

# trimmomatic修剪reads（单端测序）
tool_dir="./tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapter_dir="./tools/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
fastq_dir="./raw_data/fastq"
output_dir="./trimmed_data"
mkdir -p ${output_dir}/fastq
parallel --verbose -j 8 java -jar ${tool_dir} SE -threads 16 -phred33 ${fastq_dir}/{} ${output_dir}/fastq/{} ILLUMINACLIP:${adapter_dir}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ::: $(ls ${fastq_dir} | sort)

# fastqc + multiqc
fastq_dir="./trimmed_data/fastq"
target_dir="./trimmed_data"
mkdir -p ${target_dir}/fastqc
mkdir -p ${target_dir}/multiqc
parallel --verbose -j 8 fastqc -t 16 -o ${target_dir}/fastqc ${fastq_dir}/{} ::: $(ls ${fastq_dir} | sort )
multiqc ${target_dir}/fastqc -o ${target_dir}/multiqc

# STAR进行alignment（单端）
fastq_dir="./trimmed_data/fastq"
#构建函数
run_star(){
tool_dir="./tools/STAR_2.7.10b/Linux_x86_64_static/STAR"
fastq_dir="./trimmed_data/fastq"
fastq_suffix=".fastq"
ref_index_dir="./ref_index_GRCm39"
output_dir="./align"
name=$1
mkdir -p ${output_dir}/${name/${fastq_suffix}/}
${tool_dir} --runThreadN 16 --genomeDir ${ref_index_dir} --outFileNamePrefix ${output_dir}/${name/${fastq_suffix}/}/${name/${fastq_suffix}/}_ --outSAMtype BAM SortedByCoordinate --readFilesIn ${fastq_dir}/${name} 
# 如果是gz文件则还需要加上参数 --readFilesCommand zcat
}
export -f run_star
parallel --verbose -j 8 run_star {} ::: $(ls ${fastq_dir} | sort)

# featureCounts
align_dir="./align"
gtf_dir="./GRCm39.gtf"
dirs=$(find ${align_dir} -name "*.out.bam" | tr '\n' ' ')
featureCounts -a ${gtf_dir} -o counts.txt -T 60 ${dirs}
