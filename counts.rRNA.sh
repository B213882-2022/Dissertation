align_dir="./align"
gtf_dir="./GRCm39.rRNA.gtf"
dirs=$(find ${align_dir} -name "*.out.bam" | tr '\n' ' ')
featureCounts -a ${gtf_dir} -o counts.rRNA.txt -T 60 -t exon -g gene_id --fraction -O -M ${dirs}
