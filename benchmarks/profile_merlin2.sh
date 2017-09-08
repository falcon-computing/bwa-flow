tmp_dir=/pool/ssd1/temp
ref_dir=/space/scratch/genome/ref
input_dir=/space/scratch/genome/fastq/benchmarks
ref_genome=$ref_dir/human_g1k_v37.fasta

bwa=../bin/bwa
fpga_path=../bin/sw.aocx
bwa_options="--use_fpga \
             --max_fpga_thread=1 \
             --extra_threads=-1 \
             --fpga_path=$fpga_path \
             --t=16 \
             --chunk_size=5000"
