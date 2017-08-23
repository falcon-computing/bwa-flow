tmp_dir=/pool/local/temp
ref_dir=/pool/local/ref
input_dir=/pool/storage/fastq/benchmarks
ref_genome=$ref_dir/human_g1k_v37.fasta

bwa=../bin/bwa
fpga_path=../bin/sw_17p1.xclbin
bwa_options="--use_fpga \
             --max_fpga_thread=1 \
             --fpga_path=$fpga_path \
             --chunk_size=5000"
