tmp_dir=/pool/local/temp
ref_dir=/pool/local/ref
input_dir=/pool/storage/fastq/benchmarks
ref_genome=$ref_dir/human_g1k_v37.fasta

bwa=../bin/bwa
fpga_path=/curr/yaoh/kernels/sw_4v0_3x20_dup.xclbin
pac_path=/pool/storage/yaoh/human_g1k_v37.fasta.pac
bwa_options="--use_fpga \
             --max_fpga_thread=1 \
             --fpga_path=$fpga_path \
             --pac_path=$pac_path \
             --chunk_size=10000"
