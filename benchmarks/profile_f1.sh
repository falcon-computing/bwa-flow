tmp_dir=/local/temp
ref_dir=/local/ref
input_dir=/local/benchmarks
ref_genome=$ref_dir/human_g1k_v37.fasta

bwa=../bin/bwa
fpga_path=../bin/sw-wref.awsxclbin
pac_path=/local/ref/human_g1k_v37.fasta.fpga.pac
bwa_options="--use_fpga \
             --max_fpga_thread=1 \
             --extra_thread=0 \
             --fpga_path=$fpga_path \
             --pac_path=$pac_path \
             --chunk_size=5000"
