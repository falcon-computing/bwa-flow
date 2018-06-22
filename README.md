### Compilation
1. Use cmake to ompile. This will compile bwa, kflow and bwa-flow automatically
```
> mkdir build && cd build
> cmake ..
> make all
```
2. Setting up options

    - Build with FPGA: add `-DUSE_FPGA=On` after `cmake` command
    - Build with HTS : add `-DUSE_HTSLIB=On` after `cmake` command
    - Build in Debug mode: add `-DCMAKE_BUILD_TYPE=Debug` after `cmake` command



### Execution

```
> mv bin/run.sh.template bin/run.sh # with corresponding modifications
```
1. Normal CPU
    - use default settings in *bin/bwa* should be enough
    - add `--output_dir` for output folder name
    - Synopsis:     
    ```
    > ./bin/run.sh <path-to-ref> <path-to-fastq1> <path-to-fastq2> 
    ```
2. With FPGA
    - use settings for CPU run
    - add `--use_fpga` in the options in *bin/run.sh*
    - add `--fpga_path=<path-to-bitstream>` in the options in *bin/run.sh*
3. MPI
    - TODO
