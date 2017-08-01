### Compilation
```
> mv ./config.mk.template ./config.mk # with corresponding modifications
```
1. Setting up dependencies    

    - Build FPGA: add `BUILD_FPGA=1` in *config.mk*
    - Makefile will automatically read `$ALTERAOCLSDKROOT` and `$XILINX_SDX` and determine with platform to build
    - Build MPI: add `OPENMPI_DIR` in *config.mk*
    - Add Falcon LM license: set `FLMDIR` in *config.mk*

3. Compile, this will compile bwa, kflow and bwa-flow automatically
```
> make all
```

### Execution

```
> mv bin/run.sh.template bin/run.sh # with corresponding modifications
```
1. Normal CPU
    - use default settings in *bin/bwa* should be enough
    - add `--output_dir` for output folder name
    - Synopsis:     
    ```
    > ./bin/run.sh mem <path-to-ref> <path-to-fastq1> <path-to-fastq2> 
    ```
2. With FPGA
    - use settings for CPU run
    - add `--use_fpga` in the options in *bin/run.sh*
    - add `--fpga_path=<path-to-bitstream>` in the options in *bin/run.sh*
3. MPI
    - TODO
