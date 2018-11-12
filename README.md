### Jenkins Build-Status
AWS [![Build Status](http://us-1.falcon-computing.com:8080/buildStatus/icon?job=Falcon-Build-BWA-FLOW)](http://us-1.falcon-computing.com:8080/job/Falcon-Build-BWA-FLOW/)
HWC [![Build Status](http://us-1.falcon-computing.com:8080/buildStatus/icon?job=HWC-BWA-FLOW)](http://us-1.falcon-computing.com:8080/view/HWC/job/HWC-BWA-FLOW/)
Local [![Build Status](http://us-1.falcon-computing.com:8080/buildStatus/icon?job=Falcon-Build-BWA-FLOW-Internal)](http://us-1.falcon-computing.com:8080/job/Falcon-Build-BWA-FLOW-Internal/)

### Compilation
1. Use cmake to ompile. This will compile bwa, kflow and bwa-flow automatically
```
> mkdir build && cd build
> cmake ..
> make all
> make CTEST_OUTPUT_ON_FAILURE=1 test
```
2. Setting up options

    - Build with HTS : add `-DUSE_HTSLIB=On` after `cmake` command, by default this option is on
    - Build with FPGA: add `-DUSE_FPGA=On` after `cmake` command, by default this option is if OpenCL Runtime library found
    - BUILD with MPI : add `-DUSE_MPI=On` after `cmake` command, by default this option is off
    - Build in Debug mode: add `-DCMAKE_BUILD_TYPE=Debug` after `cmake` command (default build type)
    - Build in Release mode: add `-DCMAKE_BUILD_TYPE=Release` after `cmake` command
    - Select deployment type: add `-DDEPLOYMENT_DST=aws` or `-DDEPLOYMENT_DST=hwc` after `cmake` command



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

### Testing

Use `test/run-test.sh` to launch unit tests and sample tests. The testing script can only prepare the bitstream and data samples.
The reference file (.fasta) needs setting up manually, and `test/test-config.sh` should be adjusted accordingly. 
