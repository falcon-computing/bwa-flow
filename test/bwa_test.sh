#!/bin/bash
SCRIPT_DIR=/curr/yaoh/prog/bwa-flow/test
source $SCRIPT_DIR/../setup.sh
index=0
rm $SCRIPT_DIR/testing_files/result_bams/*.bam
rm $SCRIPT_DIR/testing_files/error_logs/*.log
rm $SCRIPT_DIR/testing_files/bwa_stderrs/*.txt
rm $SCRIPT_DIR/testing_files/temp_part_files/*

while IFS='' read -r line; do
    index=$[$index+1]
    test_bam_out=$SCRIPT_DIR/testing_files/result_bams/Result_of_line_$index.bam
    error_log_out=$SCRIPT_DIR/testing_files/error_logs/Error_log_of_line_$index.log
    bwa_stderr=$SCRIPT_DIR/testing_files/bwa_stderrs/bwa_stderr$index.txt
    golden_bam_in=`echo $line | awk '{print $NF}'`
    sort_option=`echo $line | awk '{print $1}'`
    sort_ref="--sort=true"
    echo "Executing command for line "$index
    echo $line | awk '{$NF="";print "'$SCRIPT_DIR'/../bin/bwa mem "$0"" }'
    `echo $line | awk '{$NF="";print "'$SCRIPT_DIR'/../bin/bwa mem "$0"" }'` 2>$bwa_stderr
    echo "Finished bwa computing, start processing bam"
    if [ "$sort_option" = "$sort_ref" ]; 
    then
      echo "Sort option is used so we do samtools merge"
      /curr/yaoh/prog/dnaseq/samtools-1.3/samtools merge -c -p -l 0 -@ 8 $test_bam_out $SCRIPT_DIR/testing_files/temp_part_files/part*
    else
      echo "Sort option is not used so we do samtools merge and sort"
      /curr/yaoh/prog/dnaseq/samtools-1.3/samtools merge -c -p -l 0 -@ 8 $SCRIPT_DIR/testing_files/result_bams/merged.bam $SCRIPT_DIR/testing_files/temp_part_files/part*
      cat $SCRIPT_DIR/testing_files/result_bams/merged.bam | /curr/yaoh/prog/dnaseq/samtools-1.3/samtools sort -l 0 -m 8g -@ 8  -o $test_bam_out
    fi
    rm $SCRIPT_DIR/testing_files/temp_part_files/part*
    echo "Comparing result bam with golden_bam"
    /curr/mhhuang/tools/bamUtil_1.0.13/bamUtil/bin/bam diff --in1 $test_bam_out --in2 $golden_bam_in >$error_log_out
    filesize=`wc -c <$error_log_out`
    # check the results
    if [ "$filesize" -gt 0 ]
    then
      echo "Error filesize is "$filesize "and test "$index" is failed,please check the "$SCRIPT_DIR"/testing_files/error_logs/ dir for detailed difference"
    else
      echo "Error filesize is "$filesize "and test "$index" is passed"
      #remove the bam and log files if test passed
      rm $test_bam_out
      rm $error_log_out
      rm $bwa_stderr
    fi
    rm $SCRIPT_DIR/testing_files/result_bams/merged.bam
done < "$1"

echo "All test finished"
