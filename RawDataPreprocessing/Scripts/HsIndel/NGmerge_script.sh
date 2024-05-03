#!/bin/sh 

NGmerge -d -p 0.105 -1 HI_SABnegA_5trimmed.1.fastq.gz -2 HI_SABnegA_5trimmed.2.fastq.gz -o HI_SABnegA_trimmed_merged.fastq.gz -f failedmerge -l mergelog.txt -c dovetaillog.txt