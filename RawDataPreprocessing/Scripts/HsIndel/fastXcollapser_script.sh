#!/bin/sh 

gunzip HD_SABinA_trimmed_merged_filterC.fastq.gz

fastx_collapser -v -i HD_SABinA_trimmed_merged_filterC.fastq -o HD_SABinA_collapsed.fasta