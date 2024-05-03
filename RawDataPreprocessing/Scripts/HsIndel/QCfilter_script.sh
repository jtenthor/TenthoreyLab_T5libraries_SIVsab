#!/bin/sh 

gunzip HD_SABinA_trimmed_merged.fastq.gz

fastq_quality_filter -v -q 15 -p 100 -z -i HD_SABinA_trimmed_merged.fastq -o HD_SABinA_trimmed_merged_filterA.fastq.gz

fastq_quality_filter -v -q 20 -p 100 -z -i HD_SABinA_trimmed_merged.fastq -o HD_SABinA_trimmed_merged_filterB.fastq.gz

fastq_quality_filter -v -q 25 -p 100 -z -i HD_SABinA_trimmed_merged.fastq -o HD_SABinA_trimmed_merged_filterC.fastq.gz

fastq_quality_filter -v -q 30 -p 100 -z -i HD_SABinA_trimmed_merged.fastq -o HD_SABinA_trimmed_merged_filterD.fastq.gz