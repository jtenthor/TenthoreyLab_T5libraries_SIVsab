# Script to process H2R combinatorial library sequencing data
# Experiment: H2R SIVsab GOF

## load Biostrings module
library(Biostrings)

#ExpectedSeqs = read.csv("../Data/CombiSeqList.csv") 
ExpectedSeqsAll = read.csv("../Data/CombiSeqListActual_H.csv")

## define function to turn pre-processed fasta file into table
fasta2counts = function(file) {
  fa = readDNAStringSet(file) #read in fasta file: gives width (int), names (char), and seq (not char!)
  fa_headers_split = strsplit(names(fa),"-") #split the seq names into list of 2-element list (current format as ID#-counts)
  output = data.frame(sequence = as.character(fa), #make data frame with extracted DNA seqs
                      counts = as.integer(sapply(fa_headers_split,"[[",2))) #and counts = 2nd element w/in sublist (+ change to int class)
  return(output)
}

## run fxn on each library to get table
inA = fasta2counts("../Data/H2R_SABinA_preprocessedQ25.fasta") 
colnames(inA)[2] = "inA"
inB = fasta2counts("../Data/H2R_SABinB_preprocessedQ25.fasta") 
colnames(inB)[2] = "inB"
negA = fasta2counts("../Data/H2R_SABnegA_preprocessedQ25.fasta") 
colnames(negA)[2] = "negA"
negB = fasta2counts("../Data/H2R_SABnegB_preprocessedQ25.fasta") 
colnames(negB)[2] = "negB"

## merge tables with Expected Seqs: use left join to only keep Expected Seqs' rows
library(tidyverse)
#dflist = list(ExpectedSeqs, inA, inB, negA, negB)
dflistAll = list(ExpectedSeqsAll, inA, inB, negA, negB)
#H2R = dflist %>% reduce(left_join, by='sequence')
H2RAll = dflistAll %>% reduce(left_join, by='sequence')
#head(H2R)
head(H2RAll)

## Add 1 pseudocount to all samples to avoid dividing by 0; first replace NA values with 0
H2RAll$inA[is.na(H2RAll$inA)] = 0
H2RAll$inB[is.na(H2RAll$inB)] = 0
H2RAll$negA[is.na(H2RAll$negA)] = 0
H2RAll$negB[is.na(H2RAll$negB)] = 0
H2RAll$inA = H2RAll$inA+1
H2RAll$inB = H2RAll$inB+1
H2RAll$negA = H2RAll$negA+1
H2RAll$negB = H2RAll$negB+1


## Normalize: Find sum of counts columns
#total_counts = apply(H2R[,c("inA", "inB","negA", "negB")], 2, sum, na.rm=TRUE) #apply sum fxn to all columns in this list
total_countsAll = apply(H2RAll[,c("inA", "inB","negA", "negB")], 2, sum, na.rm=TRUE) #apply sum fxn to all columns in this list
# convert total counts to millions
#total_cpm = total_counts/10^6
total_cpmAll = total_countsAll/10^6

#### Hm... I lost a lot of counts, and also the -TF sequences are basically not really there
sum(inA$inA) #2.4M
inA33 = inA[which(nchar(inA$sequence)==33),]
inA39 = inA[which(nchar(inA$sequence)==39),]
sum(inA33$inA) #126k - SHIT, THEY REALLY AREN'T THERE???
sum(inA39$inA) #1.3M
#total_counts #inA has $1.4M = 30k less than all of inA 33+39nt, not much filtered = reasonable
total_countsAll #2.1M = 88% of all inA reads; yup, these errors (seen in Sanger) account for ~all sequencing reads


## Add new columns normalized to total counts per million
#for (tempColname in names(total_cpm) ) {
#  newColName <- paste(tempColname,"_norm", sep="")
#  H2R[,newColName] <- H2R[,tempColname] / total_cpm[tempColname]
#}
for (tempColname in names(total_cpmAll) ) {
  newColName <- paste(tempColname,"_norm", sep="")
  H2RAll[,newColName] <- H2RAll[,tempColname] / total_cpmAll[tempColname]
}




## Calculate enrichment for each sort relative to input (normalized values)
#H2R[,"enrichA"] = H2R[,"negA_norm"] / H2R[,"inA_norm"]
#H2R[,"enrichB"] = H2R[,"negB_norm"] / H2R[,"inB_norm"]
#colnames(H2R)
#H2R[,"enrichAvg"] = rowMeans(H2R[21:22])

H2RAll[,"enrichA"] = H2RAll[,"negA_norm"] / H2RAll[,"inA_norm"]
H2RAll[,"enrichB"] = H2RAll[,"negB_norm"] / H2RAll[,"inB_norm"]
colnames(H2RAll)
H2RAll[,"enrichAvg"] = rowMeans(H2RAll[21:22])

## Calculate the standard deviation of enrichment scores
#H2R[,"enrichStdDev"] = apply(H2R[,21:22], 1, sd)
H2RAll[,"enrichStdDev"] = apply(H2RAll[,21:22], 1, sd)



### EXPORT to file so that I can re-analyze / have non-averaged values
#write.csv(H2R, file="../Analysis/H2R_SIVsab.csv",row.names = FALSE)
write.csv(H2RAll, file="../Analysis/H2Rall_SIVsabPseudocount.csv",row.names = FALSE)

# Save point: read back in
#H2R = read.csv("../Analysis/H2R_SIVsab.csv", sep=",")
H2RAll = read.csv("../Analysis/H2Rall_SIVsabPseudocount.csv", sep=",")



##### NEED TO FILTER OUT LOW-REPRESENTATION SEQUENCES, WHICH WILL BE NOISY

### A) Decide representation cutoffs
## Plot inputA vs. input B counts
par(mfcol=c(2,1))
#plot(H2R$inA_norm, H2R$inB_norm, log="xy")
plot(H2RAll$inA_norm, H2RAll$inB_norm, log="xy")
# mostly correlated, but definite outliers (likely not real counts) that have <~10cpm in BOTH libraries
plot(H2RAll$negA_norm, H2RAll$negB_norm, log="xy")
# well correlated above 10 cpm in each library

## Plot StdDev against input counts to choose cutoff
#plot(H2R$inA_norm, H2R$enrichStdDev, log="x")
#plot(H2R$inB_norm, H2R$enrichStdDev, log="x")
# 50cpm looks like a good cutoff, and is good for consistency's sake
plot(H2RAll$inA_norm, H2RAll$enrichStdDev, log="x")
plot(H2RAll$inB_norm, H2RAll$enrichStdDev, log="x")
#10cpm seems like it will do well

## Zoom in on the low count range to find cutoff for high std dev
#plot(H2R$inA_norm, H2R$enrichStdDev, xlim= c(1,100),ylim=c(1,20))
#plot(H2R$inB_norm, H2R$enrichStdDev, xlim= c(1,100), ylim=c(1,20))
plot(H2RAll$inA_norm, H2RAll$enrichStdDev, xlim= c(0,100),ylim=c(1,20))
plot(H2RAll$inB_norm, H2RAll$enrichStdDev, xlim= c(0,100), ylim=c(1,20))
#Yes, try 10cpm

## Plot enrichment scores against input values
#plot(H2R$inA_norm, H2R$enrichAvg, log="x")
#plot(H2R$inB_norm, H2R$enrichAvg, log="x")
#plot(H2R$inA_norm, H2R$enrichAvg, xlim= c(1,100))
#plot(H2R$inB_norm, H2R$enrichAvg, xlim= c(1,100))
plot(H2RAll$inA_norm, H2RAll$enrichAvg, log="x")
plot(H2RAll$inB_norm, H2RAll$enrichAvg, log="x")
plot(H2RAll$inA_norm, H2RAll$enrichAvg, xlim= c(1,100))
plot(H2RAll$inB_norm, H2RAll$enrichAvg, xlim= c(1,100))
#yup 10 in each input will nicely get rid of the noisiest points

plot(H2RAll$negA_norm, H2RAll$enrichStdDev, log="x")
plot(H2RAll$negB_norm, H2RAll$enrichStdDev, log="x")

### B) Filter out sequences based on these cutoffs
#H2R_cpmFilter = H2R[which(H2R$inA_norm > 50 & H2R$inB_norm > 50), ]
H2RAll_cpmFilter = H2RAll[which(H2RAll$inA_norm > 10 & H2RAll$inB_norm > 10), ]
#NB: this retains some of the -2aa insertion sequences, which are poorly represented but perhaps real?

#write.csv(H2R_cpmFilter, file="../Analysis/H2R_SIVsab_filter.csv", row.names = FALSE)
write.csv(H2RAll_cpmFilter, file="../Analysis/H2Rall_SIVsabPseudocount_filter10cpm.csv", row.names = FALSE)

