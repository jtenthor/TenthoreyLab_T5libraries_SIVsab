# Script to process R2H combinatorial library sequencing data
# Experiment: R2H SIVsab GOF

## load Biostrings module
library(Biostrings)

#ExpectedSeqs = read.csv("../Data/CombiSeqList.csv") 
ExpectedSeqsAll = read.csv("../Data/CombiSeqListActual_R.csv")

## define function to turn pre-processed fasta file into table
fasta2counts = function(file) {
  fa = readDNAStringSet(file) #read in fasta file: gives width (int), names (char), and seq (not char!)
  fa_headers_split = strsplit(names(fa),"-") #split the seq names into list of 2-element list (current format as ID#-counts)
  output = data.frame(sequence = as.character(fa), #make data frame with extracted DNA seqs
                      counts = as.integer(sapply(fa_headers_split,"[[",2))) #and counts = 2nd element w/in sublist (+ change to int class)
  return(output)
}

## run fxn on each library to get table
inA = fasta2counts("../Data/R2H_SABinA_preprocessedQ25.fasta") 
colnames(inA)[2] = "inA"
inB = fasta2counts("../Data/R2H_SABinB_preprocessedQ25.fasta") 
colnames(inB)[2] = "inB"
negA = fasta2counts("../Data/R2H_SABposA_preprocessedQ25.fasta") 
colnames(negA)[2] = "posA"
negB = fasta2counts("../Data/R2H_SABposB_preprocessedQ25.fasta") 
colnames(negB)[2] = "posB"

## merge tables with Expected Seqs: use left join to only keep Expected Seqs' rows
library(tidyverse)
#dflist = list(ExpectedSeqs, inA, inB, negA, negB)
dflistAll = list(ExpectedSeqsAll, inA, inB, negA, negB)
#R2H = dflist %>% reduce(left_join, by='sequence')
R2HAll = dflistAll %>% reduce(left_join, by='sequence')
#head(R2H)
head(R2HAll)

## Normalize: Find sum of counts columns
#total_counts = apply(R2H[,c("inA", "inB","posA", "posB")], 2, sum, na.rm=TRUE) #apply sum fxn to all columns in this list
total_countsAll = apply(R2HAll[,c("inA", "inB","posA", "posB")], 2, sum, na.rm=TRUE) #apply sum fxn to all columns in this list
# convert total counts to millions
#total_cpm = total_counts/10^6
total_cpmAll = total_countsAll/10^6

#### Hm... I lost a lot of counts, and also the -TF sequences are basically not really there
sum(inA$inA) #2.1M
inA33 = inA[which(nchar(inA$sequence)==33),]
inA39 = inA[which(nchar(inA$sequence)==39),]
sum(inA33$inA) #92k - SHIT, THEY REALLY AREN'T THERE???
sum(inA39$inA) #1.6M
total_counts #inA has $2M = 20k less than all of inA 33+39nt, not much filtered = reasonable
total_countsAll #2.05M = 97% of all inA reads; yup, these errors (seen in Sanger) account for ~all sequencing reads

## Add new columns normalized to total counts per million
#for (tempColname in names(total_cpm) ) {
#  newColName <- paste(tempColname,"_norm", sep="")
#  R2H[,newColName] <- R2H[,tempColname] / total_cpm[tempColname]
#}
for (tempColname in names(total_cpmAll) ) {
  newColName <- paste(tempColname,"_norm", sep="")
  R2HAll[,newColName] <- R2HAll[,tempColname] / total_cpmAll[tempColname]
}




## Calculate enrichment for each sort relative to input (normalized values)
#R2H[,"enrichA"] = R2H[,"posA_norm"] / R2H[,"inA_norm"]
#R2H[,"enrichB"] = R2H[,"posB_norm"] / R2H[,"inB_norm"]
#colnames(R2H)
#R2H[,"enrichAvg"] = rowMeans(R2H[21:22])

R2HAll[,"enrichA"] = R2HAll[,"posA_norm"] / R2HAll[,"inA_norm"]
R2HAll[,"enrichB"] = R2HAll[,"posB_norm"] / R2HAll[,"inB_norm"]
colnames(R2HAll)
R2HAll[,"enrichAvg"] = rowMeans(R2HAll[21:22])

## Calculate the standard deviation of enrichment scores
#R2H[,"enrichStdDev"] = apply(R2H[,21:22], 1, sd)
R2HAll[,"enrichStdDev"] = apply(R2HAll[,21:22], 1, sd)



### EXPORT to file so that I can re-analyze / have non-averaged values
#write.csv(R2H, file="../Analysis/R2H_SIVsab.csv",row.names = FALSE)
write.csv(R2HAll, file="../Analysis/R2Hall_SIVsab.csv",row.names = FALSE)

# Save point: read back in
#R2H = read.csv("../Analysis/R2H_SIVsab.csv", sep=",")
R2HAll = read.csv("../Analysis/R2Hall_SIVsab.csv", sep=",")




##### NEED TO FILTER OUT LOW-REPRESENTATION SEQUENCES, WHICH WILL BE NOISY

### A) Decide representation cutoffs
## Plot inputA vs. input B counts
par(mfcol=c(2,1))
#plot(R2H$inA_norm, R2H$inB_norm, log="xy")
plot(R2HAll$inA_norm, R2HAll$inB_norm, log="xy")
# mostly correlated, but definite outliers (likely not real counts) that have <~10cpm in BOTH libraries


## Plot StdDev against input counts to choose cutoff
#plot(R2H$inA_norm, R2H$enrichStdDev, log="x")
#plot(R2H$inB_norm, R2H$enrichStdDev, log="x")
# 50cpm looks like a good cutoff, and is good for consistency's sake
plot(R2HAll$inA_norm, R2HAll$enrichStdDev, log="x")
plot(R2HAll$inB_norm, R2HAll$enrichStdDev, log="x")
#ditto, 50cpm

## Zoom in on the low count range to find cutoff for high std dev
#plot(R2H$inA_norm, R2H$enrichStdDev, xlim= c(1,100),ylim=c(1,60))
#plot(R2H$inB_norm, R2H$enrichStdDev, xlim= c(1,100), ylim=c(1,60))
plot(R2HAll$inA_norm, R2HAll$enrichStdDev, xlim= c(1,100),ylim=c(1,60))
plot(R2HAll$inB_norm, R2HAll$enrichStdDev, xlim= c(1,100), ylim=c(1,60))
#50 norm counts will still leave some large StDev - may need to be more conservative?

## Plot enrichment scores against input values
#plot(R2H$inA_norm, R2H$enrichAvg, log="x")
#plot(R2H$inB_norm, R2H$enrichAvg, log="x")
#plot(R2H$inA_norm, R2H$enrichAvg, xlim= c(1,100))
#plot(R2H$inB_norm, R2H$enrichAvg, xlim= c(1,100))
plot(R2HAll$inA_norm, R2HAll$enrichAvg, log="x")
plot(R2HAll$inB_norm, R2HAll$enrichAvg, log="x")
plot(R2HAll$inA_norm, R2HAll$enrichAvg, xlim= c(1,100))
plot(R2HAll$inB_norm, R2HAll$enrichAvg, xlim= c(1,100))
#yup 50 in each input will nicely get rid of the noisiest points (even 20 would)


### B) Filter out sequences based on these cutoffs
#R2H_cpmFilter = R2H[which(R2H$inA_norm > 50 & R2H$inB_norm > 50), ]
R2HAll_cpmFilter10 = R2HAll[which(R2HAll$inA_norm > 10 & R2HAll$inB_norm > 10), ]
R2HAll_cpmFilter50 = R2HAll[which(R2HAll$inA_norm > 50 & R2HAll$inB_norm > 50), ]
#NB: this retains some of the -2aa insertion sequences, which are poorly represented but perhaps real?
#now down to 448 seqs (or 438 for 50cpm filter) from 513 expected

#write.csv(R2H_cpmFilter, file="../Analysis/R2H_SIVsab_filter.csv", row.names = FALSE)
write.csv(R2HAll_cpmFilter10, file="../Analysis/R2Hall_SIVsab_filter10.csv", row.names = FALSE)

## Pull out controls to calculate average and stdev of enrichment
R2H_cpmFilter_controls = R2H_cpmFilter[which(R2H_cpmFilter$aaSeq == "QaPgtLF+PSLT"),]
R2H[which(R2H$aaSeq == "QaPgtLF+PSLT"),]
#I retained both of them



### 11/17/22 update:
## Since I'm analyzing depletion, effectively, the very low counts in the sorted samples can drive big changes (as denominator)
## Require that they're at least represnted >1cpm in each output file to get rid of the biggest SD of depletion


##12/17/22 update: just do the =>1 count in both reps for 10cpm filter
R2HAll_doublefilter1 = R2HAll_cpmFilter10[which(R2HAll_cpmFilter10$posA_norm >0 & R2HAll_cpmFilter10$posB_norm >0),]
# this gets down to 434 variants

write.csv(R2HAll_doublefilter1, file="../Analysis/R2Hall_SIVsab_filter_10in1out.csv", row.names = FALSE)


