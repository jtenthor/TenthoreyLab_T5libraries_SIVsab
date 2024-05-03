# Script to process HsIndel data
# Experiment: SIVsab GOF

## load Biostrings module
library(Biostrings)

ExpectedSeqs = read.csv("../Data/HI_ExpectedSeqs.csv") 
ExpectedSeqs$sequence = toupper(ExpectedSeqs$sequence) #matching function is case-sensitive

## define function to turn pre-processed fasta file into table
fasta2counts = function(file) {
  fa = readDNAStringSet(file) #read in fasta file: gives width (int), names (char), and seq (not char!)
  fa_headers_split = strsplit(names(fa),"-") #split the seq names into list of 2-element list (current format as ID#-counts)
  output = data.frame(sequence = as.character(fa), #make data frame with extracted DNA seqs
                      counts = as.integer(sapply(fa_headers_split,"[[",2))) #and counts = 2nd element w/in sublist (+ change to int class)
  return(output)
}

## run fxn on each library to get table
plasmid = fasta2counts("../Data/HI_plasmid_collapsed.fasta")
colnames(plasmid)[2]="plasmid"
inA = fasta2counts("../Data/HI_SABinA_collapsed.fasta") 
colnames(inA)[2] = "inA"
inB = fasta2counts("../Data/HI_SABinC_collapsed.fasta") 
colnames(inB)[2] = "inB"
negA = fasta2counts("../Data/HI_SABnegA_collapsed.fasta") 
colnames(negA)[2] = "negA"
negB = fasta2counts("../Data/HI_SABnegC_collapsed.fasta") 
colnames(negB)[2] = "negB"

## merge tables with Expected Seqs: use left join to only keep Expected Seqs' rows
library(tidyverse)
dflist = list(ExpectedSeqs, plasmid, inA, inB, negA, negB)
HsIndel = dflist %>% reduce(left_join, by='sequence')
head(HsIndel)

## Avoid 0 counts in some samples to avoid dividing by 0 (or getting NA results)
min(HsIndel$inA) #NA
min(HsIndel$inB) #NA
min(HsIndel$negA) #NA
min(HsIndel$negB) #NA
## Add 1 pseudocount to all libraries; first replace NA values with 0
HsIndel$inA[is.na(HsIndel$inA)] = 0
HsIndel$inB[is.na(HsIndel$inB)] = 0
HsIndel$negA[is.na(HsIndel$negA)] = 0
HsIndel$negB[is.na(HsIndel$negB)] = 0
HsIndel$inA = HsIndel$inA+1
HsIndel$inB = HsIndel$inB+1
HsIndel$negA = HsIndel$negA+1
HsIndel$negB = HsIndel$negB+1

## Normalize: Find sum of counts columns
total_counts = apply(HsIndel[,c("plasmid","inA", "inB","negA", "negB")], 2, sum, na.rm=TRUE) #apply sum fxn to all columns in this list
# convert total counts to millions
total_cpm = total_counts/10^6

## Add new columns normalized to total counts per million
for (tempColname in names(total_cpm) ) {
  newColName <- paste(tempColname,"_norm", sep="")
  HsIndel[,newColName] <- HsIndel[,tempColname] / total_cpm[tempColname]
}

## find out how many reads were culled in left_join
sum(inA[,2]) #278,093 (this was a low read library) 
total_counts[2] #236,672 (85% retained)
sum(is.na(HsIndel$inA)) #5, now 0 by def
sum(is.na(HsIndel$inB)) #3, now 0 by def = mostly represented, some missing oh well

## Calculate enrichment for each sort relative to input (normalized values)
HsIndel[,"inAvPlasmid"] = HsIndel[,"inA_norm"] / HsIndel [,"plasmid"]
HsIndel[,"inBvPlasmid"] = HsIndel[,"inB_norm"] / HsIndel [,"plasmid"]
HsIndel[,"enrichA"] = HsIndel[,"negA_norm"] / HsIndel[,"inA_norm"]
HsIndel[,"enrichB"] = HsIndel[,"negB_norm"] / HsIndel[,"inB_norm"]
colnames(HsIndel)
HsIndel$invPlasmidAvg = rowMeans(HsIndel[17:18])
HsIndel[,"enrichAvg"] = rowMeans(HsIndel[19:20])

## Calculate the standard deviation of enrichment scores
HsIndel$invPlasmidStDev = apply(HsIndel[,17:18], 1, sd)
HsIndel[,"enrichStdDev"] = apply(HsIndel[,19:20], 1, sd)



### EXPORT to file so that I can re-analyze / have non-averaged values
write.csv(HsIndel, file="../Analysis/HsIndel_SIVsabPseudocount.csv",row.names = FALSE)

# Save point: read back in
HsIndel = read.csv("../Analysis/HsIndel_SIVsabPseudocount.csv", sep=",")




##### NEED TO FILTER OUT LOW-REPRESENTATION SEQUENCES, WHICH WILL BE NOISY

### A) Decide representation cutoffs
## Plot inputA vs. input B counts ... not that helpful, they're highly correlated
plot(HsIndel$inA_norm, HsIndel$inB_norm, log="xy")

## Plot StdDev against input counts to choose cutoff
par(mfcol=c(2,1))
plot(HsIndel$inA_norm, HsIndel$enrichStdDev, log="x")
plot(HsIndel$inB_norm, HsIndel$enrichStdDev, log="x")
# 4 high SD data points  will be easily filtered by >10cpm in both

## Zoom in on the low count range to find cutoff for high std dev
plot(HsIndel$inA_norm, HsIndel$enrichStdDev, xlim= c(1,50),ylim=c(1,50))
plot(HsIndel$inB_norm, HsIndel$enrichStdDev, xlim= c(1,50), ylim=c(1,50))
# 10cpm will leave one high SD data point and filter the rest

## Plot enrichment scores against input values
plot(HsIndel$inA_norm, HsIndel$enrichAvg, log="x")
plot(HsIndel$inB_norm, HsIndel$enrichAvg, log="x")
plot(HsIndel$inA_norm, HsIndel$enrichAvg, xlim= c(1,100))
plot(HsIndel$inB_norm, HsIndel$enrichAvg, xlim= c(1,100))
#no real pattern here


### B) Filter out sequences based on these cutoffs
HsIndel_cpmFilter = HsIndel[which(HsIndel$inA_norm > 10 & HsIndel$inB_norm > 10), ]
#now down to 402 (instead of 398 seqs at 50cpm filter) from 411
write.csv(HsIndel_cpmFilter, file="../Analysis/HsIndel_SIVsabPseudocount_filter10.csv", row.names = FALSE)

## Pull out controls to calculate average and stdev of enrichment
HsIndel_cpmFilter_controls = HsIndel_cpmFilter[which(HsIndel_cpmFilter$mutation == "WT" | HsIndel_cpmFilter$mutation == "syn" | HsIndel_cpmFilter$mutation == "stop"),]
write.csv(HsIndel_cpmFilter_controls, file="../Analysis/HsIndel_SIVsabPseudocount_filter10_controls.csv",row.names = FALSE)



#### Make a matrix of averages
## Make empty dataframe
HsIndelMatrix = data.frame()
HsIndelTransduceMatrix = data.frame()

##change the WT (not ordered) to position 324
HsIndel_cpmFilter[which(HsIndel_cpmFilter$position==0),]$position = 324
##table populates based on 1/2/3/etc. indexing, won't start at 324
##so change it back to indexing starting at 1
HsIndel_cpmFilter$position = HsIndel_cpmFilter$position - 323

## Populate dataframe
for (rowindex in 1:dim(HsIndel_cpmFilter)[1]) { #iterate over the # of rows (first value of dimensions of table)
  RowNum = HsIndel_cpmFilter[rowindex,"position"] # define the row # as the position of indel
  ColNum = HsIndel_cpmFilter[rowindex, "mutation"] # define the col # as the indel mutation
  HsIndelMatrix[RowNum,ColNum] = HsIndel_cpmFilter[rowindex, "enrichAvg"] #populate dataframe at row/col with the enrichment score from that row
}
for (rowindex in 1:dim(HsIndel_cpmFilter)[1]) { #iterate over the # of rows (first value of dimensions of table)
  RowNum = HsIndel_cpmFilter[rowindex,"position"] # define the row # as the position of indel
  ColNum = HsIndel_cpmFilter[rowindex, "mutation"] # define the col # as the indel mutation
  HsIndelTransduceMatrix[RowNum,ColNum] = HsIndel_cpmFilter[rowindex, "invPlasmidAvg"] #populate dataframe at row/col with the transduction enrichment score from that row
}

## Add column with proper position #
HsIndelMatrix[,"Position"] = c("K324","P325","Q326","I327","I328","Y329","G330","A331","R332","G333","T334","R335","Y336","Q337","T338","F339","V340","N341","F342","N343","Y344","C345")
HsIndelTransduceMatrix[,"Position"] = c("K324","P325","Q326","I327","I328","Y329","G330","A331","R332","G333","T334","R335","Y336","Q337","T338","F339","V340","N341","F342","N343","Y344","C345")

## Reorder rows into my favorite order
HsIndelMatrixOrdered = HsIndelMatrix[,c("Position","syn","stop","del1","del2","del3","del4","del5","del6","del7","del8","del9","dupe1","dupe2","dupe3","dupe4","dupe5","dupe6","dupe7","dupe8","dupe9","WT")]
HsIndelTransduceMatrixOrdered = HsIndelTransduceMatrix[,c("Position","syn","stop","del1","del2","del3","del4","del5","del6","del7","del8","del9","dupe1","dupe2","dupe3","dupe4","dupe5","dupe6","dupe7","dupe8","dupe9","WT")]

## Save this matrix as csv file
write.csv(HsIndelMatrixOrdered, file="../Analysis/HsIndel_SIVsabPseudocount_EnrichMatrix.csv",row.names = FALSE)
write.csv(HsIndelTransduceMatrixOrdered, file="../Analysis/HsIndel_SIVsabPseudocount_TransduceMatrix.csv",row.names = FALSE)
