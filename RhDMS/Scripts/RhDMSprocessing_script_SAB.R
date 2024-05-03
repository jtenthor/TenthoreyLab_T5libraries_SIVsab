# Script to process RhTRIM5 DMS data
# Experiment: SIVsab LOF


### PROCESS SAMPLE PRE-PROCESSED FASTA FILES INTO TABLE WITH READ COUNTS FOR EACH CONDITION

## load Biostrings module
library(Biostrings)

## define function to turn pre-processed fasta file into table
fasta2counts = function(file) {
  fa = readDNAStringSet(file) #read in fasta file: gives width (int), names (char), and seq (not char!)
  fa_headers_split = strsplit(names(fa),"-") #split the seq names into list of 2-element list (current format as ID#-counts)
  output = data.frame(sequence = as.character(fa), #make data frame with extracted DNA seqs
                      counts = as.integer(sapply(fa_headers_split,"[[",2))) #and counts = 2nd element w/in sublist (+ change to int class)
  return(output)
}

## run fxn on inputA library to get table
inA = fasta2counts("../Data/RD_SABinA_collapsed.fasta") 
head(inA)
colnames(inA)[2] = "inA"

## get table of all libraries for this experiment
inB = fasta2counts("../Data/RD_SABinB_collapsed.fasta")
colnames(inB)[2] = "inB"

posA = fasta2counts("../Data/RD_SABposA_collapsed.fasta")
colnames(posA)[2] = "posA"

posB = fasta2counts("../Data/RD_SABposB_collapsed.fasta")
colnames(posB)[2] = "posB"

plasmid = fasta2counts("../Data/RD_plasmid_collapsed.fasta")
colnames(plasmid)[2] = "plasmid"

## merge counts tables based on sequence
library(tidyverse)
dflist = list(plasmid, inA, inB, posA, posB)
RhDMS = dflist %>% reduce(full_join, by='sequence')
head(RhDMS)



### FILTER SEQUENCES BASED ON CODON CHANGES FROM WT

## FILTER 1. 39nt seqs only (expected length of RhT5v1)
RhDMS$ntLen = nchar(RhDMS$sequence)
RhDMS_filter1 = RhDMS[RhDMS$ntLen==39,]
#filters from 6.8M to 66k unique sequences

## translate sequences
AAs = translate(DNAStringSet(RhDMS_filter1$sequence), no.init.codon = TRUE) #note avoid Biostring forcing first codon to be a Met for ATG or CTG or TTG
RhDMS_filter1$aa_seq = as.character(AAs)

## Load seqinr package for the following functions
### don't do early on lest it mask Biostrings translate function
library(seqinr)

## Define function to count # of codon differences b/w given seq and WT
CountCodonDiffs2 = function(Testnt) {
  WTnt = "CAGGCACCAGGGACATTATTTACGTTTCCGTCACTCACG" #RsT5 WT v1, 39nt
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  CodonDiffs = sum(TestCodons != WTCodons) #finds number of differences from WT for each seq
  return(CodonDiffs)
}

## Apply fxn to whole table, make new column w/ # of codon differences
RhDMS_filter1[,"dist_codon"] = sapply(RhDMS_filter1[,"sequence"], CountCodonDiffs2) 

## FILTER 2. remove seqs with more than 1 codon difference
RhDMS_filter2 = RhDMS_filter1[which(RhDMS_filter1[,"dist_codon"] <= 1), ] #keep all rows that meet this criterion & all columns
#in this instance, down to 820 unique seqs, should really reduce processing time

## Define function to return mutated codon
ReturnMutCodon = function(Testnt) {
  WTnt = "CAGGCACCAGGGACATTATTTACGTTTCCGTCACTCACG" #RhT5 WT v1
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  FirstMutCodon = (TestCodons[FirstMutPosition]) #return this codon's sequence
  return(FirstMutCodon)
}

## Apply fxn to whole table, make new column w/ codon sequence at mutated site
RhDMS_filter2[,"mut_codon"] = sapply(RhDMS_filter2[,"sequence"], ReturnMutCodon)

## Make new column with just the wobble position
RhDMS_filter2[,"mut_codonWob"] = sapply(strsplit(RhDMS_filter2[,"mut_codon"],''), '[', 3)

## FILTER 3. Filter seqs where mutated codon doesn't end in C or G (non-NNS)
## NB: | indicates OR statements in R
RhDMS_filter3 = RhDMS_filter2[which(RhDMS_filter2[,"mut_codonWob"] == "C" | RhDMS_filter2[,"mut_codonWob"] == "G" | is.na(RhDMS_filter2[,"mut_codonWob"])), ]
# can check that this worked and returned only C/G (or NA)
unique(RhDMS_filter3[,"mut_codonWob"])
# now down to 411 unique seqs (expected max = 13 positions * 32 codons = 416 seqs)

## Define function to return mutated position (index) based on codon, not amino acid, comparison
ReturnMutIndex = function(Testnt) {
  WTnt = "CAGGCACCAGGGACATTATTTACGTTTCCGTCACTCACG" #RhT5 WT v1 nt seq; 39nt
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  return(FirstMutPosition)
}

## Apply fxn to whole table, make new column
RhDMS_filter3[,"mutcod_index"] = sapply(RhDMS_filter3[,"sequence"], ReturnMutIndex)
#NB: now have unique identifier for WT-synonymous variants!

## Define function to translate mutated codon
TranslateCodon = function(NtSequence) {
  if (is.na(NtSequence)) {
    ProteinSeq = NA #return null value for empty sets
  } 
  else {
    ProteinSeq = translate(s2c(NtSequence))
  }
  return(ProteinSeq)
}

## Apply function to whole table, make new column
RhDMS_filter3[,"mut_aa"] = sapply(RhDMS_filter3[,"mut_codon"], TranslateCodon)

## Now change index to actual amino acid position
RhDMS_filter3[,"mutcod_index"] = RhDMS_filter3[,"mutcod_index"] + 331


## FILTER 4. Must be represented at least 1 count in both input libraries (else calculations problematic)
RhDMS_filter4 = RhDMS_filter3[which(RhDMS_filter3[,"inA"] > 0 & RhDMS_filter3$inB > 0), ]
dim(RhDMS_filter4)
#still have 411 sequences
tail(RhDMS_filter4) #good representation even at tail end!



#### MAKE NORMALIZED COUNTS AND ENRICHMENT SCORES

## Find sum of columns
total_counts = apply(RhDMS_filter4[,c("plasmid","inA", "inB","posA", "posB")], 2, sum, na.rm=TRUE) #apply sum fxn to all columns in this list
# convert total counts to millions
total_cpm = total_counts/10^6

## Add new columns normalized to total counts per million
for (tempColname in names(total_cpm) ) {
  newColName <- paste(tempColname,"_norm", sep="")
  RhDMS_filter4[,newColName] <- RhDMS_filter4[,tempColname] / total_cpm[tempColname]
}

## Calculate enrichment for each sort relative to input (normalized values)
RhDMS_filter4[,"inAvPlasmid"] = RhDMS_filter4[,"inA_norm"] / RhDMS_filter4 [,"plasmid"]
RhDMS_filter4[,"inBvPlasmid"] = RhDMS_filter4[,"inB_norm"] / RhDMS_filter4 [,"plasmid"]
RhDMS_filter4[,"enrichA"] = RhDMS_filter4[,"posA_norm"] / RhDMS_filter4[,"inA_norm"]
RhDMS_filter4[,"enrichB"] = RhDMS_filter4[,"posB_norm"] / RhDMS_filter4[,"inB_norm"]
colnames(RhDMS_filter4)
RhDMS_filter4$invPlasmidAvg = rowMeans(RhDMS_filter4[19:20])
RhDMS_filter4[,"enrichAvg"] = rowMeans(RhDMS_filter4[21:22])

## Calculate the standard deviation of enrichment scores
RhDMS_filter4$invPlasmidStDev = apply(RhDMS_filter4[,19:20], 1, sd)
RhDMS_filter4[,"enricRhtdDev"] = apply(RhDMS_filter4[,21:22], 1, sd)



### EXPORT to file so that I can re-analyze / have non-averaged values
write.csv(RhDMS_filter4, file="../Analysis/RhDMS_SIVsab_filter4.csv",row.names = FALSE)

### Save point: read in my file again here if I'd like to change the cutoffs below
RhDMS_filter4 = read.delim("../Analysis/RhDMS_SIVsab_filter4.csv", sep=",")



##### NEED TO FILTER OUT LOW-REPRESENTATION SEQUENCES, WHICH WILL BE NOISY

### A) Decide representation cutoffs
## Plot inputA vs. input B counts ... not that helpful, they're highly correlated
plot(RhDMS_filter4$inA_norm, RhDMS_filter4$inB_norm, log="xy")

## Plot StdDev against input counts to choose cutoff
par(mfcol=c(2,1))
plot(RhDMS_filter4$inA_norm, RhDMS_filter4$enrichStdDev, log="x")
plot(RhDMS_filter4$inB_norm, RhDMS_filter4$enrichStdDev, log="x")
# Actually, generally pretty low Stdev, cutoff somewhere less than 500cpm

## Zoom in on the low count range to find cutoff for high std dev
plot(RhDMS_filter4$inA_norm, RhDMS_filter4$enrichStdDev, xlim= c(1,500))
plot(RhDMS_filter4$inB_norm, RhDMS_filter4$enrichStdDev, xlim= c(1,500))
# use >100cpm in each? Will only filter 2 points, I think

## Plot enrichment scores against input values
plot(RhDMS_filter4$inA_norm, RhDMS_filter4$enrichAvg, log="x")
plot(RhDMS_filter4$inB_norm, RhDMS_filter4$enrichAvg, log="x")
plot(RhDMS_filter4$inA_norm, RhDMS_filter4$enrichAvg, xlim= c(1,100))
plot(RhDMS_filter4$inB_norm, RhDMS_filter4$enrichAvg, xlim= c(1,100))
#no real pattern I see here, apparently this isn't a big problem


### B) Filter out sequences based on these cutoffs
RhDMS_filter5 = RhDMS_filter4[which(RhDMS_filter4$inA_norm > 100 & RhDMS_filter4$inB_norm > 100), ]
#now down to 410 seqs from 411



##### AVERAGE ACROSS ALL NT SEQS FOR EACH AA MUTATION

dim(RhDMS_filter5) # check how many unique sequence we have; I expect 32 x 13 = 416; I have 410
table(table(RhDMS_filter5[,"aa_seq"])) # determine how many aa seqs have 1/2/3/etc. unique nt seqs
# NOTE HERE THAT ALL WT SEQS ARE GROUPED INSTEAD OF KEPT SEPARATE, hence 21 WT nt seqs
length(unique(RhDMS_filter5[,"aa_seq"])) # determine # unique aa seqs I have; I expect 20x13 +1 = 261; I have 260 (missing 1 non-WT)

## Pull out WT seqs to calculate average and stdev of enrichment
RhDMS_filter5_WT = RhDMS_filter5[which(RhDMS_filter5$aa_seq == "QAPGTLFTFPSLT"),]
write.csv(RhDMS_filter5_WT, file="../Analysis/RhDMS_SIVsab_filter5_WT.csv",row.names = FALSE)

## Concatenate mutcodon index with aa seq to differentiate WT seqs
RhDMS_filter5[,"mutindex+seq"] = paste(RhDMS_filter5$mutcod_index, RhDMS_filter5$aa_seq, sep="_")
length(unique(RhDMS_filter5$`mutindex+seq`)) # determine # unique aa seqs I have; I expect 21 x 13 = 273; I have 272

## Split whole dataframe into list of dataframes for eqch unique aa seq
RhDMS_filter5Split <- split(RhDMS_filter5, RhDMS_filter5[,"mutindex+seq"]) #2nd argument = what to split based upon
head(RhDMS_filter5Split,3)

## Define function to average values; also keeps only useful columns
get_averages = function(datf){  #named function that will act on dataframe
  datf[,"num_nt_seqs"] = dim(datf)[1] #identify # unique nt seqs (rows) for each aa seq (subsetted dataframe)
  datf[,"inAvPlasmid"] = mean(datf[,"inAvPlasmid"]) #add column w/ mean (across synonymous seqs only)
  datf[,"inBvPlasmid"] = mean(datf[,"inBvPlasmid"])
  datf[,"enrichA_avg"] = mean(datf[,"enrichA"]) 
  datf[,"enrichB_avg"] = mean(datf[,"enrichB"])
  datf[,"inVplasmid_Avg_RepAndAA"] = mean(datf[,"invPlasmidAvg"]) #add column with the mean (across reps & synonymous codons)
  datf[,"inVplasmid_SD_RepAndAA"] = sd(c(datf[,"inAvPlasmid"],datf[,"inBvPlasmid"]))
  datf[,"enrich_Avg_RepAndAA"] = mean(datf[,"enrichAvg"]) #add column with the mean (across reps & synonymous codons)
  datf[,"enrich_SD_RepAndAA"] = sd(c(datf[,"enrichA"],datf[,"enrichB"]))
  columnsToKeep <- c("dist_codon","aa_seq","mutcod_index","mut_aa", "num_nt_seqs", "inAvPlasmid", "inBvPlasmid","enrichA_avg", "enrichB_avg", "inVplasmid_Avg_RepAndAA", "inVplasmid_SD_RepAndAA", "enrich_Avg_RepAndAA","enrich_SD_RepAndAA")
  datf <- datf[1, columnsToKeep] # overwrite dataframe (w/in fxn only) with the first row & only these columns
  return(datf)
}

## Perform this averaging function on all unique aa seqs
RhDMS_filter5SplitAvgs = lapply (RhDMS_filter5Split, get_averages)
head(RhDMS_filter5SplitAvgs,3)

## Bring all these values back into single table
RhDMS_filter5Avgs = do.call("rbind", RhDMS_filter5SplitAvgs) #row bind each of these rows
rownames(RhDMS_filter5Avgs) = NULL #gets rid of useless row names
head(RhDMS_filter5Avgs,3)

## SAVE OUTPUT (enrichment scores averaged across synonymous codons) AS CSV FILE
write.csv(RhDMS_filter5Avgs, file="../Analysis/RhDMS_SIVsab_filter5_avg.csv",row.names = FALSE)




##### GENERATE A POSITION X AMINO ACID MATRIX

##null values are a problem for populating the table
## replace single null value in index position (WT nt seq)
## Only Q332 has no synonymous codons -> put it there
RhDMS_filter5Avgs$mutcod_index[is.na(RhDMS_filter5Avgs$mutcod_index)] = 332
RhDMS_filter5Avgs$mut_aa[is.na(RhDMS_filter5Avgs$mut_aa)] = "Q"

##table populates based on 1/2/3/etc. indexing, won't start at 332 (has 331 empty rows)
##so change it back to indexing starting at 1
RhDMS_filter5Avgs$mutcod_index = RhDMS_filter5Avgs$mutcod_index - 331

## Make empty dataframe
RhDMSmatrix = data.frame()
RhDMSmatrixTransduce = data.frame()

## Populate dataframe
for (rowindex in 1:dim(RhDMS_filter5Avgs)[1]) { #iterate over the # of rows (first value of dimensions of table)
  RowNum = RhDMS_filter5Avgs[rowindex,"mutcod_index"] # define the row # as the index position
  ColNum = RhDMS_filter5Avgs[rowindex, "mut_aa"] # define the col # as the amino acid
  RhDMSmatrix[RowNum,ColNum] = RhDMS_filter5Avgs[rowindex, "enrich_Avg_RepAndAA"] #populate dataframe at row/col with the enrichment score from that row
}

for (rowindex in 1:dim(RhDMS_filter5Avgs)[1]) { #iterate over the # of rows (first value of dimensions of table)
  RowNum = RhDMS_filter5Avgs[rowindex,"mutcod_index"] # define the row # as the index position
  ColNum = RhDMS_filter5Avgs[rowindex, "mut_aa"] # define the col # as the amino acid
  RhDMSmatrixTransduce[RowNum,ColNum] = RhDMS_filter5Avgs[rowindex, "inVplasmid_Avg_RepAndAA"] #populate dataframe at row/col with the enrichment score from that row
}

## Add column with proper position #
RhDMSmatrix[,"Position"] = c("Q332","A333","P334","G335","T336","L337","F338", "T339","F340","P341","S342","L343","T344")
RhDMSmatrixTransduce[,"Position"] = c("Q332","A333","P334","G335","T336","L337","F338", "T339","F340","P341","S342","L343","T344")

## Reorder rows into my favorite order
RhDMSmatrixOrdered = RhDMSmatrix[,c("Position","G","P","A","I","L","V","F","W","Y","M","C","S","T","N","Q","H","K","R","D","E","*")]
RhDMSmatrixTransduceOrdered = RhDMSmatrixTransduce[,c("Position","G","P","A","I","L","V","F","W","Y","M","C","S","T","N","Q","H","K","R","D","E","*")]

## Save this matrix as csv file
write.csv(RhDMSmatrixOrdered, file="../Analysis/RhDMS_SIVsab_filter5_EnrichMatrix.csv",row.names = FALSE)
write.csv(RhDMSmatrixTransduceOrdered, file="../Analysis/RhDMS_SIVsab_filter5_TransduceMatrix.csv",row.names = FALSE)
