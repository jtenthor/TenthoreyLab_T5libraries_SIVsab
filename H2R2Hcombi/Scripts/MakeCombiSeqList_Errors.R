#Script to generate possible sequence list for combinatorial libaries

## set up list object that has the types of sequence we see at each "site"
SeqsToCombine <- list()
SeqsToCombine[[1]] <- c("GGG","CAG") #HsG330, RhQ332 
SeqsToCombine[[2]] <- "GCA" #HsA331/RhA333
SeqsToCombine[[3]] <- c("CGA","CCA","CA") #HsR332, RhP334 (both synonymous recoded), include -1fs
SeqsToCombine[[4]] <- "GGGACA" #HsG332/T333; RhG334/T335
SeqsToCombine[[5]] <- c("CGA","CTA") #HsR335 (recode), RhL336 (recode)
SeqsToCombine[[6]] <- c("TAT","TTT") #HsY336 (recode), RhF338
SeqsToCombine[[7]] <- c("","ACGTTT") #+/- 2aa (RhT339/F340)
SeqsToCombine[[8]] <- c("CAG","CCG") #HsQ337 (recode), RhP341 (recode)
SeqsToCombine[[9]] = c("ACA","TCA") #HsT338 (recode), RhS342 (recode)
SeqsToCombine[[10]] = c("TTT","TTG") #HsF339 (recode), RhL343 (recode)
SeqsToCombine[[11]] = c("GTG","ACG","ACTG") #HsV340, RhT341, include error I made in 1 oligo order

## go through the list

# initialize an empty object
sequenceCombos <- ""
expectedNumCombos <- 1

for(i in 1:length(SeqsToCombine)) {
  thisBitOptions <- SeqsToCombine[[i]]
  howManyOptions <- length(thisBitOptions)
  expectedNumCombos <- expectedNumCombos * howManyOptions
  howManyCombosExistAlready <- length(sequenceCombos)
  cat("round",i,"repeating the current list",howManyOptions,"times and adding these choices",thisBitOptions,"\n")
  
  # multiply the existing combos
  sequenceCombos <- rep(sequenceCombos, howManyOptions)
  
  # repeat the current site options enough times to add it to each one of the strings we are building up
  bitsToAdd <- rep(thisBitOptions, each=howManyCombosExistAlready)
  
  # check we didn't mess up
  if (length(sequenceCombos) != length(bitsToAdd)) {
    stop("ERROR - something is wrong with the code\n")
  }
  sequenceCombos <- paste(sequenceCombos, bitsToAdd, sep="")
}

# take a look
sequenceCombos

length(unique(sequenceCombos))

# expected length
expectedNumCombos

# Turn into dataframe
CombiSeqList = data.frame(variant = 1:expectedNumCombos, sequence = sequenceCombos)



### Now add columns for amino acids

## set up list object that has the types of sequence we see at each "site"
AAsToCombine <- list()
AAsToCombine[[1]] <- c("G","Q") #HsG330, RhQ332 
AAsToCombine[[2]] <- "a" #HsA331/RhA333
AAsToCombine[[3]] <- c("R","P","x") #HsR332, RhP334 (both synonymous recoded), include -1fs
AAsToCombine[[4]] <- "gt" #HsG332/T333; RhG334/T335
AAsToCombine[[5]] <- c("R","L") #HsR335 (recode), RhL336 (recode)
AAsToCombine[[6]] <- c("Y","F") #HsY336 (recode), RhF338
AAsToCombine[[7]] <- c("-","+") #+/- 2aa (RhT339/F340)
AAsToCombine[[8]] <- c("Q","P") #HsQ337 (recode), RhP341 (recode)
AAsToCombine[[9]] = c("T","S") #HsT338 (recode), RhS342 (recode)
AAsToCombine[[10]] = c("F","L") #HsF339 (recode), RhL343 (recode)
AAsToCombine[[11]] = c("V","T","x") #HsV340, RhT341, include error I made in 1 oligo order

## go through the list

# initialize an empty object
AACombos <- ""
expectedNumAACombos <- 1

for(i in 1:length(AAsToCombine)) {
  thisBitOptions <- AAsToCombine[[i]]
  howManyOptions <- length(thisBitOptions)
  expectedNumAACombos <- expectedNumAACombos * howManyOptions
  howManyCombosExistAlready <- length(AACombos)
  cat("round",i,"repeating the current list",howManyOptions,"times and adding these choices",thisBitOptions,"\n")
  
  # multiply the existing combos
  AACombos <- rep(AACombos, howManyOptions)
  
  # repeat the current site options enough times to add it to each one of the strings we are building up
  bitsToAdd <- rep(thisBitOptions, each=howManyCombosExistAlready)
  
  # check we didn't mess up
  if (length(AACombos) != length(bitsToAdd)) {
    stop("ERROR - something is wrong with the code\n")
  }
  AACombos <- paste(AACombos, bitsToAdd, sep="")
}

library(Biostrings)

AACombos[1]
class(AACombos)[1]

CombiSeqList$aaSeq = AACombos

splitAAseqs = function(AAsequence) {
  AA_split = strsplit(AAsequence,"") #split the seq into individual characters
  output = data.frame(Pos1 = sapply(AA_split,"[[",1),
                      Pos3 = sapply(AA_split,"[[",3),
                      Pos6 = sapply(AA_split,"[[",6),
                      Pos7 = sapply(AA_split,"[[",7),
                      Pos8ins2aa = sapply(AA_split,"[[",8),
                      Pos9 = sapply(AA_split,"[[",9),
                      Pos10 = sapply(AA_split,"[[",10),
                      Pos11 = sapply(AA_split,"[[",11),
                      Pos12 = sapply(AA_split,"[[",12))
  return(output)
}

## run fxn on AAcombos
AAdiv = splitAAseqs(AACombos)

Combis = cbind(CombiSeqList,AAdiv)
head(Combis)
unique(nchar(Combis$sequence)) #now can have -1 or +1fs on 33 and 39 = 6 choices

write.csv(Combis, file="../Data/CombiSeqList_errors.csv", row.names = FALSE)
