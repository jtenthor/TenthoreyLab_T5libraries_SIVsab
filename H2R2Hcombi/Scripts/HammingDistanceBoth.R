### Calculate Hamming distance for each mutant from (1) WT or (2) functional variants

H2R = read.csv("../Analysis/H2Rall_SIVsabPseudocount_filter10cpm.csv")
R2H = read.csv("../Analysis/R2Hall_SIVsab_filter_10in1out.csv")

library(DescTools)

#Test with WT sequence first
HsWT = "GaRgtRY-QTFV"
RhWT = "QaPgtLF+PSLT"

H2R$WT_HammingDist = sapply(H2R$aaSeq, function(x){StrDist(x,HsWT,method="hamming")})
unique(H2R$WT_HammingDist) #0-9 differences - looks good!
R2H$WT_HammingDist = sapply(R2H$aaSeq, function(x){StrDist(x,RhWT,method="hamming")})
unique(R2H$WT_HammingDist) #0-9 differences - looks good!

##Do this with a for loop for all gain-of-fxn (human) or retention-of-fxn (rhesus) sequences
H2Rgof = read.csv("../Analysis/H2R_GOF.csv")
colnames(H2Rgof) = NULL
H2Rgof = as.list(H2Rgof)
H2Rgof

for(AAseq in H2Rgof) {
  newColname = paste(AAseq,"HammingDist", sep="_")
  H2R[,newColname] = sapply(H2R$aaSeq, function(x){StrDist(x,AAseq,method="hamming")})
}

H2R_HD = H2R[which(H2R$Pos3 != "x"),]
H2R_HD2 = H2R_HD[which(H2R_HD$Pos12 != "x"),]
write.csv(H2R_HD2, file="../Analysis/H2Rall_SIVsabPseudocount_filter10cpm_nonFS_Hamming.csv")


R2Hrof = read.csv("../Analysis/R2H_ROF.csv")
R2Hrof
colnames(R2Hrof) = NULL
roflist = as.list(R2Hrof)
roflist

for(AAseq in roflist) {
  newColname = paste(AAseq,"HammingDist", sep="_")
  R2H[,newColname] = sapply(R2H$aaSeq, function(x){StrDist(x,AAseq,method="hamming")})
}

R2H_HD = R2H[which(R2H$Pos3 != "x"),]
write.csv(R2H_HD, file="../Analysis/R2Hall_SIVsab_10in1out_nonFS_Hamming.csv")
