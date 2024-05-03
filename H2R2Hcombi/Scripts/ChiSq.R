## Run chi squared test on H2R gof and R2H rof variants

H2R330 = matrix(c(5,12,172,160), ncol=2)
colnames(H2R330) = c("GOF","notGOF")
rownames(H2R330) = c("Human","Rhesus")
chisq.test(H2R330) #p = 0.12
chisq.test(H2R330)$expected #yes, it calculates this based on total frequencies (50.7% human)

H2R332 = matrix(c(1,16,175,157), ncol=2)
chisq.test(H2R332)$p.value

H2R335 = matrix(c(2,15,185,147), ncol=2)
chisq.test(H2R335)$p.value

H2R336 = matrix(c(0,17,163,169), ncol=2)
chisq.test(H2R336)$p.value

H2RTF = matrix(c(0,17,93,239), ncol=2)
chisq.test(H2RTF)

H2R337 = matrix(c(2,15,168,164), ncol=2)
chisq.test(H2R337)$p.value

H2R338 = matrix(c(6,11,159,173), ncol=2)
chisq.test(H2R338)$p.value

H2R339 = matrix(c(10,7,163,169), ncol=2)
chisq.test(H2R339)$p.value

H2R340 = matrix(c(6,11,197,135), ncol=2)
chisq.test(H2R340)$p.value

R2H332 = matrix(c(11,7,131,158), ncol=2)
colnames(R2H332) = c("ROF","LOF")
rownames(R2H332) = c("rhesus","human")
chisq.test(R2H332)
chisq.test(R2H332)$expected

R2H334 = matrix(c(10,8,140,149), ncol=2)
chisq.test(R2H334)$p.value

R2H337 = matrix(c(18,0,134,155), ncol=2)
chisq.test(R2H337)$p.value

R2H338 = matrix(c(18,0,145,144), ncol=2)
chisq.test(R2H338)$p.value

R2HTF = matrix(c(16,2,237,52), ncol=2)
chisq.test(R2HTF)$p.value

R2H341 = matrix(c(17,1,134,154), ncol=2)
chisq.test(R2H341)$p.value

R2H342 = matrix(c(13,5,150,139), ncol=2)
chisq.test(R2H342)$p.value

R2H343 = matrix(c(7,11,150,139), ncol=2)
chisq.test(R2H343)$p.value

R2H344 = matrix(c(6,12,158,131), ncol=2)
chisq.test(R2H344)$p.value

