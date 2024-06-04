library(microDecon)

example <- cbind.data.frame(c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6"),
                            c(0,200,1000,50,0,25),
                            c(0,220,800,30,0,10),
                            c(0,180,1300,70,0,30),
                            c(60,660,1440,70,2400,30),
                            c(64,520,1000,48,1900,20),
                            c(40,480,700,35,2100,15),
                            c("K_Bacteria; P_Actinobacteria","K_Bacteria; P_Proteobacteria","K_Bacteria; P_Proteobacteria","K_Bacteria; P_Bacteroidetes","K_Bacteria","K_Bacteria"))
colnames(example) <- c("OTU_ID","Blank1","Blank2","Blank3","Pop1_Sample1","Pop1_Sample2","Pop2_Sample3","Taxa")


dat2 = data.frame(t(dat))
dat2$OTU_ID = c(paste0('OTU', 1:nrow(dat2)))
dat2$taxa = rownames(dat2)
rownames(dat2) = NULL
dat2 = dat2 %>%
  relocate(OTU_ID, .before = New_trimmed_1)

decontaminated = decon(dat2, numb.blanks = 6, numb.ind = c(90, 96), taxa = T)
