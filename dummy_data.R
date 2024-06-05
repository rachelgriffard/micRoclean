### dummy data for development
dat = read.csv("Level6_Genus.csv", header=T,row.name=1)
batch = dat$Batch
group = dat$Groups

index = grep(".*.g__*", colnames(dat)) # keep if have IDed datus
dat = dat[,index]
comp = data.frame(colnames(dat))
comp$compare = sub(".*.g__", "", colnames(dat)) # subset to only datus
index2 = which(comp$compare=="")
comp = comp[-index2,]

index = grep(".*.g__*", colnames(dat)) # keep if have IDed datus
dat = dat[,index]
colnames(dat) = sub(".*.g__", "", colnames(dat)) # subset to only datus
index2 = which(names(dat)=="")
dat = dat[,-index2] # remove if no value

## dummy res data for pipeline2
res = read.csv('res_toydata.csv')

### dummy data for FL function
removed = sample(rownames(FL), 100)

### dummy meta data
control = group
control = control == "Negative Control"

sample = group
sample[!sample == "Control"] = "Plasma"

dat = as.matrix(dat)

meta = data.frame("is_control" = control,
                  "sample" = sample,
                  "batch" = batch)
rownames(meta) = rownames(dat)

### dummy technical replicates (p2s3)
technical_replicates = data.frame("Batch1" = c("Old_trimmed_2", "Old_trimmed_86",
                             "Old_trimmed_85", "Old_trimmed_49",
                             "Old_trimmed_38", "Old_trimmed_3",
                             "Old_trimmed_13", "Old_trimmed_26"),
                "Batch2" = c("New_trimmed_29", "New_trimmed_35",
                             "New_trimmed_41", "New_trimmed_47",
                             "New_trimmed_53", "New_trimmed_59",
                             "New_trimmed_65", "New_trimmed_71"))
