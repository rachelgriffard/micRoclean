### dummy data for development
dat = read.csv("Level6_datus.csv", header=T,row.name=1)
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

res = read.csv('res_toydata.csv')

### dummy data for FL function
removed = sample(rownames(FL), 100)