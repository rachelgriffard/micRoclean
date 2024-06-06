# MB low biomass pipeline
# Adapted from Zozoya et al. 2021
# Adapted by Rachel Griffard

# session >> set working directory >> to source file location

### Require older version of ANCOMBC due to removal of pairwise comparison of newer ancombc() function
# Use R 4.1.1
# BiocManager::install("ANCOMBC")

# load libraries
library("phyloseq") # object wrapper for pipeline
library("tidyverse") # stringr, dplyr, ggplot2, tidyr
library("reshape2")
library("scales")
library("decontam") # method for step 3 of pipeline
library("metagenomeSeq")
library("metagMisc")
library("knitr")
library("irr")
library("vegan")
library("Biostrings")
library("compositions")
library("factoextra")
library("ANCOMBC") # appropriate for compositional data as in microbiome studies diff abundance analysis
library("mia")
library("microbiome")
library("dendextend")
library("biomeUtils") # for comparing phyloseq objects
library("DT")
library("ggVennDiagram") # for comparing those removed

dir.create("Results") # create directory for Results

# Zozoya et al. 2021 functions
low.count.removal = function(
    data, # OTU count data frame of size p (OTU) x n (sample); (rows x columns)
    percent=0.01 # cutoff chosen
){
  OTU_percent_abund = rowSums(data)*100/(sum(rowSums(data)))
  keep.otu = which(OTU_percent_abund > percent)
  data.filter = data[keep.otu,]
  return(list(data.filter = data.filter, OTU_percent_abund = OTU_percent_abund[keep.otu]))
}

## Modification of the low.count.removal function that filters based on a range of abundances
abund.range.filter = function(
    data, # OTU count data frame of size p (OTU) x n (sample); (rows x columns)
    min_percent=0.01,
    max_percent=100
){
  OTU_percent_abund = rowSums(data)*100/(sum(rowSums(data)))
  keep.otu = which(OTU_percent_abund > min_percent & OTU_percent_abund < max_percent)
  data.filter = data[keep.otu,]
  return(list(data.filter = data.filter, OTU_percent_abund = OTU_percent_abund[keep.otu]))
}

# Import blocklist (Eisenhofer et al., 2019) and genus data
blocklist = read.csv("contaminant-blocklist.csv", header = F)
gen = read.csv("Level6_Genus.csv", header=T,row.name=1)

# Extract meta data
batch = gen$Batch
group = gen$Groups
age = gen$Age
race = gen$Race
gen = gen[,1:(ncol(gen)-5)] # remove non-count cols

index = grep(".*.g__*", colnames(gen)) # keep if have IDed genus
gen = gen[,index]
comp = data.frame(colnames(gen))
comp$compare = sub(".*.g__", "", colnames(gen)) # subset to only genus
index2 = which(comp$compare=="")
comp = comp[-index2,]

# Adjust names to genus alone
index = grep(".*.g__*", colnames(gen)) # keep if have IDed genus
gen = gen[,index]
colnames(gen) = sub(".*.g__", "", colnames(gen)) # subset to only genus
index2 = which(names(gen)=="")
gen = gen[,-index2] # remove if no value
## removes 650 based on lack of genus alone

# Create genus phyloseq object with all meta data
d.genus = cbind(batch, group, age, race, group, gen) # remove non MB, move groups to front

colnames(d.genus)[1] = "Batch"
d.genus$Batch[d.genus$Batch=="Old"] = "1. Old" # to set old as reference
d.genus$Batch[d.genus$Batch=="New"] = "2. New"
table(d.genus$Batch)

colnames(d.genus)[2] = "Groups"
d.genus$Groups[d.genus$Groups=="Negative Control"] = "0.Negative Control"
d.genus$Groups[d.genus$Groups=="Control"] = "1.Control"
d.genus$Groups[d.genus$Groups=="Benign"] = "2.Sample"
d.genus$Groups[d.genus$Groups=="Non_EOC"] = "2.Sample"
d.genus$Groups[d.genus$Groups=="EOC"] = "2.Sample"
table(d.genus$Groups)

colnames(d.genus)[3] = "Age"
table(d.genus$Age)

colnames(d.genus)[4] = "Race"
table(d.genus$Race)

colnames(d.genus)[5] = "Ex_Groups"
d.genus$Groups[d.genus$Groups=="Negative Control"] = "0.Negative Control"
d.genus$Groups[d.genus$Groups=="Control"] = "1.Control"
d.genus$Groups[d.genus$Groups=="Benign"] = "2.Sample"
d.genus$Groups[d.genus$Groups=="Non_EOC"] = "3.Sample"
d.genus$Groups[d.genus$Groups=="EOC"] = "4.Sample"
table(d.genus$Groups)

d.genus.ancomBC = t(d.genus[,6:ncol(d.genus)]) # non meta data
d.genus.ancomBC = round(d.genus.ancomBC)

all(colnames(d.genus.ancomBC)==rownames(d.genus))
otu_mat = d.genus.ancomBC
meta = data.frame(group=d.genus$Groups,
                  batch=d.genus$Batch,
                  # age=d.genus$Age,
                  # race=d.genus$Race,
                  ex_groups=d.genus$Ex_Groups,
                  row.names=colnames(d.genus.ancomBC))
tax_mat = matrix(rownames(otu_mat),nrow=nrow(otu_mat),ncol=1)
rownames(tax_mat) = rownames(otu_mat)
colnames(tax_mat) = c("genus")

OTU_g = otu_table(otu_mat, taxa_are_rows = TRUE)
META_g = sample_data(meta)
TAX_g = tax_table(tax_mat)
genus = phyloseq(OTU_g, META_g, TAX_g)

####################################################################################
# Identify abundance categories within samples
# Get OTUs that fall into the low-abundance category (below 0.1%) 
sum(taxa_sums(genus))*0.001 # OTU abundance max-threshold in absolute numbers
genus_low = genus
df = abund.range.filter(otu_table(genus), min_percent = 0, max_percent = 0.1)
otu_table(genus_low) = otu_table(df$data.filter)
ntaxa(genus_low)  # Number of OTUs in the "low abund category"
(ntaxa(genus_low) / ntaxa(genus)) * 100 

# Get OTUs that fall into the medium-abundance category (0.1% - 1%) 
sum(taxa_sums(genus))*0.01 # OTU abundance max-threshold in absolute numbers
genus_med = genus
df = abund.range.filter(otu_table(genus), min_percent = 0.1, max_percent = 1.0)
otu_table(genus_med) = otu_table(df$data.filter)
ntaxa(genus_med)  # Number of OTUs in the "medium abund category"
(ntaxa(genus_med) / ntaxa(genus)) * 100 

# Get OTUs that fall into the high-abundance category (> 1%) 
sum(taxa_sums(genus))*0.1 # OTU abundance max-threshold in absolute numbers
genus_high = genus
df = abund.range.filter(otu_table(genus), min_percent = 1.0)
otu_table(genus_high) = otu_table(df$data.filter)
ntaxa(genus_high)  # Number of OTUs in the "high abund category"
(ntaxa(genus_high) / ntaxa(genus)) * 100 

####################################################################################
### Step 1: Remove features that showed different abundance in different batches ###
####################################################################################
# Equivalence test more appropriate?

# run ANCOMBC
# https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
genus_s1_res = ancombc(phyloseq = genus, assay_name = "counts", 
                       group = "batch", p_adj_method = "BH",  lib_cut = 0,
                       formula = "batch", 
                       struc_zero = TRUE, neg_lb = FALSE,
                       tol = 1e-5, max_iter = 100, conserve = FALSE,
                       alpha = 0.05, global = TRUE) 

# make phyloseq into TreeSummarizedExperiment (approp for ANCOMBC2)
#tse1 = mia::makeTreeSummarizedExperimentFromPhyloseq(genus)

# run ANCOMBC2 - not peer-reviewed yet
#genus_s1_res = ancombc2(data=tse1, assay_name = "counts", tax_level = "genus",
#                        fix_formula = "batch", rand_formula = NULL,
#                        p_adj_method = "BH", group = "batch", 
#                        alpha = 0.05,
#                        global = TRUE,
#                        iter_control = list(tol = 1e-5, max_iter = 20, 
#                                            verbose = FALSE),
#                        em_control = list(tol = 1e-5, max_iter = 100))

# genus_s11_res = genus_s1_res$res %>%
#   mutate_if(is.numeric, function(x) round(x, 2))
# table(genus_s11_res$`diff_batch2. New`)

genus_s11_res = do.call(cbind, genus_s1_res$res)

s1_rem = rownames(genus_s11_res[genus_s11_res$`diff_abn.batch2. New`==TRUE,])

index3 = grep(paste(s1_rem,collapse="$|"), colnames(gen))
gen_s1 = gen[,-index3]

# recreate phyloseq object without batch 1
d.genus_s1 = cbind(batch, group, age, race, gen_s1) # remove non MB, move groups to front

d.genus_s1 = d.genus_s1[d.genus_s1$batch=="New",] # subset to only batch 2 (neg controls)
table(d.genus_s1$batch)

colnames(d.genus_s1)[2] = "Groups"
d.genus_s1$Groups[d.genus_s1$Groups=="Negative Control"] = "0.Negative Control"
d.genus_s1$Groups[d.genus_s1$Groups=="Control"] = "1.Control"
d.genus_s1$Groups[d.genus_s1$Groups=="Benign"] = "2.Benign"
d.genus_s1$Groups[d.genus_s1$Groups=="Non_EOC"] = "3.Non_EOC"
d.genus_s1$Groups[d.genus_s1$Groups=="EOC"] = "4.EOC"
table(d.genus_s1$Groups)

colnames(d.genus_s1)[3] = "Age"
table(d.genus_s1$Age)

colnames(d.genus_s1)[4] = "Race"
table(d.genus_s1$Race)

d.genus_s1.ancomBC = t(d.genus_s1[,5:ncol(d.genus_s1)]) # non meta data
d.genus_s1.ancomBC = round(d.genus_s1.ancomBC)

all(colnames(d.genus_s1.ancomBC)==rownames(d.genus_s1))
otu_mat = d.genus_s1.ancomBC
meta = data.frame(group=d.genus_s1$Groups,
                  age=d.genus_s1$Age,
                  race=d.genus_s1$Race,
                  row.names=colnames(d.genus_s1.ancomBC))
tax_mat = matrix(rownames(otu_mat),nrow=nrow(otu_mat),ncol=1)
rownames(tax_mat) = rownames(otu_mat)
colnames(tax_mat) = c("genus")

OTU_g = otu_table(otu_mat, taxa_are_rows = TRUE)
META_g = sample_data(meta)
TAX_g = tax_table(tax_mat)
genus_s1 = phyloseq(OTU_g, META_g, TAX_g)

# compare original to subset after step 1
ps.list = c("Original" = genus,
            "Step 1" = genus_s1)
biomeUtils::comparePhyloseq(ps.list)

####################################################################################
#### Step 2: Remove features that are differentially abund in negative controls ####
####################################################################################
d.genus_s2 = cbind(batch, group, age, race, gen) # remove non-MB, move groups to front

d.genus_s2 = d.genus_s2[d.genus_s2$batch=="New",] # subset to only batch 2 (has neg controls)
table(d.genus_s2$batch)

colnames(d.genus_s2)[2] = "Groups"
d.genus_s2$Groups[d.genus_s2$Groups=="Negative Control"] = "0.Negative Control"
d.genus_s2$Groups[d.genus_s2$Groups=="Control"] = "1.Control"
d.genus_s2$Groups[d.genus_s2$Groups=="Benign"] = "2.Sample"
d.genus_s2$Groups[d.genus_s2$Groups=="Non_EOC"] = "2.Sample"
d.genus_s2$Groups[d.genus_s2$Groups=="EOC"] = "2.Sample"
table(d.genus_s2$Groups)

colnames(d.genus_s2)[3] = "Age"
table(d.genus_s2$Age)

colnames(d.genus_s2)[4] = "Race"
table(d.genus_s2$Race)

d.genus_s2.ancomBC = t(d.genus_s2[,5:ncol(d.genus_s2)]) # non meta data
d.genus_s2.ancomBC = round(d.genus_s2.ancomBC)

all(colnames(d.genus_s2.ancomBC)==rownames(d.genus_s2))
otu_mat = d.genus_s2.ancomBC
meta = data.frame(group=d.genus_s2$Groups,
                  age=d.genus_s2$Age,
                  race=d.genus_s2$Race,
                  row.names=colnames(d.genus_s2.ancomBC))
tax_mat = matrix(rownames(otu_mat),nrow=nrow(otu_mat),ncol=1)
rownames(tax_mat) = rownames(otu_mat)
colnames(tax_mat) = c("genus")

OTU_g = otu_table(otu_mat, taxa_are_rows = TRUE)
META_g = sample_data(meta)
TAX_g = tax_table(tax_mat)
genus_s2 = phyloseq(OTU_g, META_g, TAX_g)

sample_data(genus_s2)$is.neg = rep(NA)
sample_data(genus_s2)$is.neg = sample_data(genus_s2)$group == "0.Negative Control"
table(sample_data(genus_s2)$is.neg)

decontam_output = isContaminant(genus_s2, method="prevalence", neg="is.neg", threshold=0.5)
table(decontam_output$contaminant)
index4 = which(decontam_output$contaminant)

# Prevalence of contaminant plot
ps.pa = transform_sample_counts(genus_s2, function(abund) 1*(abund>0))
ps.pa.neg = prune_samples(sample_data(ps.pa)$group == "0.Negative Control", ps.pa)
ps.pa.pos = prune_samples(sample_data(ps.pa)$group != "0.Negative Control", ps.pa)

df.pa = data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=decontam_output$contaminant)
library(viridis)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) +
  geom_point(alpha=0.3, size=.8) +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  theme_bw() +
  # scale_color_viridis_d("Contaminant", option="cividis") +
  ggtitle("Prevalence of contaminants",
          subtitle = "Negative controls v positive samples") +
  scale_color_manual("Contaminant", values=c("red", "darkblue"))
ggsave("./Results/ContaminantPrevPlot.png", width = 6, height = 4)

# get removal list
s2_rem = rownames(decontam_output[decontam_output$contaminant==TRUE,])


genus_s2 = prune_taxa(!decontam_output$contaminant, genus_s2)


# compare original, step 1, and subset of step 2
ps.list2 = c("Original" = genus,
             "Step 1" = genus_s1,
             "Step 2" = genus_s2)
biomeUtils::comparePhyloseq(ps.list2)

####################################################################################
########### Step 3: Remove if DA in diff batches for technical replicates ##########
####################################################################################

# identify resequenced samples
rs = data.frame("Batch1" = c("Old_trimmed_2", "Old_trimmed_86",
                              "Old_trimmed_85", "Old_trimmed_49",
                              "Old_trimmed_38", "Old_trimmed_3",
                              "Old_trimmed_13", "Old_trimmed_26"),
                "Batch2" = c("New_trimmed_29", "New_trimmed_35",
                              "New_trimmed_41", "New_trimmed_47",
                              "New_trimmed_53", "New_trimmed_59",
                              "New_trimmed_65", "New_trimmed_71"))

# tab = otu_table(genus)
# class(tab) = "matrix"
# dat = vegan::rarecurve(t(tab), step=50, cex=0.5, label = F, tidy = T)

# ind = grep("New*", dat$Site)
# dat$NO = ifelse(grepl("New*", dat$Site), "New", "Old")
# ggplot(dat, aes(Sample, Species, color = NO)) +
# # geom_line(alpha = 0.3, linewidth = 0.5) +
#  geom_line(alpha = 0.5, linewidth = 0.2, aes(color = Site)) +
#  theme_bw() +
#  xlim(0,1500) +
#  theme(legend.position = "none") +
#  ggtitle("Rarefaction curves")

# rarefy (debatable)
# set.seed(711) # reproducibility
# genus_rf = rarefy_even_depth(genus, rngseed = 711, replace=FALSE, trimOTUs = TRUE)

# create dataframes for shared samples
genus_new_sh = prune_samples(rs$Batch1, genus)
genus_old_sh = prune_samples(rs$Batch2, genus)

genus_new.df = otu_table(genus_new_sh) # Extract OTU table from phyloseq object
genus_new.df.PA = transform_sample_counts(genus_new.df, function(abund) 1*(abund>0)) # Set as present/absence

genus_old.df = otu_table(genus_old_sh) # Extract OTU table from phyloseq object
genus_old.df.PA = transform_sample_counts(genus_old.df, function(abund) 1*(abund>0))

# calculate kappa-statistic for all ASVs shared btwn DNA extraction batches
kappa_results = data.frame(value = numeric(), statistic = numeric(), p.value = numeric())
for(i in 1:nrow(genus_new.df.PA)){
  k = kappa2(t(rbind(genus_new.df.PA[i,], genus_old.df.PA[i,])), "unweighted")
  kappa_results[i,"value"] = k$value
  kappa_results[i,"statistic"] = k$statistic
  kappa_results[i,"p.value"] = k$p.value
}
row.names(kappa_results) = taxa_names(genus_new.df.PA)

kappa_results.no_NA = subset(kappa_results, !is.na(value) & !is.na(p.value))
kappa_results_sig = subset(kappa_results.no_NA, p.value < 0.05)
kappa_res_keep = subset(kappa_results.no_NA, p.value < 0.05 & value > 0.4)
kappa_res_keep = cbind(tax_table(genus)[row.names(kappa_res_keep), "genus"], kappa_res_keep)
kappa_res_keep[order(kappa_res_keep$value, decreasing = TRUE),]
kappa_res_remove = subset(kappa_results.no_NA, p.value >= 0.05 | value <= 0.4)

s3_rem = rownames(kappa_res_remove)

allTaxa = taxa_names(genus)
keepTaxa = allTaxa[!(allTaxa %in% rownames(kappa_res_remove))]
genus_s3 = prune_taxa(keepTaxa, genus)

####################################################################################
############ Step 4: Remove previously published 'contaminants-blocklist' ##########
####################################################################################
blocklist = as.vector(blocklist$V1)
allTaxa = taxa_names(genus)

s4_rem = allTaxa[(allTaxa %in% blocklist)]

myTaxa = allTaxa[!(allTaxa %in% blocklist)]
genus_s4 = prune_taxa(myTaxa, genus)

# compare original, step 1, and subset of step 2
ps.list4 = c("Original" = genus,
             "Step 1" = genus_s1,
             "Step 2" = genus_s2,
             "Step 3" = genus_s3,
             "Step 4" = genus_s4)
biomeUtils::comparePhyloseq(ps.list4)

####################################################################################
#################### Prune failed taxa from final phyloseq object ##################
####################################################################################
# compare those removed
x = list(s1_rem, s2_rem, s3_rem, s4_rem)
ggVennDiagram(x, stroke.size=1,
              category.names = c("Step 1", "Step 2", "Step 3", "Step 4"),
              edge_lty="solid", set_size=6,
              label_alpha = 0.5, label_percent_digit = 1) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  ggplot2::scale_color_grey(start=0, end=0) +
  scale_fill_distiller(direction=1) +
  labs(title="Taxa Removal by Step") +
  theme(legend.position="none", plot.title=element_text(size=25, hjust = 0.5)) +
  scale_fill_distiller(palette = "Spectral")
ggsave("./Results/VennTaxaRem.png", height = 6, width = 8)

# remove all that failed test
rem = unique(c(s1_rem, s2_rem, s3_rem, s4_rem))
df = data.frame("Genus"=rem, "Step1"=rep(NA),
                "Step2"=rep(NA), "Step3"=rep(NA),
                "Step4"=rep(NA))
df$Step1 = ifelse(df$Genus %in% s1_rem, "Yes", "No")
df$Step2 = ifelse(df$Genus %in% s2_rem, "Yes", "No")
df$Step3 = ifelse(df$Genus %in% s3_rem, "Yes", "No")
df$Step4 = ifelse(df$Genus %in% s4_rem, "Yes", "No")
write.csv(df, "./Results/RemovalTaxa.csv", row.names=F)

index_fin = grep(paste(rem,collapse="$|"), colnames(gen))
keep = colnames(gen[,-index_fin])

genus_final = prune_taxa(keep, genus)

ps.list5 = c("Original" = genus,
             "Step 1" = genus_s1,
             "Step 2" = genus_s2,
             "Step 3" = genus_s3,
             "Step 4" = genus_s4,
             "Final" = genus_final)
biomeUtils::comparePhyloseq(ps.list5)
res = taxa_names(genus_final)
length(res)
write.csv(res, "./Results/HighConfidenceTaxa.csv", row.names=F)

# find extra taxa - resolved with addition of $ to end of grep statement 
taxa = colnames(gen)
total = c(rem, res)
diff = taxa[!(taxa %in% total)]
diff



####################################################################################
########################### Run diff abdn btwn groups ##############################
####################################################################################
# remove batch 1 samples from genus_final object
# remove negative controls from genus_final object
neg.controls = c("New_trimmed_30",
                  "New_trimmed_36",
                  "New_trimmed_42",
                  "New_trimmed_48",
                  "New_trimmed_54",
                  "New_trimmed_60",
                  "New_trimmed_66",
                  "New_trimmed_72",
                  "New_trimmed_77",
                  "New_trimmed_82",
                  "New_trimmed_87",
                  "New_trimmed_92"
)
remove = c(neg.controls, rownames(gen)[97:192]) # combine removal of negative controls and old samples
samples = rownames(gen)
pos.samples = samples[!(samples %in% remove)]
genus_final = prune_samples(pos.samples, genus_final)
ps.list6 = c("Original" = genus,
             "Final" = genus_final)
comparePhyloseq(ps.list6)

# create d.genus object
d.genus_final = cbind(batch, group, gen)

d.genus_final = d.genus_final[d.genus_final$batch=="New",] # subset to only batch 2 (neg controls)
table(d.genus_final$batch)

colnames(d.genus_final)[2] = "Groups"
d.genus_final = subset(d.genus_final, Groups!="Negative Control") # remove negative controls
d.genus_final$Groups[d.genus_final$Groups=="Control"] = "1.Control"
d.genus_final$Groups[d.genus_final$Groups=="Benign"] = "2.Benign"
d.genus_final$Groups[d.genus_final$Groups=="Non_EOC"] = "3.Non_EOC"
d.genus_final$Groups[d.genus_final$Groups=="EOC"] = "4.EOC"
table(d.genus_final$Groups)

Groups = d.genus_final$Groups
d.genus_final = d.genus_final[colnames(d.genus_final) %in% keep]

d.genus_final.ancomBC = t(d.genus_final[,3:ncol(d.genus_final)]) # non meta data
d.genus_final.ancomBC = round(d.genus_final.ancomBC)

all(colnames(d.genus_final.ancomBC)==rownames(d.genus_final))
otu_mat = d.genus_final.ancomBC
meta = data.frame(group=Groups,
                  row.names=colnames(d.genus_final.ancomBC))
tax_mat = matrix(rownames(otu_mat),nrow=nrow(otu_mat),ncol=1)
rownames(tax_mat) = rownames(otu_mat)
colnames(tax_mat) = c("genus")

write.csv(otu_mat, "./DC_Deliverable/otu.csv")
OTU_final = otu_table(otu_mat, taxa_are_rows = TRUE)
write.csv(meta, "./DC_Deliverable/meta.csv")
META_final = sample_data(meta)
write.csv(tax_mat, "./DC_Deliverable/tax.csv")
TAX_final = tax_table(tax_mat)
genus_final = phyloseq(OTU_final, META_final, TAX_final)

ps.list7 = c("Original" = genus,
             "Final" = genus_final)
comparePhyloseq(ps.list7)

# run ancombc2 comparing across "sample"
# tse2 = mia::makeTreeSummarizedExperimentFromPhyloseq(genus_final)

# run ANCOMBC2
#genus_res = ancombc2(data=tse2, assay_name = "counts", tax_level = "genus",
#                        fix_formula = "group", rand_formula = NULL,
#                        p_adj_method = "BH", group = "group", 
#                        alpha = 0.05,
#                        global = TRUE,
#                        iter_control = list(tol = 1e-5, max_iter = 20, 
#                                            verbose = FALSE),
#                        em_control = list(tol = 1e-5, max_iter = 100))

genus_res_ancombc = ancombc(phyloseq = genus_final, formula = "group",   # Covariates can be added here on top of group
                            p_adj_method = "BH", lib_cut = 0,
                            group = "group", struc_zero = TRUE, neg_lb = FALSE, # group includes main groups to compare
                            tol = 1e-5, max_iter = 100, conserve = FALSE,
                            alpha = 0.05, global = TRUE)
d.ancombc.result = do.call(cbind,genus_res_ancombc$res)

d.ancombc.result = d.ancombc.result[order(d.ancombc.result$p_val.group2.Sample),]

write.csv(d.ancombc.result, "ancombc_result.csv")



