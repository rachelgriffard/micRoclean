library(tidyverse)

# data
blocklist = read.csv("contaminant-blocklist.csv", header = F)
gen = read.csv("Level6_Genus.csv", header=T,row.name=1)

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

index = grep(".*.g__*", colnames(gen)) # keep if have IDed genus
gen = gen[,index]
colnames(gen) = sub(".*.g__", "", colnames(gen)) # subset to only genus
index2 = which(names(gen)=="")
gen = gen[,-index2] # remove if no value

dat = gen

# subset
set.seed(42)

features = colnames(gen)
sub = sample(features, 282) # randomly select half

dat = select(gen, all_of(sub))


########################### PERFect ###############################################
# devtools::install_github("cxquy91/PERFect")

# Attempt to find statistic defining diversity of community 

library(PERFect)

# "This algorithm's objective is to remove noise taxa while identifying and 
# retaining signal taxa. It takes an OTU table and a test critical value alpha
# as inputs and produces a reduced OTU table with less taxa."
res_sim = PERFect_sim(X = dat)
dim(res_sim$filtX)      
colnames(res_sim$filtX) # signal taxa
p = pvals_Plots(PERFect = res_sim, X = dat, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

library(ggplot2)
pl = p$plot + ggtitle("Simultanenous Filtering")

library(plotly)
ggplotly(pl)

DFL = data.frame("DFL" = res_sim$DFL[order(res_sim$DFL)])

# more robust than simulation filtering as above, according to package authors
# res_perm = PERFect_perm(X = gen, Order = "pvals", pvals_sim = res_sim, algorithm = "full")
# res_perm2 = PERFect_perm(X = dat, Order = "pvals", pvals_sim = res_sim, algorithm = "fast")
# p1 = pvals_Plots(res_perm, gen)
# p1 = p1$plot + ggtitle("Full Algorithm")
# p2 = pvals_Plots(res_perm2, dat)
#p2 = p2$plot + ggtitle("Fast Algorithm")
# p2
# ggpubr::ggarrange(p1,p2,ncol = 2)

########################### SCRuB ###############################################
# devtools::install_github("shenhav-and-korem-labs/SCRuB")
library(SCRuB)

control = group
control = control == "Negative Control"

sample = group
sample[!sample == "Control"] = "Plasma"

dat = as.matrix(dat)

meta = data.frame("is_control" = control,
                  "sample" = sample)
rownames(meta) = rownames(dat)

scr_out = SCRuB(dat, 
                metadata = meta)
sc_dat = data.frame(scr_out$decontaminated_samples)

sc_sim = PERFect_sim(X = sc_dat)
dim(sc_sim$filtX)      
colnames(sc_sim$filtX) # signal taxa
p = pvals_Plots(PERFect = sc_sim, X = sc_dat, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

library(ggplot2)
sc_pl = p$plot + ggtitle("Post-SCRuB Simultanenous Filtering")
ggplotly(sc_pl)

cowplot::plot_grid(pl, sc_pl)

sc_DFL = data.frame("DFL" = sc_sim$DFL[order(sc_sim$DFL)])

library(dplyr)
library(tibble)
DFL = rownames_to_column(DFL, "feature")
sc_DFL = rownames_to_column(sc_DFL, "feature")
DFL_com = DFL %>%
  left_join(sc_DFL, by = c('feature' = 'feature'))

DFL_pl = DFL_com %>%
  pivot_longer(!feature, names_to = 'Group', values_to = 'DFL') %>%
  mutate(Group = ifelse(Group == 'DFL.x', 'Pre-SCRuB', 'Post-SCRuB')) %>%
  mutate(DFL = ifelse(is.na(DFL)==TRUE, 0, DFL)) %>%
  mutate(order(DFL)) %>%
  ggplot(aes(fct_reorder(feature, DFL), DFL, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge2') +
  xlab('Feature') + ggtitle('DFL Comparison Pre- and Post-SCRuB Method') +
  coord_flip()
ggplotly(DFL_pl)

DFL_pl_sub = DFL_com %>%
  pivot_longer(!feature, names_to = 'Group', values_to = 'DFL') %>%
  filter(DFL > 1e-3) %>%
  mutate(Group = ifelse(Group == 'DFL.x', 'Pre-SCRuB', 'Post-SCRuB')) %>%
  mutate(DFL = ifelse(is.na(DFL)==TRUE, 0, DFL)) %>%
  mutate(order(DFL)) %>%
  ggplot(aes(fct_reorder(feature, DFL), DFL, fill = Group)) +
  geom_bar(stat = 'identity', position = 'dodge2') +
  xlab('Feature') + ggtitle('DFL Comparison Pre- and Post-SCRuB Method') +
  coord_flip()
ggplotly(DFL_pl_sub)

features = DFL[,1]

library(randomcoloR)
color = data.frame('feature' = features,
                   'color' = randomColor(count = length(features),
                                         hue = 'random'))
DFL = merge(DFL, color, by = 'feature')
sc_DFL = merge(sc_DFL, color, by = 'feature')

#focus on highest DFL values, those above 1e-3
DFL_pl = DFL %>%
  filter(DFL > 1e-3) %>%
  ggplot(aes(fct_reorder(feature, DFL), DFL)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  ggtitle('Pre-Scrub') + xlab('Feature') +
  theme_bw() +
  theme(legend.position = 'None')
sc_pl_DFL = sc_DFL %>%
  filter(DFL > 1e-3) %>%
  ggplot(aes(fct_reorder(feature, DFL), DFL)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  ggtitle('Post-Scrub') + xlab('Feature') +
  theme_bw() +
  theme(legend.position = 'None')
cowplot::plot_grid(DFL_pl, sc_pl_DFL)

########################### microDecon ###############################################
# devtools::install_github("donaldtmcknight/microDecon")
library(microDecon)

sum(control)

ind = which(control)

dat = t(dat)
dat = tibble::rownames_to_column(dat, "OTU ID")

decontaminated = decon(data = dat,numb.blanks=12,numb.ind=ind)
