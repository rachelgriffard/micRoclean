# Flesh out well2well function
# Rachel Griffard
# 20240312

library(tidyverse)
library(kohonen) # SOM
library(aweSOM)

# load  data from pipeline.R file within same directory

set.seed(42)

# compare within batches, requires user to have put in metadata regarding batches

rownames(meta) = rownames(gen)

head(meta)

w2w = dat %>%
  mutate(dat[order(row.names(dat))]) %>%
  as.matrix() %>%
  scale()

som_grid = somgrid(xdim = 3, ydim = 3, topo = 'hexagonal')

som_mod = som(w2w, grid=som_grid, rlen=1000, alpha=c(1, .01), keep.data=T)

plot(som_mod,
     type = 'mapping')

plot(som_mod, type='changes')

plot(som_mod, type='count')

plot(som_mod, type='dist.neighbours')

plot(som_mod, type='codes')

aweSOMplot(som_mod, type = 'Cloud')

aweSOMplot(som_mod, type = 'UMatrix')

aweSOMsmoothdist(som_mod)
