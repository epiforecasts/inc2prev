library(data.table)
library(tidyverse)

child_age_groups = c('2-10', '11-15')
estimates_c = readRDS('outputs/estimates_kids_combined_age_ab.rds')
estimates_a = readRDS('outputs/estimates_age_ab.rds')

estimates_c = estimates_c %>% filter(variable %in% child_age_groups)
estimates_a = estimates_a %>% filter(!(variable %in% child_age_groups))

estiamtes_all = rbind(estimates_c, estimates_a)

saveRDS(estiamtes_all, 'outputs/estimates_combined_age_ab_A.rds')


samples_c = readRDS('outputs/samples_kids_combined_age_ab.rds')
samples_a = readRDS('outputs/samples_age_ab.rds')

samples_c = samples_c %>% filter(variable %in% child_age_groups)
samples_a = samples_a %>% filter(!(variable %in% child_age_groups))

samples_all = rbind(samples_c, samples_a)

saveRDS(samples_all, 'outputs/samples_combined_age_ab_A.rds')
