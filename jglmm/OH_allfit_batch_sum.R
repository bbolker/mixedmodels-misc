library(tidyverse)
library(broom.mixed)

load("OH_allfit_batch.rda")
aa_sum <- (map_dfr(aa, tidy, .id="optimizer", effects="fixed", conf.int=TRUE)
    %>% select(optimizer, term, estimate, std.error, conf.low, conf.high)
)
saveRDS(aa_sum, file="OH_allfit_batch_sum.rds")
