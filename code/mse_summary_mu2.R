library(dplyr)

dir <- './data/cluster/aggregate_mu2/'

#---------------------------

#Fully Sparse, mu2

#---------------------------


# min 2 max 10



mse_df <- c()

fn <- 'n10_muf2_mint2_maxt10'
load(paste0(dir, fn, '.Rdata'))

data.frame(t(apply(Result$SASE, 2, mean)))


mse_df <- rbind(mse_df, 
                data.frame(n = 10, mint = 2, maxt = 10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))




fn <- 'n20_muf2_mint2_maxt10'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 20, mint = 2, maxt = 10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n30_muf2_mint2_maxt10'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 30, mint = 2, maxt = 10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))



fn <- 'n40_muf2_mint2_maxt10'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 40, mint = 2, maxt = 10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n50_muf2_mint2_maxt10'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 50, mint = 2, maxt = 10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))





# min 6 max 14


fn <- 'n10_muf2_mint6_maxt14'
load(paste0(dir, fn, '.Rdata'))

data.frame(t(apply(Result$SASE, 2, mean)))


mse_df <- rbind(mse_df, 
                data.frame(n = 10, mint = 6, maxt = 14) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))




fn <- 'n20_muf2_mint6_maxt14'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 20, mint = 6, maxt = 14) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n30_muf2_mint6_maxt14'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 30, mint = 6, maxt = 14) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))



fn <- 'n40_muf2_mint6_maxt14'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 40, mint = 6, maxt = 14) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n50_muf2_mint6_maxt14'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 50, mint = 6, maxt = 14) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


save(mse_df, file = './data/cluster/mse_fs_mu2.Rdata')


















#---------------------------

# Sparsified, mu2

#---------------------------


# sparsity .06



mse_df <- c()

fn <- 'n10_muf2_sparsity0.06'
load(paste0(dir, fn, '.Rdata'))

data.frame(t(apply(Result$SASE, 2, mean)))


mse_df <- rbind(mse_df, 
                data.frame(n = 10, sparsity = 0.06) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))




fn <- 'n20_muf2_sparsity0.06'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 20, sparsity = 0.06) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n30_muf2_sparsity0.06'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 30, sparsity = 0.06) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))



fn <- 'n40_muf2_sparsity0.06'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 40, sparsity = 0.06) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n50_muf2_sparsity0.06'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 50, sparsity = 0.06) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))





# sparsity 0.10


fn <- 'n10_muf2_sparsity0.1'
load(paste0(dir, fn, '.Rdata'))

data.frame(t(apply(Result$SASE, 2, mean)))


mse_df <- rbind(mse_df, 
                data.frame(n = 10, sparsity = 0.10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))




fn <- 'n20_muf2_sparsity0.1'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 20, sparsity = 0.10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n30_muf2_sparsity0.1'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 30, sparsity = 0.10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))



fn <- 'n40_muf2_sparsity0.1'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 40, sparsity = 0.10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))


fn <- 'n50_muf2_sparsity0.1'
load(paste0(dir, fn, '.Rdata'))

mse_df <- rbind(mse_df, 
                data.frame(n = 50, sparsity = 0.10) 
                %>% cbind(data.frame(t(apply(Result$SASE, 2, mean)))))




save(mse_df, file = './data/cluster/mse_sf_mu2.Rdata')
