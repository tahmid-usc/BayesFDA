# summarize mse Mu1
library(dplyr)
library(xtable)

load("./data/cluster/mse_fs_mu1.Rdata")
xtable(mse_df[, -(2:3)])

load("./data/cluster/mse_sf_mu1.Rdata")
xtable(mse_df[, -2])






# summarize mse Mu2

load("./data/cluster/mse_fs_mu2.Rdata")
xtable(mse_df[, -(2:3)])

load("./data/cluster/mse_sf_mu2.Rdata")
xtable(mse_df[, -2])