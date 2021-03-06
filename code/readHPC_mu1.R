#---------------------------------------------------------

# Fully Sparse, Mu1

#--------------------------------------------------------




dir_fs_mu1 <- './data/cluster/results/fully_sparse/mu1/'
file_list <- list.files(dir_fs_mu1)


### mint 2 maxt 10


# n 10 

nm <- "n10_muf1_mint2_maxt10"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))




# n 20 

nm <- "n20_muf1_mint2_maxt10"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))



# n 30

nm <- "n30_muf1_mint2_maxt10"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))



# n 40

nm <- "n40_muf1_mint2_maxt10"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 50 

nm <- "n50_muf1_mint2_maxt10"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))









### mint 6 maxt 14


# n 10 

nm <- "n10_muf1_mint6_maxt14"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))




# n 20 

nm <- "n20_muf1_mint6_maxt14"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))



# n 30 

nm <- "n30_muf1_mint6_maxt14"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))



# n 40 

nm <- "n40_muf1_mint6_maxt14"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 50 

nm <- "n50_muf1_mint6_maxt14"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_fs_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_fs_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))

#------------------------------------------------------------------









#---------------------------------------------------------------

# Sparsified, Mu1

#---------------------------------------------------------------


# sparsity .06



dir_sf_mu1 <- './data/cluster/results/sparsified/mu1/'
file_list <- list.files(dir_sf_mu1)


### mint 2 maxt 10


# n 10 

nm <- "n10_muf1_sparsity0.06"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 20 

nm <- "n20_muf1_sparsity0.06"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 30 

nm <- "n30_muf1_sparsity0.06"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))




# n 40 

nm <- "n40_muf1_sparsity0.06"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 50 

nm <- "n50_muf1_sparsity0.06"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))







# sparsity .10

# n 10 

nm <- "n10_muf1_sparsity0.1"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 20 

nm <- "n20_muf1_sparsity0.1"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 30 

nm <- "n30_muf1_sparsity0.1"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))




# n 40 

nm <- "n40_muf1_sparsity0.1"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))





# n 50 

nm <- "n50_muf1_sparsity0.1"

sel <- grep(nm, file_list, value = T)

load(paste0(dir_sf_mu1, sel[1]))

Result <- c()
for (i in sel) {
  load(paste0(dir_sf_mu1, i))
  Result$SASE <- rbind(Result$SASE, result$SASE)
  Result$MC <- rbind(Result$MC, result$MC)
}

save(Result, file = paste0('./data/cluster/', nm, '.Rdata'))