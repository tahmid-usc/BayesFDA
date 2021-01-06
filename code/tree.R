# tree data

library(dplyr)
library(ggplot2)


tree <- read.csv("E:/R Project/BayesFDA/data/tree/tree.csv")
tree <- tree %>% dplyr::select(Site, ID.Code, Tree, Rep, Sp, Year, Height) %>% filter(!is.na(Height)) %>% 
  mutate(id = paste0(Site, ID.Code, Rep), logHeight = log(Height))
View(tree)

summary(tree)

table(tree$Site)
table(tree$Sp)
table(tree$Site, tree$Sp, dnn = c('Site', 'Species'))

tree %>% group_by(id) %>% summarize(n())

ggplot(data = tree, aes(x = Year, y = Height, col = Sp, group = id)) +
  geom_line(size = 1.5, alpha = .2) + facet_wrap(~Site)


ggplot(data = tree, aes(x = Year, y = Height, col = Sp, group = id)) +
  geom_line(size = 1.5, alpha = .5) 

#split data on species

splitData <- split(tree, tree$Sp)

# train GP model

library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)

source("code/functions.R")
source("code/RBF.R")

gp_Sp <- lapply(splitData, function(L) {feature(id = L$id, x = L$Year, y = L$Height)})
gp_Sp_log <- lapply(splitData, function(L) {feature(id = L$id, x = L$Year, y = L$logHeight)})

# gp_BS <- feature(id = splitData$BS$id, x = splitData$BS$Year, y = splitData$BS$Height)  
# gp_EL <- feature(id = splitData$EL$id, x = splitData$EL$Year, y = splitData$EL$logHeight)  
# gp_HL <- feature(id = splitData$HL$id, x = splitData$HL$Year, y = splitData$HL$logHeight)  
# gp_JL <- feature(id = splitData$JL$id, x = splitData$JL$Year, y = splitData$JL$logHeight)  
# gp_JP <- feature(id = splitData$JP$id, x = splitData$JP$Year, y = splitData$JP$logHeight)  
# gp_RP <- feature(id = splitData$RP$id, x = splitData$RP$Year, y = splitData$RP$logHeight)  
# gp_TL <- feature(id = splitData$TL$id, x = splitData$TL$Year, y = splitData$TL$logHeight)  
#  





testTime <- seq(0, 30, length.out = 100)

# postDist_BS <- postDist(testTime, gp_Sp$BS)
# postDist_EL <- postDist(testTime, gp_Sp$EL)
# postDist_HL <- postDist(testTime, gp_Sp$HL)
# postDist_JL <- postDist(testTime, gp_Sp$JL)
# postDist_JP <- postDist(testTime, gp_Sp$JP)
# postDist_RP <- postDist(testTime, gp_Sp$BS)
# postDist_TL <- postDist(testTime, gp_Sp$BS)


postDist_Sp <- lapply(gp_Sp, function(L){postDist(testTime, L)})
postDist_Sp_log <- lapply(gp_Sp_log, function(L){postDist(testTime, L)})


plotGP(gp_Sp$BS, postDist_Sp$BS, trans = .05)
plotGP(gp_Sp$EL, postDist_Sp$EL, trans = .05)
plotGP(gp_Sp$HL, postDist_Sp$HL, trans = .05)
plotGP(gp_Sp$JL, postDist_Sp$JL, trans = .05)
plotGP(gp_Sp$JP, postDist_Sp$JP, trans = .05)
plotGP(gp_Sp$RP, postDist_Sp$RP, trans = .05)
plotGP(gp_Sp$TL, postDist_Sp$TL, trans = .05)
plotGP(gp_Sp$WS, postDist_Sp$WS, trans = .05)




plotGP(gp_Sp_log$BS, postDist_Sp_log$BS, trans = .05)
plotGP(gp_Sp_log$EL, postDist_Sp_log$EL, trans = .05)
plotGP(gp_Sp_log$HL, postDist_Sp_log$HL, trans = .05)
plotGP(gp_Sp_log$JL, postDist_Sp_log$JL, trans = .05)
plotGP(gp_Sp_log$JP, postDist_Sp_log$JP, trans = .05)
plotGP(gp_Sp_log$RP, postDist_Sp_log$RP, trans = .05)
plotGP(gp_Sp_log$TL, postDist_Sp_log$TL, trans = .05)
plotGP(gp_Sp_log$WS, postDist_Sp_log$WS, trans = .05)
