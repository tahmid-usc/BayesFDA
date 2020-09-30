source("./code/functions_regular.R")

source("./code/rbf_regular.R")
# source("./code/laplace_regular.R")
# source("./code/matern52_regular.R")
# source("./code/matern32_regular.R")



yeast <- read.delim("data/yeast.txt")
gene <- read_excel("data/gene.xlsx")
names(gene)[1] <- "X"
mdt <- merge(yeast, gene, by = "X")
phase <- ifelse(mdt$Peak == 'G1', 'G1', 'NonG1')
mdt <- cbind(mdt, phase)
mdt <- mdt[,-(2:7)]
mdt <- mdt[, -(20:78)]
mdt <- na.omit(mdt)


train <- split(mdt, mdt$phase, drop = T)

x <- seq(0,119, length.out = 18)

plot(x, train$G1[1,2:19], type = 'b', ylim = c(min(train$G1[,2:19]), max(train$G1[,2:19])), col = 'gray30', cex = .5,
     xlab = 'Time (mins)', ylab = 'Expression', cex.lab = 1.5)
apply(train$G1[,2:19], 1, function(y) {lines(x, y, type = 'b', col = 'gray30', cex = .5)})

plot(x, train$G1[1,2:19], type = 'b', ylim = c(min(train$G1[,2:19]), max(train$NonG1[,2:19])), col = 'gray30', cex = .5,
     xlab = 'Time (mins)', ylab = 'Expression', cex.lab = 1.5)
apply(train$NonG1[,2:19], 1, function(y) {lines(x, y, type = 'b', col = 'gray30', cex = .5)})



G1 <- train$G1[,-c(1,20)]

gpFit.g1 <- feature(G1)

y <- seq(0,1, length.out = 100)
postDist.g1 <- postDist(y, gpFit.g1)

y <- seq(0,119, length.out = 100)


#plot
plot(x, G1[1,], type = 'b', ylim = c(min(G1), max(G1)), col = rgb(0,0,0,.1), cex = .5,
     xlab = 'Time (mins)', ylab = 'Expression', cex.lab = 1.5)
apply(G1, 1, function(y) {lines(x, y, type = 'b', col = rgb(0,0,0,.1), cex = .5)})
lines(y, postDist.g1$mu, lwd = 4, col = 2)
# lines(y, postDist.g1$ul, lwd = 4, lty = 2, col = 2)
# lines(y, postDist.g1$ll, lwd = 4, lty = 2, col = 2)
polygon(c(y,rev(y)), c(postDist.g1$ll,rev(postDist.g1$ul)), col= rgb(1,0,0,.2),border=NA)



nonG1 <- train$NonG1[,-c(1,20)]

gpFit.nong1 <- feature(nonG1)

y <- seq(0,1, length.out = 100)
postDist.nong1 <- postDist(y, gpFit.nong1)

y <- seq(0,119, length.out = 100)


#plot
plot(x, nonG1[1,], type = 'b', ylim = c(min(nonG1), max(nonG1)), col = rgb(0,0,0,.1), 
     cex = .5, xlab = 'Time (mins)', ylab = 'Expression', cex.lab = 1.5)
apply(nonG1, 1, function(y) {lines(x, y, type = 'b', col = rgb(0,0,0,.1), cex = .5)})
lines(y, postDist.nong1$mu, lwd = 4, col = 2)
# lines(y, postDist.nong1$ul, lwd = 2, lty = 2)
# lines(y, postDist.nong1$ll, lwd = 2, lty = 2)
polygon(c(y,rev(y)), c(postDist.nong1$ll,rev(postDist.nong1$ul)), col= rgb(1,0,0,.2),border=NA)



#save(list = c('gpFit.g1', 'gpFit.nong1'), file = './data/gene/m32.Rdata')

