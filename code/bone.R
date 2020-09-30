rm(list=ls())

source("code/functions.R")


bone <- read_csv("data/spnbmd.csv")
bone <- as.data.frame(bone)


# data without discarding subjects with single observation

q <- table(bone$idnum)
q <- row.names(q)[as.numeric(table(bone$idnum))>1]
bone <- bone[bone$idnum %in% q,]


# Split data by class labels

train <- split(bone, bone$sex) # for classifying sex

plot(1, type="n", xlab = 'Time', ylab = 'Spinal Bone Mineral Density', 
     xlim = c(min(bone$age), max(bone$age)),
     ylim = c(min(bone$spnbmd), max(bone$spnbmd)),
     cex.main = 1.5, cex.lab = 1.5)
for(i in bone$idnum){
  Cr <- ifelse(bone[bone$idnum == i,4] == 'mal', 1, 0)
  x <- as.numeric(bone[bone$idnum == i,3])
  y <- as.numeric(bone[bone$idnum == i,5])
  lines(x, y, type = 'b', 
        lwd = 4, col = rgb(Cr, 0, 0, .3))
}
legend('bottomright', c('Male', 'Female'), lty = 1, lwd = 2, bty = 'n', 
       col = c(rgb(1, 0, 0, 1),col = rgb(0, 0, 0, 1)))










#choose kernel
source("code/RBF.R")
# source("code/laplace.R")
# source("code/matern52.R")
# source("code/matern32.R")



# Fit the model
gpFit.male <- feature(id = train$mal$idnum, x = train$mal$age, y = train$mal$spnbmd)
gpFit.female <- feature(id = train$fem$idnum, x = train$fem$age, y = train$fem$spnbmd)

#save(list = c('gpFit.male', 'gpFit.female'), file = './data/bone/rbf.Rdata')

testTime <- seq(min(bone$age), max(bone$age), length.out = 100)

postDist.male <- postDist(testTime, gpFit.male)
postDist.female <- postDist(testTime, gpFit.female)


# plot
plot(1, type="n", xlab = 'Age (Years)', ylab = 'Spinal Bone Mineral Density', 
     xlim = c(min(bone$age), max(bone$age)),
     ylim = c(min(bone$spnbmd), max(bone$spnbmd)),
     cex.main = 1.5, cex.lab = 1.5)
for(i in bone$idnum){
  Cr <- ifelse(bone[bone$idnum == i,4] == 'mal', 1, 0)
  x <- as.numeric(bone[bone$idnum == i,3])
  y <- as.numeric(bone[bone$idnum == i,5])
  lines(x, y, type = 'b', cex = .5, lwd = .1, col = rgb(Cr, 0, 0 ,0.1))
}
polygon(c(testTime,rev(testTime)), c(postDist.male$ll,rev(postDist.male$ul)), 
        col= rgb(1,0,0,.2),border=NA)
polygon(c(testTime,rev(testTime)), c(postDist.female$ll,rev(postDist.female$ul)), 
        col= rgb(0,0,0,.2),border=NA)
lines(testTime, postDist.male$mu, lwd = 4, col= rgb(1,0,0,1))
lines(testTime, postDist.female$mu, lwd = 4, col= rgb(0,0,0,1))
legend('bottomright', c('Male', 'Female'), lty = 1, lwd = 3, bty = 'n', 
       col = c(2,1))


#subject specific

sub_mu <- postDist_sub(testTime, gpFit.male, sub_id = 1)  


#plot
par(cex.main = 1.5, cex.lab=1.5, cex.axis=1.5, oma = c(1,1,1,1), mar = c(6,6,3,3)) 

#----
plotBone <- function() {
plot(bone$age, bone$spnbmd, type="n", xlab = 'Age (Years)', 
     ylab = 'Spinal Bone Mineral Density')

for(i in bone$idnum){
  Cr <- ifelse(bone[bone$idnum == i,4] == 'mal', 1, 0)
  x <- as.numeric(bone[bone$idnum == i,3])
  y <- as.numeric(bone[bone$idnum == i,5])
  lines(x, y, type = 'b', 
        lwd = 1, col = rgb(Cr, 0, 0, .05))
}
}



unique(train$fem$idnum)
unique(train$mal$idnum)


sub_mu <- postDist_sub(testTime, gpFit.male, sub_id = 1)  
plotBone()
lines(testTime, sub_mu$mu, lwd = 4, col = 2)
lines(as.numeric(bone[bone$idnum == 1,3]), as.numeric(bone[bone$idnum == 1,5]), 
      lwd = 4, col = 2, type = 'b')
sub_mu <- postDist_sub(testTime, gpFit.female, sub_id = 4) 
lines(testTime, sub_mu$mu, lwd = 4, col = 1)
lines(as.numeric(bone[bone$idnum == 4,3]), as.numeric(bone[bone$idnum == 4,5]), 
      lwd = 4, col = 1, type = 'b')
legend('bottomright', c('Male', 'Female'), lty = 1, lwd = 2, bty = 'n', 
       col = c(rgb(1, 0, 0, 1),col = rgb(0, 0, 0, 1)))



sub_mu <- postDist_sub(testTime, gpFit.male, sub_id = 2)
plotBone()
lines(testTime, sub_mu$mu, lwd = 4, col = 2)
lines(as.numeric(bone[bone$idnum == 2,3]), as.numeric(bone[bone$idnum == 2,5]), 
      lwd = 4, col = 2, type = 'b')
sub_mu <- postDist_sub(testTime, gpFit.female, sub_id = 5)
lines(testTime, sub_mu$mu, lwd = 4, col = 1)
lines(as.numeric(bone[bone$idnum == 5,3]), as.numeric(bone[bone$idnum == 5,5]), 
      lwd = 4, col = 1, type = 'b')
legend('bottomright', c('Male', 'Female'), lty = 1, lwd = 2, bty = 'n', 
       col = c(rgb(1, 0, 0, 1),col = rgb(0, 0, 0, 1)))




