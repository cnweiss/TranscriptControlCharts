#Control charts for chemical process data of O'Donovan, T. M. (1983),
#investigate effect of possible jittering.

data <- c(40,54,48,52,41,52,38,56,48,45,66,17,62,50,38,59,51,55,48,51,50,52,44,65,40,65,41,64,53,48,53,43,66,48,52,42,44,56,44,58,41,54,51,56,38,56,49,52,33,52,59,34,57,39,60,40,52,44,65,43,48,44,49,44,49,69,40,54,58,49)
T <- length(data)
T #70





##########################
#Ordinal pattern analysis:
##########################

#Statistics discussed in Section 2 of the manuscript:

#Delta^2:
Delta2 <- function(p){
	puni <- 1/length(p)
	sum((p-puni)^2)
}

#beta for m=3:
beta3 <- function(p){
	p[1]-p[6]
}

#tau for m=3:
tau3 <- function(p){
	p[1]+p[6]-1/3
}





getPerms <- function(x) {
    if (length(x) == 1) {
        return(x)
    }
    else {
        res <- matrix(nrow = 0, ncol = length(x))
        for (i in seq_along(x)) {
            res <- rbind(res, cbind(x[i], Recall(x[-i])))
        }
        return(res)
    }
}





#OPs of order 3:
m <- 3

S3 <- getPerms(1:m)

#Same as vector of strings:
sS3 <- apply(S3, 1, function(x) paste(x, collapse=""))
#sS3
#"123" "132" "213" "231" "312" "321"

np3 <- length(sS3)
bincodes <- diag(1,np3)


#For transcripts, OPs use permutation representation corresponding to "order()".


#Cayley table of permutations:
sCayley <- array("", c(factorial(m),factorial(m)))
for(i in 1:factorial(m)){
	sCayley[i,] <- apply(S3, 1, function(x) paste(x[S3[i,]], collapse=""))
}


#Change table into "index version":
iCayley <- array(0, c(factorial(m),factorial(m)))
for(i in 1:factorial(m)){
for(j in 1:factorial(m)){
	iCayley[i,j] <- which(sS3==sCayley[i,j])
}}

#Neutral element is "123" (always first element in lexicographic order), so inverses of sS3 are:
inv.iS3 <- rep(0, factorial(m))
for(i in 1:factorial(m)){
	inv.iS3[i] <- which(sS3[1]==sCayley[i,])
}
inv.sS3 <- sS3[inv.iS3]
inv.sS3 #"123" "132" "213" "312" "231" "321"
inv.iS3 #1 2 3 5 4 6


#Transcript tau(a,b) = b * inv(a):
transcript <- function(a,b){
	#a,b are index (position) within sS3
	iCayley[b, inv.iS3[a]]
}

#Cayley tables of transcript:
iTrans <- array(0, c(factorial(m),factorial(m)))
sTrans <- array("", c(factorial(m),factorial(m)))
for(i in 1:factorial(m)){
for(j in 1:factorial(m)){
	iTrans[i,j] <- transcript(i,j)
	sTrans[i,j] <- sS3[iTrans[i,j]]
}}


#Pairs (OP1,OP2) leading to transcript k:
iPairs <- vector(mode="list", length=factorial(m))
for(i in 1:factorial(m)){
for(j in 1:factorial(m)){
	iPairs[[iTrans[i,j]]] <- rbind(iPairs[[iTrans[i,j]]], c(i,j))
}}


#Kendall distances from Amigo & Dale (2025):
dKendall <- rbind(
c(0, 1, 1, 2, 2, 3),
c(1, 0, 2, 3, 1, 2),
c(1, 2, 0, 1, 3, 2),
c(2, 3, 1, 0, 2, 1),
c(2, 1, 3, 2, 0, 1),
c(3, 2, 2, 1, 1, 0))

#Relation between transcripts and Kendall distance:
dK2tr <- list(c(1), c(2,3), c(4,5), c(6)) #0-3

#Resulting transformation matrix:
TK <- array(0, c(length(dK2tr), length(sS3)))
for(i in 1:length(dK2tr)) TK[i, dK2tr[[i]]] <- 1


#Delta_tau or Delta_K:
Delta <- function(p, p0){
	sum((p-p0)^2/p0)
}

#Re-scaled Kendall distance:
dK <- function(p){
	sum((0:3)*p) - 11/6
}





n <- T-(m-1) #length of OP-series

p.op.iid <- rep(1/np3, np3)

p.tr.iid <- c(2,2,2,7,7,4)/24
p.K.iid <- c(1,2,7,2)/12
mu.K.iid <- 11/6


#If rate of ties considered as not being negligible,
#jittering has to be used.
#Analyze how jittering affects the chart decisions.
#Data are integers, so U(0,1)-jittering.

reps <- 1e4

#Collect time of first alarm, where no alarm labeled as "NA":
alarm.times.op <- 3:T
first.alarm.op <- array(NA, c(reps,2*3)) #two charts, three lambda

alarm.times.tr <- 4:T
first.alarm.tr <- array(NA, c(reps,3*3)) #three charts, three lambda




set.seed(123)

for(r in 1:reps){
data.jitt <- data + runif(T)


#For technical reasons, define array of lagged time series:
data0 <- data.jitt[1:n]
for(k in 1:(m-1)){
	data0 <- cbind(data0, data.jitt[(k+1):(k+n)])
}


#Generate ordinal patterns as rank vectors:
opts <- t(matrix(c(apply(data0, 1, order)), nrow=m))

rm(data0) #data0 not needed anymore

#transform into vector of strings:
sopts <- apply(opts, 1, function(x) paste(x, collapse=""))

rm(opts) #opts not needed anymore

#So sopts is categorical time series with range sS3.

#Numeric coding by 1-m!:
datanum <- match(sopts, sS3)

#Binarization (vectors Y_t of Section 2):
databin <- bincodes[datanum,] #assign row vectors







##############################################
#Construct EWMA control chart based on databin
##############################################

#EWMA smoothing of OP-probabilities with lambda=0.25:

lam <- 0.25
p <- p.op.iid

tabp <- array(NA, c(n, np3))
for(t in 1:n){
	p <- lam*databin[t,] + (1-lam)*p
	tabp[t,] <- p
} #for t


#Delta-chart:
stats <- apply(tabp, 1, Delta2)
alarms <- (stats>0.3338)
if(sum(alarms>0)) first.alarm.op[r,1] <- min(alarm.times.op[alarms])

#tau-chart:
stats <- apply(tabp, 1, tau3)
alarms <- (abs(stats)>0.4253)
if(sum(alarms>0)) first.alarm.op[r,4] <- min(alarm.times.op[alarms])





#EWMA smoothing of OP-probabilities with lambda=0.10:

lam <- 0.10
p <- p.op.iid

tabp <- array(NA, c(n, np3))
for(t in 1:n){
	p <- lam*databin[t,] + (1-lam)*p
	tabp[t,] <- p
} #for t


#Delta-chart:
stats <- apply(tabp, 1, Delta2)
alarms <- (stats>0.1115)
if(sum(alarms>0)) first.alarm.op[r,2] <- min(alarm.times.op[alarms])

#tau-chart:
stats <- apply(tabp, 1, tau3)
alarms <- (abs(stats)>0.2529)
if(sum(alarms>0)) first.alarm.op[r,5] <- min(alarm.times.op[alarms])





#EWMA smoothing of OP-probabilities with lambda=0.05:

lam <- 0.05
p <- p.op.iid

tabp <- array(NA, c(n, np3))
for(t in 1:n){
	p <- lam*databin[t,] + (1-lam)*p
	tabp[t,] <- p
} #for t


#Delta-chart:
stats <- apply(tabp, 1, Delta2)
alarms <- (stats>0.05125)
if(sum(alarms>0)) first.alarm.op[r,3] <- min(alarm.times.op[alarms])

#tau-chart:
stats <- apply(tabp, 1, tau3)
alarms <- (abs(stats)>0.16775)
if(sum(alarms>0)) first.alarm.op[r,6] <- min(alarm.times.op[alarms])





###################
#Transcript charts:
###################



#OP series in numerical coding in datanum, being of length n.
#Resulting transcript series in numerical coding:

trOP <- rep(NA, n-1)
for(t in 2:n){
	trOP[t-1] <- transcript(datanum[t-1],datanum[t])
}

#Binarization (vectors Z_t of Section 4):
trbin <- bincodes[trOP,] #assign row vectors





############################################
#Construct EWMA control chart based on trbin
############################################


#EWMA smoothing of tr-probabilities with lambda=0.25:

lam <- 0.25
ptr <- p.tr.iid

tabptr <- array(NA, c(n-1, np3))
for(t in 1:(n-1)){
	ptr <- lam*trbin[t,] + (1-lam)*ptr
	tabptr[t,] <- ptr
} #for t
#matplot(1:n, tabptr, type="l")

tabpK <- tabptr %*% t(TK)



#Delta.tr-chart:
stats <- apply(tabptr, 1, Delta, p0=p.tr.iid)
alarms <- (stats>3.225)
if(sum(alarms>0)) first.alarm.tr[r,1] <- min(alarm.times.tr[alarms])

#Delta.K-chart:
stats <- apply(tabpK, 1, Delta, p0=p.K.iid)
alarms <- (stats>3.19)
if(sum(alarms>0)) first.alarm.tr[r,4] <- min(alarm.times.tr[alarms])

#mu.K-chart:
stats <- apply(tabpK, 1, dK)
alarms <- (abs(stats)>1.0188)
if(sum(alarms>0)) first.alarm.tr[r,7] <- min(alarm.times.tr[alarms])




#EWMA smoothing of tr-probabilities with lambda=0.10:

lam <- 0.10
ptr <- p.tr.iid

tabptr <- array(NA, c(n-1, np3))
for(t in 1:(n-1)){
	ptr <- lam*trbin[t,] + (1-lam)*ptr
	tabptr[t,] <- ptr
} #for t
#matplot(1:n, tabptr, type="l")

tabpK <- tabptr %*% t(TK)



#Delta.tr-chart:
stats <- apply(tabptr, 1, Delta, p0=p.tr.iid)
alarms <- (stats>0.9685)
if(sum(alarms>0)) first.alarm.tr[r,2] <- min(alarm.times.tr[alarms])

#Delta.K-chart:
stats <- apply(tabpK, 1, Delta, p0=p.K.iid)
alarms <- (stats>0.8078)
if(sum(alarms>0)) first.alarm.tr[r,5] <- min(alarm.times.tr[alarms])

#mu.K-chart:
stats <- apply(tabpK, 1, dK)
alarms <- (abs(stats)>0.5827)
if(sum(alarms>0)) first.alarm.tr[r,8] <- min(alarm.times.tr[alarms])




#EWMA smoothing of tr-probabilities with lambda=0.05:

lam <- 0.05
ptr <- p.tr.iid

tabptr <- array(NA, c(n-1, np3))
for(t in 1:(n-1)){
	ptr <- lam*trbin[t,] + (1-lam)*ptr
	tabptr[t,] <- ptr
} #for t
#matplot(1:n, tabptr, type="l")

tabpK <- tabptr %*% t(TK)



#Delta.tr-chart:
stats <- apply(tabptr, 1, Delta, p0=p.tr.iid)
alarms <- (stats>0.4328)
if(sum(alarms>0)) first.alarm.tr[r,3] <- min(alarm.times.tr[alarms])

#Delta.K-chart:
stats <- apply(tabpK, 1, Delta, p0=p.K.iid)
alarms <- (stats>0.3229)
if(sum(alarms>0)) first.alarm.tr[r,6] <- min(alarm.times.tr[alarms])

#mu.K-chart:
stats <- apply(tabpK, 1, dK)
alarms <- (abs(stats)>0.3785)
if(sum(alarms>0)) first.alarm.tr[r,9] <- min(alarm.times.tr[alarms])

} #for reps

apply(first.alarm.op, 2, table, useNA="ifany")
# [[1]]
 # <NA> 
# 10000 

# [[2]]
 # <NA> 
# 10000 

# [[3]]
  # 48   49   50   58   61   62   63 
 # 651  626 1848  920  636 5162  157 

# [[4]]
 # <NA> 
# 10000 

# [[5]]
   # 26 
# 10000 

# [[6]]
   # 24 
# 10000

apply(first.alarm.tr, 2, table, useNA="ifany")
# [[1]]
 # <NA> 
# 10000 

# [[2]]
  # 26   27   28 
# 4940 2514 2546 

# [[3]]
  # 25   26   28 
# 4987 2467 2546 

# [[4]]
 # <NA> 
# 10000 

# [[5]]
   # 25 
# 10000 

# [[6]]
  # 24   25 
# 4987 5013 

# [[7]]
 # <NA> 
# 10000 

# [[8]]
  # 24   25 
# 4987 5013 

# [[9]]
  # 23   24 
# 4987 5013

