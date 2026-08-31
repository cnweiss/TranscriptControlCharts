#Control charts for chemical process data of O'Donovan, T. M. (1983) Short Term Forecasting: An Introduction to the Box–Jenkins Approach. London: Wiley.
#Time series 4.1 (printed in Appendix A.3, discussed in Section 4.5).

data <- c(40,54,48,52,41,52,38,56,48,45,66,17,62,50,38,59,51,55,48,51,50,52,44,65,40,65,41,64,53,48,53,43,66,48,52,42,44,56,44,58,41,54,51,56,38,56,49,52,33,52,59,34,57,39,60,40,52,44,65,43,48,44,49,44,49,69,40,54,58,49)
T <- length(data)
T #70

plot(data, type="b", xlab="Batch", ylab=expression("Batch yields  x"[t]), pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,70))

hist(data, xlab="Batch yields")

acf(data, lag.max=10, ylim=c(-1,1), lwd=3, ci.col=grey(0.5), main="Batch yields")


acf(data)[[1]][-1] #-0.58842235  0.27151541 -0.07750074  0.01834453  0.09646180 -0.14868820 ...
c(pacf(data)[[1]]) #-0.58842235 -0.11430120  0.05156247  0.03136143  0.15755861 -0.03330803 ...




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





m <- 3 #order of OPs
n <- T-(m-1) #length of OP-series


#For technical reasons, define array of lagged time series:
data0 <- data[1:n]
for(k in 1:(m-1)){
	data0 <- cbind(data0, data[(k+1):(k+n)])
}

#Check for ties:
ties <- apply(data0, 1, function(x) (max(table(x))>1))
data0[ties,]
# [1,]    52 41 52
# [2,]    65 40 65
# [3,]    53 48 53
# [4,]    44 56 44
# [5,]    56 38 56
# [6,]    52 33 52
# [7,]    44 49 44
# [8,]    49 44 49



#First, proceed without jittering.
#Jittering experiment in ChemProcess_Jittering.r.


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
p <- rep(1/np3, np3)

tabp <- array(NA, c(n, np3))
for(t in 1:n){
	p <- lam*databin[t,] + (1-lam)*p
	tabp[t,] <- p
} #for t
#matplot(1:n, tabp, type="l")



#Delta-chart:
plot(c(3:T), apply(tabp, 1, Delta2), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,0.4), xlab="t", ylab=expression(Delta[pi]-chart))
abline(h=0.3338) #no alarm
abline(h=0, lty=2)

#tau-chart:
plot(c(3:T), apply(tabp, 1, tau3), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(-0.45,0.45), xlab="t", ylab=expression(tau-chart))
abline(h=0.4253*c(-1,1)) #no alarm
abline(h=0, lty=2)





#EWMA smoothing of OP-probabilities with lambda=0.10:

lam <- 0.10
p <- rep(1/np3, np3)

tabp <- array(NA, c(n, np3))
for(t in 1:n){
	p <- lam*databin[t,] + (1-lam)*p
	tabp[t,] <- p
} #for t
matplot(1:n, tabp, type="l")



#Delta-chart:
plot(c(3:T), apply(tabp, 1, Delta2), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,0.15), xlab="t", ylab=expression(Delta[pi]-chart))
abline(h=0.1115) #no alarm
abline(h=0, lty=2)

#tau-chart:
plot(c(3:T), apply(tabp, 1, tau3), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(-0.35,0.35), xlab="t", ylab=expression(tau-chart))
abline(h=0.2529*c(-1,1)) #alarm at 26
abline(h=0, lty=2)
abline(v=26, lty=3)






#EWMA smoothing of OP-probabilities with lambda=0.05:

lam <- 0.05
p <- rep(1/np3, np3)

tabp <- array(NA, c(n, np3))
for(t in 1:n){
	p <- lam*databin[t,] + (1-lam)*p
	tabp[t,] <- p
} #for t
matplot(1:n, tabp, type="l")



#Delta-chart:
plot(c(3:T), apply(tabp, 1, Delta2), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,0.075), xlab="t", ylab=expression(Delta[pi]-chart))
abline(h=0.05125) #alarm at 62
abline(h=0, lty=2)
abline(v=62, lty=3)

#tau-chart:
plot(c(3:T), apply(tabp, 1, tau3), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(-0.3,0.3), xlab="t", ylab=expression(tau-chart))
abline(h=0.16775*c(-1,1)) #alarm at 24
abline(h=0, lty=2)
abline(v=24, lty=3)





###################
#Transcript charts:
###################


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

p.tr.iid <- c(2,2,2,7,7,4)/24
p.K.iid <- c(1,2,7,2)/12
mu.K.iid <- 11/6

colMeans(trbin)
#0.00000000 0.01492537 0.02985075 0.25373134 0.28358209 0.41791045
p.tr.iid
#0.08333333 0.08333333 0.08333333 0.29166667 0.29166667 0.16666667

colMeans(trbin %*% t(TK))
#0.00000000 0.04477612 0.53731343 0.41791045
p.K.iid
#0.08333333 0.16666667 0.58333333 0.16666667

mean(trbin %*% t(TK) %*% (0:3)) #2.373134
mu.K.iid #1.833333




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
plot(c(4:T), apply(tabptr, 1, Delta, p0=p.tr.iid), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,3.5), xlab="t", ylab=expression(Delta[tau]-chart))
abline(h=3.225) #no alarm
abline(h=0, lty=2)

#Delta.K-chart:
plot(c(4:T), apply(tabpK, 1, Delta, p0=p.K.iid), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,3.5), xlab="t", ylab=expression(Delta[K]-chart))
abline(h=3.19) #no alarm
abline(h=0, lty=2)

#mu.K-chart:
plot(c(4:T), apply(tabpK, 1, dK), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(-1.2,1.2), xlab="t", ylab=expression(mu[K]-chart))
abline(h=1.0188*c(-1,1)) #no alarm
abline(h=0, lty=2)




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
plot(c(4:T), apply(tabptr, 1, Delta, p0=p.tr.iid), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,1.5), xlab="t", ylab=expression(Delta[tau]-chart))
abline(h=0.9685)
abline(h=0, lty=2)
abline(v=26, lty=3)

#Delta.K-chart:
plot(c(4:T), apply(tabpK, 1, Delta, p0=p.K.iid), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,1.5), xlab="t", ylab=expression(Delta[K]-chart))
abline(h=0.8078)
abline(h=0, lty=2)
abline(v=25, lty=3)

#mu.K-chart:
plot(c(4:T), apply(tabpK, 1, dK), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(-0.8,0.8), xlab="t", ylab=expression(mu[K]-chart))
abline(h=0.5827*c(-1,1))
abline(h=0, lty=2)
abline(v=24, lty=3)




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
plot(c(4:T), apply(tabptr, 1, Delta, p0=p.tr.iid), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,1), xlab="t", ylab=expression(Delta[tau]-chart))
abline(h=0.4328)
abline(h=0, lty=2)
abline(v=25, lty=3)

#Delta.K-chart:
plot(c(4:T), apply(tabpK, 1, Delta, p0=p.K.iid), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(0,0.9), xlab="t", ylab=expression(Delta[K]-chart))
abline(h=0.3229)
abline(h=0, lty=2)
abline(v=24, lty=3)

#mu.K-chart:
plot(c(4:T), apply(tabpK, 1, dK), type="b", pch=19, cex=0.5, xlim=c(0,T), ylim=c(-0.5,0.7), xlab="t", ylab=expression(mu[K]-chart))
abline(h=0.3785*c(-1,1))
abline(h=0, lty=2)
abline(v=23, lty=3)



