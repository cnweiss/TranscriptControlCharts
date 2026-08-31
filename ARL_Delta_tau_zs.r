#Simulate out-of-control zero-state ARL for Delta_tau chart:



#In what follows, p is the vector of transcript or distance frequencies.

#i.i.d. transcript marginal:
pTr.iid <- c(2,2,2,7,7,4)/24

#i.i.d. Cayley marginal:
pC.iid <- c(1,4,7)/12

#i.i.d. Kendall marginal:
pK.iid <- c(1,2,7,2)/12


#Delta_tau or Delta_K:
Delta <- function(p, p0){
	sum((p-p0)^2/p0)
}

#Centered mean Kendall distance:
dK <- function(p){
	sum((0:3)*p) - 11/6
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
np3 <- 6
bincodes <- diag(1,np3)

S3 <- getPerms(1:m)

#Same as vector of strings:
sS3 <- apply(S3, 1, function(x) paste(x, collapse=""))
#sS3
#"123" "132" "213" "231" "312" "321"


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





#Cayley and Kendall distances from Amigo & Dale (2025):
dCayley <- rbind(
c(0, 1, 1, 2, 2, 1),
c(1, 0, 2, 1, 1, 2),
c(1, 2, 0, 1, 1, 2),
c(2, 1, 1, 0, 2, 1),
c(2, 1, 1, 2, 0, 1),
c(1, 2, 2, 1, 1, 0))

dKendall <- rbind(
c(0, 1, 1, 2, 2, 3),
c(1, 0, 2, 3, 1, 2),
c(1, 2, 0, 1, 3, 2),
c(2, 3, 1, 0, 2, 1),
c(2, 1, 3, 2, 0, 1),
c(3, 2, 2, 1, 1, 0))



#Relation between transcripts and distances:
dC2tr <- list(c(1), c(2,3,6), c(4,5)) #0-2
dK2tr <- list(c(1), c(2,3), c(4,5), c(6)) #0-3

#Resulting transformation matrices:
TC <- array(0, c(length(dC2tr), length(sS3)))
for(i in 1:length(dC2tr)) TC[i, dC2tr[[i]]] <- 1

TK <- array(0, c(length(dK2tr), length(sS3)))
for(i in 1:length(dK2tr)) TK[i, dK2tr[[i]]] <- 1





#Run length simulations
reps <- 1e5

# lam <- 0.25
# lam <- 0.1
lam <- 0.05


#Delta_tau chart:
cl <- 3.225 #370.09066	367.548915687176
cl <- 0.9685 #369.75601	360.548284657597
cl <- 0.4328 #370.24635	354.504668879382


#Out-of-control model: AR(1)
taba1 <- c(-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8)

# #Out-of-control model: AR(1)-type model X_t = a1*abs(X_{t-1})+eps[t]
# prerun <- 100
# taba1 <- c(1:4)/5

# #Out-of-control model: AR(1)-type model X_t = a1*X_{t-1}^2+eps[t]
# prerun <- 100
# taba1 <- c(3:6)/20

# #Out-of-control model: TEAR(1)
# prerun <- 100
# taba1 <- c(1:6)/10

# #Out-of-control model: QMA(1)
# taba1 <- c(1:4)/5


set.seed(1111)

results <- c()
for(a1 in taba1){
rls <- rep(NA, reps)
for(r in 1:reps){
	p <- pTr.iid
	# pK <- pK.iid

	stat <- Delta(p, pTr.iid)
	
	rl <- 1 #start with 1 to make it comparable to OP-EWMA charts
	
	#Initial OP as index:
	
	#AR(1):
	X <- rnorm(1,0,sqrt(1/(1-a1^2)))
	win <- rnorm(m)
	for(t in 1:m){
		X <- win[t] <- win[t] + a1*X
	}
	
	# #AAR(1), pre-run:
	# eps <- rnorm(prerun)
	# X <- 0
	# for(t in 1:prerun) X <- a1*abs(X) + eps[t]
	
	# win <- rnorm(m)
	# for(t in 1:m){
		# X <- win[t] <- a1*abs(X) + win[t]
	# }
	
	# #QAR(1), pre-run:
	# eps <- rnorm(prerun)
	# X <- 0
	# for(t in 1:prerun) X <- a1*X^2 + eps[t]
	
	# win <- rnorm(m)
	# for(t in 1:m){
		# X <- win[t] <- a1*X^2 + win[t]
	# }
	
	# #TEAR(1), pre-run:
	# eps <- rexp(prerun, 1)
	# tabR <- rbinom(prerun, 1, a1)
	# X <- 1
	# for(t in 1:prerun) X <- (1-a1)*eps[t] + tabR[t]*X
	
	# win <- rexp(m, 1)
	# tabR <- rbinom(m, 1, a1)
	# for(t in 1:m){
		# X <- win[t] <- (1-a1)*win[t] + tabR[t]*X
	# }
	
	# #QMA(1):
	# eps <- rnorm(2)
	# win <- rep(NA, m)
	# for(t in 1:m){
		# win[t] <- eps[2]+a1*eps[1]^2
		# eps <- c(eps[-1], rnorm(1))
	# }

	iop1 <- match(paste(order(win), collapse=""), sS3)
	
	while(stat<cl){
		rl <- rl+1
		
		#Next OP as index:
		
		#AR(1):
		X <- rnorm(1) + a1*X
		
		# #AAR(1):
		# X <- a1*abs(X) + rnorm(1)
		
		# #QAR(1):
		# X <- a1*X^2 + rnorm(1)
		
		# #TEAR(1):
		# X <- (1-a1)*rexp(1, 1) + rbinom(1, 1, a1)*X
		
		# #QMA(1):
		# X <- eps[2]+a1*eps[1]^2
		# eps <- c(eps[-1], rnorm(1))

		win <- c(win[-1], X)
		iop <- match(paste(order(win), collapse=""), sS3)
		
		#Resulting transcript (numerical index or distance value):
		trOP <- transcript(iop1, iop)
		
		#EWMA-smoothed transcript frequencies:
		p <- lam*bincodes[trOP,] + (1-lam)*p
		
		stat <- Delta(p, pTr.iid)
		
		iop1 <- iop
	}
	rls[r] <- rl
} #for reps

results <- rbind( results, c(a1, mean(rls), sd(rls)) )
} #for a1

print(results)


