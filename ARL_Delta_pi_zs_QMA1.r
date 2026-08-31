#Simulate out-of-control zero-state ARL for Delta_pi chart:



#Delta^2:
Delta2 <- function(p){
	puni <- 1/length(p)
	sum((p-puni)^2)
}

#beta for m=3:
beta3 <- function(p){
	p[6]-p[1]
}

#tau for m=3:
tau3 <- function(p){
	p[6]+p[1]-1/3
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


#If m <- 3
S3 <- rbind( c(3,2,1), c(3,1,2), c(2,3,1), c(1,3,2), c(2,1,3), c(1,2,3))
np3 <- dim(S3)[1]

#Same as vector of strings:
sS3 <- apply(S3, 1, function(x) paste(x, collapse=""))
#sS3
#[1] "321" "312" "231" "132" "213" "123"



#Run length simulations
m <- 3
np3 <- 6
bincodes <- diag(1,np3)
reps <- 1e5


# #In-control designs from Weiß & Testik (2023), Table 1:

# lam <- 0.25

# clDelta2 <- 0.3338 #369.3
# clbeta3 <- 0.6437 #369.7
# cltau3 <- 0.4253 #368.6


# lam <- 0.1

# clDelta2 <- 0.1115 #369.5
# clbeta3 <- 0.3638 #369.9
# cltau3 <- 0.2529 #369.8


lam <- 0.05

clDelta2 <- 0.05125 #370.3
clbeta3 <- 0.233 #370.0
cltau3 <- 0.16775 #369.4


#Out-of-control model: QMA(1)
taba1 <- c(1:4)/5


set.seed(4567)

results <- c()
for(a1 in taba1){
rls <- rep(NA, reps)
for(r in 1:reps){
	p <- rep(1/np3, np3)
	stat <- Delta2(p)
	
	rl <- 0
	eps <- rnorm(2)
	win <- rep(NA, m)
	for(t in 1:m){
		win[t] <- eps[2]+a1*eps[1]^2
		eps <- c(eps[-1], rnorm(1))
	}
	while(stat<clDelta2){
		rl <- rl+1
		
		#binarized ordinal pattern:
		p <- lam*bincodes[match(paste(rank(win, ties.method="first"), collapse=""), sS3),] + (1-lam)*p
		stat <- Delta2(p)
		
		win <- c(win[-1], eps[2]+a1*eps[1]^2)
		eps <- c(eps[-1], rnorm(1))
	}
	rls[r] <- rl
} #for reps

results <- rbind( results, c(a1, mean(rls), sd(rls)) )
} #for ooc

print(results)


