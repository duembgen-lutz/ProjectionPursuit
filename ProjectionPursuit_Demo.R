source("ProjectionPursuit.R")

# Test of ICS plus LocalPP in dimension 2 ----

## Generate white noise data matrix:
n <- 500
p <- 10
X0 <- matrix(rnorm(n*p),n,p)

## Generate interesting structure (choose one):

# Interesting structure 1:
X1 <- X0
R2 <- sqrt(2/rowSums(X1[,1:2]^2))
X1[,1:2] <- R2*X1[,1:2]
cov(X1)
sqrt(mean((cov(X1) - diag(p))^2))
pairs(X1)

# Interesting structure 2:
X1 <- X0
b <- 0.5
delta <- 5
X1[,1] <- X1[,1] + delta*(runif(n) <= b)
X1[,1] <- X1[,1]/sqrt(1 + delta^2*b*(1 - b))
sqrt(mean((cov(X1) - diag(p))^2))
pairs(X1)

## Hide this structure:
C <- HaarMatrix(p)
Xraw <- X1 %*% t(C)
pairs(Xraw)

## Search for this structure:

# Basic procedure 1: Pre-whiten the data
Xpre <- PreWhiten(Xraw)
pairs(Xpre)

# Basic procedure 2: ICS
Xics <- ICS2(Xpre,nu=0,exponent=1)
pairs(Xics)

# Basic procedure 3: LocalPP

# Bandwidth for empirical entropy:
h <- 0.5

# BP3a: Look for promising 2-dim. standard projections
GlobalPP2_pre(Xics,h=h)
# Rearranged data (only if standard projection onto components
#  i, j rather than 1, 2 is wanted!):
i <- ?
j <- ?
Xics <- cbind(Xics[,c(i,j)],Xics[,(1:p)[-c(i,j)]])

Xfinal <- LocalPP(Xics,d=2,h,showsteps=FALSE)
# Show up to 40 intermediate steps by clicking on the plot:
Xfinal <- LocalPP(Xics,d=2,h,showsteps=TRUE,maxiter=40)

# How do we get from Xraw to Xfinal:
res <- RecoverTransformation(Xraw,Xfinal)
res$a
res$B
# Xfinal[,j] = res$a + t(res$B) %*% Xraw[,j]

# Just checking:
Error <- Xfinal - rep(1,n) %*% t(res$a) - Xraw %*% res$B
sqrt(sum(Error^2))


# Test of ICS plus LocalPP in dimension 3 ----

## Generate white noise data matrix:
n <- 500
p <- 10
X0 <- matrix(rnorm(n*p),n,p)

# Generate interesting structure (choose one):

# Structure 1 (four clusters):
Tetraeder <- matrix(0,4,3)
Tetraeder[1,1] <- 1
Tetraeder[2:4,1] <- -1/3
Tetraeder[2,2] <- sqrt(8/9)
Tetraeder[3:4,2] <- - sqrt(2/9)
Tetraeder[3:4,3] <- sqrt(2/3)*c(1,-1)
Tetraeder
X1 <- X0
b <- 1.65
for (i in 1:n){
	X1[i,1:3] <- sqrt(1 - b^2/3)*X1[i,1:3] +
		b*Tetraeder[ceiling(4*runif(1)),]
}
mean((cov(X1) - diag(p))^2)
C <- HaarMatrix(p)
Xraw <- X1 %*% t(C)
pairs(Xraw)

# Structure 2 (data on a 2-dim. manifold in a 3-dim. subspace):
Y <- X0[,1]*X0[,2]
C <- HaarMatrix(p)
Xraw <- cbind(X0 %*% t(C),Y)
p <- p+1

pairs(Xraw)

## Search for this structure:

# Basic procedure 1: Pre-whiten the data
Xpre <- PreWhiten(Xraw)
pairs(Xpre)

# Basic procedure 2: ICS
Xics <- ICS2(Xpre,nu=0,exponent=1)
pairs(Xics)

# Basic procedure 3: LocalPP

# Bandwidth for empirical entropy:
h <- 0.5

# BP3a: Look for promising 2-dim. standard projections
GlobalPP3_pre(Xics,h=h)
# Rearranged data (only if standard projection onto components
#  i, j, k rather than 1, 2, 3 is wanted!):
i <- 8
j <- 9
k <- 11
Xics <- cbind(Xics[,c(i,j,k)],Xics[,(1:p)[-c(i,j,k)]])

Xfinal <- LocalPP(Xics,d=3,h,showsteps=FALSE)
# Show up to 40 intermediate steps by clicking on the plot:
Xfinal <- LocalPP(Xics,d=3,h,showsteps=TRUE,maxiter=40)

# Explore the 3-dim. projection:
source("Interactive3D.R")
tmp <- Interactive3D(Xfinal,center.X=FALSE)


# How do we get from Xraw to Xfinal:
res <- RecoverTransformation(Xraw,Xfinal)
res$a
res$B
# Xfinal[,j] = res$a + t(res$B) %*% Xraw[,j]

# Just checking:
Error <- Xfinal - rep(1,n) %*% t(res$a) - Xraw %*% res$B
sqrt(sum(Error^2))
