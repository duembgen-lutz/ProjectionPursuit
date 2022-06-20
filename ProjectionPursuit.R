# Various functions in relation to Projection Pursuit
# with empirical entropy as target function to be
# minimized.
# 
# Lutz Duembgen, Katrin Gysel, March 2022

source("Interactive3D.R")

# Prewhitening and ICS2 ----

PreWhiten <- function(X0)
	# For a numerical data matrix or data frame X0,
	# this procedure computes a pre-whitended data
	# matrix X. That means, the rows of X are affine
	# transformations of the rows of X0 such that
	# the sample means equals 0 and the sample
	# covariance matrix is the identity matrix.
{
	if (is.data.frame(X0)){
		X <- as.matrix(X0)
	}else{
		X <- X0
	}
	n <- dim(X)[1]
	q <- dim(X)[2]
	# Centering:
	X <- X - rep(1,n) %*% t(colMeans(X))
	# Multivariate rescaling:
	S <- t(X) %*% X / (n-1)
	tmp <- eigen(S)
	X <- X %*% tmp$vectors
	X <- X %*% diag(1/sqrt(tmp$values))
	return(X)
}	

OneStepMsymm <- function(X,nu=1,exponent=1)
	# Computes a one-step symmetrized M-estimator of scatter,
	# assuming prewhitened data X.
{
	n <- dim(X)[1]
	p <- dim(X)[2]
	tmp <- X[n,] - X[n-1,]
	S <- tmp %*% (t(tmp) / (nu + sum(tmp^2))^exponent)
	for (i in 1:(n-2)){
		tmp <- t(X[(i+1):n,]) - X[i,]
		S <- S +
			tmp %*% (t(tmp) / (nu + colSums(tmp^2))^exponent)
	}
	S <- (2*(p+nu)^exponent/n/(n-1))*S
	return(S)
}

ICS2 <- function(X,nu=1,exponent=1)
	# Computes an ICS transformation of prewhitened
	# data X by means of a one-step
	# symmetrized M-estimator of scatter.
{
	S <- OneStepMsymm(X,nu,exponent)
	V <- eigen(S)$vectors
	Xics <- X %*% V
	return(Xics)
}

# Projection pursuit based on empirical entropy ----

StandardEntropy <- function(h=0.5,d=2)
	# If we estimate the differential Shannon
	# entropy with a Gaussian kernel density
	# estimator with bandwidth h, the estimator
	# is close to StandardEntropy(h,d) in case
	# of data from a standard Gaussian distribution.
{
	SE <- (d/2)*(1/(1 + h^2) + log(1 + h^2) + log(2*pi))
	return(SE)
}

EmpEntropy <- function(Y,h=0.5)
	# For a data matrix Y consisting of n vectors
	# Y[i,], 1 <= i <= n, this procedure computes
	# the estimated differential Shannon entropy
	# of the underlying distribution.
{
	if (is.vector(Y)){
		n <- length(Y)
		d <- 1
	}else{
		n <- dim(Y)[1]
		d <- dim(Y)[2]
	}
	lsf <- d*(log(2*pi)/2 + log(h))
	K <- Y %*% t(Y) / h^2
	K <- exp(K -
			diag(K) %*% t(rep(1,n))/2 -
			rep(1,n) %*% t(diag(K))/2)
	gh <- colMeans(K)
	entropy <- - mean(log(gh)) + lsf
	return(entropy)
}

GradEmpEntropy <- function(X,d,h=0.5)
	# For a data matrix X consisting of n vectors
	# X[i,], 1 <= i <= n, this procedure computes
	# the estimated Shannon entropy of the
	# underlying distribution of the subvectors
	# X[i,1:d] and the negative Gradient of this
	# function if the X[i,] are rotated slightly.
{
	n <- dim(X)[1]
	q <- dim(X)[2]
	lsf <- d*(log(2*pi)/2 + log(h))
	K <- X[,1:d] %*% t(X[,1:d]) / h^2
	K <- exp(K -
			diag(K) %*% t(rep(1,n))/2 -
			rep(1,n) %*% t(diag(K))/2)
	ngh <- colSums(K)
	H <- log(n) - mean(log(ngh)) + lsf
	K <- K/ngh
	C <- matrix(0,q-d,d)
	for (i in 1:n){
		Xi <- t(t(X) - X[i,])
		C <- C + t(Xi[,(d+1):q]) %*% (K[i,]*Xi[,1:d])
	}
	C <- C/(n*h^2)
	return(list(H=H,C=C))
}

Exp <- function(V,sigma,W)
	# computes exp(D) for the skew-symmetric matrix D
	# with blocks
	#    0, - t(C),
	#    C, 0,
	# where C = W %*% diag(sigma) %*% t(V) with
	# matrices V,W having orthonormal columns and
	# a vector sigma of singular values.
{
	m <- dim(V)[2]
	d <- dim(V)[1]
	qmd <- dim(W)[1]
	q <- d + qmd
	ii <- 1:d
	jj <- (d+1):q
	M <- matrix(0,q,q)
	M[jj,ii] <- W %*% diag(sin(sigma),nrow=m) %*% t(V)
	M[ii,jj] <- - t(M[jj,ii])
	cv1 <- cos(sigma)
	if (d != qmd){
		cv2 <- 2*sin(sigma/2)^2
	}
	if (d == m){
		M[ii,ii] <- V %*% diag(cv1,nrow=m) %*% t(V)
	}else{
		M[ii,ii] <- diag(d) -
			V %*% diag(cv2,nrow=m) %*% t(V)
	}
	if (qmd == m){
		M[jj,jj] <- W %*% diag(cv1,nrow=m) %*% t(W)
	}else{
		M[jj,jj] <- diag(qmd) -
			W %*% diag(cv2,nrow=m) %*% t(W)
	}
	return(M)
}


LocalPP <- function(X,d=2,h=0.5,start=1:d,
	delta0=10^(-11),
	display=TRUE,
	showsteps=TRUE,
	maxiter=500)
	# Performs a local minimization of the empirical
	# entropy of d-dimensional projections of the data X.
	# Starting point is the projection on the components
	# specified by start.
	# Input:
	# -  X : data matrix (n x q)
	# -  d : dimension of projection (default 2)
	# -  h : bandwidth of kernel densitiy estimator
	# -  start : tuple of d indices specifying the
	#        starting projection
	# -  delta0 : threshold for directional derivative
	# -  display : if TRUE, the scatter plot matrix of
	#        the first min(d,3) components is shown
	#        initially and after the local search.
	# -  showsteps : if TRUE, the scatter plot matrix
	#        of the first min(d,3) components is shown
	#        after each iteration.
	#        (To get to the next iteration, one has to
	#         click on the scatter plot matrix.)
	# Output: A list of two items:
	# -  Xtrans : the transformed data matrix
	# -  H : the final value of the empirical entropy.
{
	n <- dim(X)[1]
	q <- dim(X)[2]
	# Rotate (permute) the data to the starting
	# point:
	X <- X[,c(start,(1:q)[-start])]
	tmp <- GradEmpEntropy(X,d,h)
	H <- tmp$H
	C <- tmp$C
	delta <- sum(C^2)
	Hv <- rep(0,maxiter+1)
	iter <- 0
	Hv[1] <- H
	if (display){
		if (d <= 2){
			par(cex=1.2,mai=c(1,1,0.5,0.01))
			plot(X[,1],X[,2],
				xlab=expression(italic(y[list(i,1)])),
				ylab=expression(italic(y[list(i,2)])),
				main=paste('H = ',round(H,7),
					', delta = ',round(delta,abs(log10(delta0))),
					sep=''))
		}else{
			pairs(X[,1:3],labels=c('y_1','y_2','y_3'),
				 main=paste('H = ',round(H,7),
				 		   ', delta = ',round(delta,abs(log10(delta0))),
				 		   sep=''))
		}
	}
	while (delta >= delta0 & iter < maxiter){
		tmp <- svd(C)
		W <- tmp$u
		sigma <- tmp$d
		V <- tmp$v
		U <- Exp(V,sigma,W)
		Xtmp <- X %*% t(U)
		Htmp <- EmpEntropy(Xtmp[,1:d],h)
		iter2 <- 0
		while (H - Htmp < delta/3 & iter2 < 20){
			iter2 <- iter2+1
			delta <- delta/2
			sigma <- sigma/2
			U <- Exp(V,sigma,W)
			Xtmp <- X %*% t(U)
			Htmp <- EmpEntropy(Xtmp[,1:d],h)
		}
		X <- Xtmp
		tmp <- GradEmpEntropy(X,d,h)
		H <- tmp$H
		C <- tmp$C
		delta <- sum(C^2)
		iter <- iter+1
		Hv[iter+1] <- H
		if (showsteps){
			tmp <- locator(1)
			if (d <= 2){
				plot(X[,1],X[,2],
					xlab=expression(italic(y[list(i,1)])),
					ylab=expression(italic(y[list(i,2)])),
					main=paste('H = ',round(H,7),
						', delta = ',round(delta,abs(log10(delta0))),
						sep=''))
			}else{
				pairs(X[,1:3],labels=c('y_1','y_2','y_3'),
					 main=paste('H = ',round(H,7),
					 		   ', delta = ',round(delta,abs(log10(delta0))),
					 		   sep=''))
			}
		}
	}
	if (display){
		if (d <= 2){
			plot(X[,1],X[,2],
				xlab=expression(italic(y[list(i,1)])),
				ylab=expression(italic(y[list(i,2)])),
				main=paste('H = ',round(H,7),
					', delta = ',round(delta,abs(log10(delta0))),
					sep=''))
		}else{
			pairs(X[,1:3],labels=c('y_1','y_2','y_3'),
				 main=paste('H = ',round(H,7),
				 		   ', delta = ',round(delta,abs(log10(delta0))),
				 		   sep=''))
		}
	}
	print(paste('After',iter,'iterations:'))
	print(paste("Empirical entropy = ",H," (h = ",h,")",sep=''))
	return(X)
}


GlobalPP1_rand <- function(X,h=0.5,mcsim=10)
	# Given a data matrix X with n rows and q > 2 colums,
	# this procedure computes the empirical entropy for
	# the q standard 1-dim. projections of X and of mcsim-1
	# random orthogonal transformations (row-wise) of X.
	# The output is an orthogonal transformation of
	# X such that the projection onto the first variable
	# has the smallest entropy seen in all mcsim*q
	# projections which have been evaluated.
{
	mcsim <- max(mcsim,2)
	q <- dim(X)[2]
	H <- Inf
	for (s in 1:mcsim){
		print(paste("Round",s))
		if (s > 1){
			U <- svd(matrix(rnorm(q*q),q,q))$u
			X <- Xtrans %*% U
		}
		for (i in 1:q){
			sigmatmp <- c(i,(1:q)[-i])
			Xtmp <- X[,sigmatmp]
			res <- GradEmpEntropy(Xtmp,d=1,h)
			if (res$H < H){
				H <- res$H
				Xtrans <- Xtmp
			}
		}
	}
	print(paste("Empirical entropy = ",H," (h = ",h,")",sep=''))
	return(Xtrans)
}


GlobalPP2_pre <- function(X,h=0.5,NrOfMinima=3)
	# Given a data matrix X with n rows and q > 2 colums,
	# this procedure computes the empirical entropy and
	# its gradient for the choose(q,2) possible standard
	# 2-dim. projections.
	# The output PartTable provides a list of
	# promising starting points for a local PP.
{
	q <- dim(X)[2]
	Table <- matrix(0,q*(q-1)/2,6)
	dimnames(Table)[[2]] <- c('i','j',
							  'H.init','dir.deriv.',
							  'R(H.init)','R(dir.deriv)')
	k <- 0
	for (i in 1:(q-1)){
		print(paste('(i,.) = (',i,',.)',sep=''))
		for (j in (i+1):q){
			k <- k+1
			sigmatmp <- c(i,j,(1:q)[-c(i,j)])
			Xtmp <- X[,sigmatmp]
			res <- GradEmpEntropy(Xtmp,d=2,h)
			Table[k,1:2] <- c(i,j)
			Table[k,3] <- res$H
			Table[k,4] <- -sum(res$C^2)
		}
	}
	Table[,5] <- rank(Table[,3])
	Table[,6] <- rank(Table[,4])
	tmp <- (Table[,5] <= NrOfMinima | Table[,6] <= NrOfMinima)
	Table <- Table[tmp,]
	return(Table)
}


GlobalPP2_rand <- function(X,h=0.5,mcsim=10)
	# Given a data matrix X with n rows and q > 2 colums,
	# this procedure computes the empirical entropy
	# for the choose(q,2) possible standard 2-dim. projections
	# of X and of mcsim-1 random orthogonal transformations
	# (row-wise) of X.
	# The output is an orthogonal transformation of X such
	# that the projection onto the first 2 variables has
	# the smallest entropie seen in all mcsim*choose(q,2)
	# projections which have been evaluated.
{
	mcsim <- max(mcsim,2)
	q <- dim(X)[2]
	H <- Inf
	for (s in 1:mcsim){
		print(paste("Round",s))
		if (s > 1){
			U <- svd(matrix(rnorm(q*q),q,q))$u
			X <- Xtrans %*% U
		}
		for (i in 1:(q-1)){
			print(paste('(i,.) = (',i,',.)',sep=''))
			for (j in (i+1):q){
				sigmatmp <- c(i,j,(1:q)[-c(i,j)])
				Xtmp <- X[,sigmatmp]
				res <- GradEmpEntropy(Xtmp,d=2,h)
				if (res$H < H){
					H <- res$H
					Xtrans <- Xtmp
				}
			}
		}
	}
	print(paste("Empirical entropy = ",H," (h = ",h,")",sep=''))
	return(Xtrans)
}

GlobalPP3_pre <- function(X,h=0.5,NrOfMinima=3)
	# Given a data matrix X with n rows and q colums,
	# this procedure computes the empirical entropy
	# and its gradient for the choose(q,3) possible
	# 3-dim. standard projections.
	# The output PartTable provides a list of
	# promising starting points for a local PP.
{
	q <- dim(X)[2]
	Table <- matrix(0,q*(q-1)*(q-2)/6,7)
	dimnames(Table)[[2]] <- c('i','j','k',
							  'H.init','dir.deriv.',
							  'R(H.init)','R(dir.deriv)')
	ell <- 0
	for (i in 1:(q-2)){
		for (j in (i+1):(q-1)){
			print(paste('(i,j,.) = (',i,',',j,',.)',sep=''))
			for (k in (j+1):q){
				ell <- ell+1
				sigmatmp <- c(i,j,k,(1:q)[-c(i,j,k)])
				Xtmp <- X[,sigmatmp]
				res <- GradEmpEntropy(Xtmp,d=3,h)
				Table[ell,1:3] <- c(i,j,k)
				Table[ell,4] <- res$H
				Table[ell,5] <- - sum(res$C^2)
			}
		}
	}
	Table[,6] <- rank(Table[,4])
	Table[,7] <- rank(Table[,5])
	tmp <- (Table[,6] <= NrOfMinima | Table[,7] <= NrOfMinima)
	Table <- Table[tmp,]
	return(Table)
}

GlobalPP3_rand <- function(X,h=0.5,mcsim=10)
	# Given a data matrix X with n rows and q > 3 colums,
	# this procedure computes the empirical entropy
	# and its gradient for the choose(q,3) possible
	# 3-dim. standard projections of X and mcsim-1 random
	# orthogonal transformations (row-wise) of X.
	# The output is an orthogonal transformation of X such
	# such that the projection onto the first 3 variables
	# has the smallest entropie seen in all mcsim*choose(q,3)
	# projections which have been evaluated.
{
	q <- dim(X)[2]
	H <- Inf
	for (s in 1:mcsim){
		print(paste("Round",s))
		if (s > 1){
			U <- svd(matrix(rnorm(q*q),q,q))$u
			X <- Xtrans %*% U
		}
		for (i in 1:(q-2)){
			for (j in (i+1):(q-1)){
				print(paste('(i,j,.) = (',i,',',j,',.)',sep=''))
				for (k in (j+1):q){
					sigmatmp <- c(i,j,k,(1:q)[-c(i,j,k)])
					Xtmp <- X[,sigmatmp]
					res <- GradEmpEntropy(Xtmp,d=3,h)
					if (res$H < H){
						H <- res$H
						Xtrans <- Xtmp
					}
				}
			}
		}
	}
	print(paste("Empirical entropy = ",H," (h = ",h,")",sep=''))
	return(Xtrans)
}


# Auxiliary functions ----

RecoverTransformation <- function(X,Y)
	# Input are data matrices X (n x q) and Y (n x d),
	# where d <= q and q < n. (X could also be a data frame.)
	# Under the assumption that each observation in Y
	# is an affine function of the corresponding
	# observation in X, this procedure recovers the
	# vector a in R^d and the matrix B in R^{q x d}
	# such that
	#    Y[i,] = a + t(B) %*% X[i,]
	# for 1 <= i <= n.
	# The third output argument is the resulting
	# root-mean-squared approximation error, just for
	# detecting situations in which Y is NOT such a
	# function of X.
{
	if (is.data.frame(X)){
		X <- as.matrix(X)
	}
	n <- dim(X)[1]
	q <- dim(X)[2]
	if (is.matrix(Y)){
		d <- dim(Y)[2]
	}else{
		d <- 1
	}
	C <- qr.solve(cbind(1,X),Y)
	if (d > 1){
		a <- C[1,]
		B <- C[2:(q+1),]
	}else{
		a <- C[1]
		B <- C[2:(q+1)]
	}
	Ytest <- rep(1,n)%*%t(a) + X %*% B
	rmse <- sqrt(mean(rowSums((Y - Ytest)^2)))
	return(list(a=a,B=B,rmse=rmse))
}


# Auxiliary function for numerical examples
# (to hide interesting structures in components 1 or 2
#  in p-dimensional space)

HaarMatrix <- function(p)
	# Computes an orthogonal pxp matrix such that
	# its first column is constant and its second
	# column has only values in {0, +- c} with
	# c = sqrt(1/floor(p/2)).
{
	if (p == 2){
		tmp <- sqrt(0.5)
		return(matrix(c(tmp,tmp,tmp,-tmp),2,2))
	}
	if (p == 3){
		tmp0 <- sqrt(1/3)
		tmp1 <- sqrt(0.5)
		tmp2 <- sqrt(0.5/3)
		return(matrix(c(tmp0,tmp0,tmp0,
						tmp1,0,-tmp1,
						-tmp2,2*tmp2,-tmp2),3,3))
	}
	M <- matrix(0,p,p-1)
	p2 <- p %/% 2
	tmp <- sqrt(1/(2*p2))
	M[1:p2,1] <- tmp
	M[(p-p2+1):p,1] <- -tmp
	Tmp <- HaarMatrix(p2)[,2:p2]
	if (p != 2*p2){
		tmp <- sqrt(1/p/(p-1))
		M[,2] <- -tmp
		M[p2+1,2] <- (p-1)*tmp
		M[1:p2,2 + (1:(p2-1))] <- Tmp
		M[p2 + 1 + (1:p2),1 + p2 + (1:(p2-1))] <- Tmp
	}else{
		M[1:p2,1 + (1:(p2-1))] <- Tmp
		M[p2 + (1:p2),p2 + (1:(p2-1))] <- Tmp
	}
	return(cbind(sqrt(1/p),M))
}
