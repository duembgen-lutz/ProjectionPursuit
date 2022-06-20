Interactive3D <- function(X0,
	pchX=1,colX='black',
	col.outside='gray',
	rmax=1.1*max(sqrt(rowSums(X0[,1:3]^2))),
	observer=3.5,pch=1,
	center.X=TRUE,
	Edges=NULL,
	col.edges='gray',
	show.axes=TRUE,
	col.axes='blue',
	names.axes=1:3,
	AddAxes0=NULL,
	col.addaxes='red',
	names.addaxes=NULL,
	static=FALSE)
	# For a data matrix X0 with at least 3 columns, this
	# procedure visualizes the empirical distribution of
	# the first three columns in a 3D rotation plot.
	# By clicking on the plot, you can turn the data set,
	# and clicking in the four corners, you can zoom in (ZI),
	# zoom out (ZO), look at the points from direction (1,1,1),
	# or terminate the procedure (End).
	# If the additional parameter Edges is given, its rows
	# are index pairs (i,j) specifying edges of a graph with
	# vertices X0[i,] to be displayed as well.
	# If show.axes==TRUE, the three main coordinate axes are
	# displayed in color col.axes.
	# If the additional parameter AddAxes0 is given (e.g. for
	# BiPlots in PCA), the first three columns of AddAxes0 are
	# properly rescaled and displayed as well, in color
	# col.addaxes.
	#
	# The ouptput is a data matrix Xproj with three columns,
	# representing the last view of the data. By means of
	# RecoverTransformation(X0,Xproj), one can find out how
	# to generate Xproj from the original data matrix X0 via an affine transformation.
	# 
	# Lutz Duembgen, March 8, 2022
{
	n <- dim(X0)[1]
	q <- dim(X0)[2]
	if (length(colX)==1){
		colX=rep(colX,n)
	}
	if (length(pchX)==1){
		pchX=rep(pchX,n)
	}
	if (q < 3){
		return("Too few variables (>= 3 are needed)!")
	}
	X1 <- X0[,1:3]
	if (center.X){
		xbar <- colMeans(X1)
		X1 <- X1 - rep(1,n) %*% t(xbar)
	}else{
		xbar <- rep(0,3)
	}
	Axes1 <- rbind(diag(3),xbar)
	if (!is.null(AddAxes0)){
		if (is.matrix(AddAxes0) && dim(AddAxes0)[1] > 1){
			rA <- max(sqrt(rowSums(AddAxes0[,1:3]^2)))
			Axes1 <- rbind(Axes1,AddAxes0[,1:3]/rA)
			nraddaxes <- dim(AddAxes0)[1]
		}else{
			Axes1 <- rbind(Axes1,
						  AddAxes0[1:3]/
						  	sqrt(sum(AddAxes0[1:3]^2)))
			nraddaxes <- 1
		}
		if (is.null(names.addaxes)){
			names.addaxes <- 1:nraddaxes
		}
	}else{
		nraddaxes <- 0
	}
	X <- X1/rmax
	Inside <- (rowSums(X^2) <= 1)
	Outside <- !Inside
	Axes <- Axes1
	par(cex=1,mai=c(0.01,0.01,0.01,0.01))
	plot(0,0,col='white',
		 xlab='',ylab='',
		 xlim=c(-1,1),ylim=c(-1,1))
	magfac <- 1 - X[,3]/observer
	if (!is.null(Edges)){
		for (i in 1:dim(Edges)[1]){
			lines(X[Edges[i,],1]/magfac[Edges[i,]],
				  X[Edges[i,],2]/magfac[Edges[i,]],
				  col=col.edges)
		}
	}
	tmp <- Outside & (magfac > 0)
	if (!is.character(pchX)){
		points(X[tmp,1]/magfac[tmp],X[tmp,2]/magfac[tmp],
			   pch=pchX[tmp],col=col.outside,
			   cex=pmin(1/magfac[tmp],1/(1 - 1/observer)))
		points(X[Inside,1]/magfac[Inside],
			   X[Inside,2]/magfac[Inside],
			   pch=pchX[Inside],col=colX[Inside],
			   cex=1/magfac[Inside])
	}else{
		if (sum(tmp)>0){
			text(X[tmp,1]/magfac[tmp],
				 X[tmp,2]/magfac[tmp],
				 pchX[tmp],
				 col=col.outside,
				 cex=pmin(1/magfac[tmp],1/(1 - 1/observer)))
		}
		if (sum(Inside)>0){
			text(X[Inside,1]/magfac[Inside],
				 X[Inside,2]/magfac[Inside],
				 pchX[Inside],
				 col=colX[Inside],
				 cex=1/magfac[Inside])
		}
	}
	if (show.axes){
		arrows(x0=rep(0,3),
			   x1=0.9*Axes[1:3,1]/
			   	(1 - 0.9*Axes[1:3,3]/observer),
			   y0=rep(0,3),
			   y1=0.9*Axes[1:3,2]/
			   	(1 - 0.9*Axes[1:3,3]/observer),
			   length=0,
			   col=col.axes)
		text(Axes[1:3,1]/(1 - Axes[1:3,3]/observer),
			 Axes[1:3,2]/(1 - Axes[1:3,3]/observer),
			 names.axes,
			 cex=1/(1 - Axes[1:3,3]/observer),
			 col=col.axes)
		if (nraddaxes > 0){
			jj <- 4 + (1:nraddaxes)
			arrows(x0=rep(0,nraddaxes),
				   x1=0.9*Axes[jj,1]/
				   	(1 - 0.9*Axes[jj,3]/observer),
				   y0=rep(0,nraddaxes),
				   y1=0.9*Axes[jj,2]/
				   	(1 - 0.9*Axes[jj,3]/observer),
				   length=0,
				   col=col.addaxes)
			text(Axes[jj,1]/(1 - Axes[jj,3]/observer),
				 Axes[jj,2]/(1 - Axes[jj,3]/observer),
				 names.addaxes,
				 cex=1/(1 - Axes[jj,3]/observer),
				 col=col.addaxes)
		}
	}
	if (static){
		return()
	}
	abline(a=1.9,b=-1,col='gray')
	text(0.995,0.995,'ZI')
	abline(a=1.9,b=1,col='gray')
	text(-0.995,0.995,'ZO')
	abline(a=-1.9,b=-1,col='gray')
	text(-0.995,-0.995,'St')
	abline(a=-1.9,b=1,col='gray')
	text(0.995,-0.995,'End')
	tmp <- locator(1)
	while (tmp$y - tmp$x >= -1.9){
		if (tmp$x + tmp$y > 1.9){
			X <- X*2^(0.25)
			rmax <- rmax*2^(-0.25)
			Inside <- (rowSums(X^2) <= 1)
			Outside <- !Inside
		}
		if (tmp$y - tmp$x > 1.9){
			X <- X*2^(-0.25)
			rmax <- rmax*2^(0.25)
			Inside <- (rowSums(X^2) <= 1)
			Outside <- !Inside
		}
		if (tmp$x + tmp$y < -1.9){
			U <- matrix(c(-sqrt(1/2),-sqrt(1/6),sqrt(1/3),
						  sqrt(1/2),-sqrt(1/6),sqrt(1/3),
						  0,sqrt(2/3),sqrt(1/3)),3,3)
			X <- X1 %*% t(U)
			X <- X/rmax
			Axes <- Axes1 %*% t(U)
		}
		if (tmp$x + tmp$y <= 1.9 &
			tmp$y - tmp$x <= 1.9 &
			tmp$x + tmp$y >= -1.9){
			t <- sqrt(tmp$x^2 + tmp$y^2)
			a <- tmp$x/t
			b <- tmp$y/t
			U <- matrix(0,3,3)
			U[3,1] <- 1
			U[1:2,2] <- c(a,b)
			U[1:2,3] <- c(-b,a)
			theta <- atan(t/2)
			ct <- cos(theta)
			st <- sin(theta)
			R <- matrix(0,3,3)
			R[3,3] <- 1
			R[1:2,1] <- c(ct,-st)
			R[1:2,2] <- c(st,ct)
			X <- X %*% U
			X <- X %*% R
			X <- X %*% t(U)
			Axes <- Axes %*% U
			Axes <- Axes %*% R
			Axes <- Axes %*% t(U)
		}
		plot(0,0,col='white',
			 xlab='',ylab='',
			 xlim=c(-1,1),ylim=c(-1,1))
		magfac <- 1 - X[,3]/observer
		if (!is.null(Edges)){
			for (i in 1:dim(Edges)[1]){
				lines(X[Edges[i,],1]/magfac[Edges[i,]],
					  X[Edges[i,],2]/magfac[Edges[i,]],
					  col=col.edges)
			}
		}
		tmp <- Outside & (magfac > 0)
		if (!is.character(pchX)){
			points(X[tmp,1]/magfac[tmp],X[tmp,2]/magfac[tmp],
				   pch=pchX[tmp],col=col.outside,
				   cex=pmin(1/magfac[tmp],1/(1 - 1/observer)))
			points(X[Inside,1]/magfac[Inside],
				   X[Inside,2]/magfac[Inside],
				   pch=pchX[Inside],col=colX[Inside],
				   cex=1/magfac[Inside])
		}else{
			if (sum(tmp)>0){
				text(X[tmp,1]/magfac[tmp],
					 X[tmp,2]/magfac[tmp],
					 pchX[tmp],
					 col=col.outside,
					 cex=pmin(1/magfac[tmp],1/(1 - 1/observer)))
			}
			if (sum(Inside)>0){
				text(X[Inside,1]/magfac[Inside],
					 X[Inside,2]/magfac[Inside],
					 pchX[Inside],
					 col=colX[Inside],
					 cex=1/magfac[Inside])
			}
		}
		if (show.axes){
			arrows(x0=rep(0,3),
				   x1=0.9*Axes[1:3,1]/
				   	(1 - 0.9*Axes[1:3,3]/observer),
				   y0=rep(0,3),
				   y1=0.9*Axes[1:3,2]/
				   	(1 - 0.9*Axes[1:3,3]/observer),
				   length=0,
				   col=col.axes)
			text(Axes[1:3,1]/(1 - Axes[1:3,3]/observer),
				 Axes[1:3,2]/(1 - Axes[1:3,3]/observer),
				 names.axes,
				 cex=1/(1 - Axes[1:3,3]/observer),
				 col=col.axes)
			if (nraddaxes > 0){
				jj <- 4 + (1:nraddaxes)
				arrows(x0=rep(0,nraddaxes),
					   x1=0.9*Axes[jj,1]/
					   	(1 - 0.9*Axes[jj,3]/observer),
					   y0=rep(0,nraddaxes),
					   y1=0.9*Axes[jj,2]/
					   	(1 - 0.9*Axes[jj,3]/observer),
					   length=0,
					   col=col.addaxes)
				text(Axes[jj,1]/(1 - Axes[jj,3]/observer),
					 Axes[jj,2]/(1 - Axes[jj,3]/observer),
					 names.addaxes,
					 cex=1/(1 - Axes[jj,3]/observer),
					 col=col.addaxes)
			}
		}
		abline(a=1.9,b=-1,col='gray')
		text(0.995,0.995,'ZI')
		abline(a=1.9,b=1,col='gray')
		text(-0.995,0.995,'ZO')
		abline(a=-1.9,b=-1,col='gray')
		text(-0.995,-0.995,'St')
		abline(a=-1.9,b=1,col='gray')
		text(0.995,-0.995,'End')
		tmp <- locator(1)
	}
	plot(0,0,col='white',
		 xlab='',ylab='',
		 xlim=c(-1,1),ylim=c(-1,1))
	magfac <- 1 - X[,3]/observer
	if (!is.null(Edges)){
		for (i in 1:dim(Edges)[1]){
			lines(X[Edges[i,],1]/magfac[Edges[i,]],
				  X[Edges[i,],2]/magfac[Edges[i,]],
				  col=col.edges)
		}
	}
	tmp <- Outside & (magfac > 0)
	if (!is.character(pchX)){
		points(X[tmp,1]/magfac[tmp],X[tmp,2]/magfac[tmp],
			   pch=pchX[tmp],col=col.outside,
			   cex=pmin(1/magfac[tmp],1/(1 - 1/observer)))
		points(X[Inside,1]/magfac[Inside],
			   X[Inside,2]/magfac[Inside],
			   pch=pchX[Inside],col=colX[Inside],
			   cex=1/magfac[Inside])
	}else{
		if (sum(tmp)>0){
			text(X[tmp,1]/magfac[tmp],
				 X[tmp,2]/magfac[tmp],
				 pchX[tmp],
				 col=col.outside,
				 cex=pmin(1/magfac[tmp],1/(1 - 1/observer)))
		}
		if (sum(Inside)>0){
			text(X[Inside,1]/magfac[Inside],
				 X[Inside,2]/magfac[Inside],
				 pchX[Inside],
				 col=colX[Inside],
				 cex=1/magfac[Inside])
		}
	}
	if (show.axes){
		arrows(x0=rep(0,3),
			   x1=0.9*Axes[1:3,1]/
			   	(1 - 0.9*Axes[1:3,3]/observer),
			   y0=rep(0,3),
			   y1=0.9*Axes[1:3,2]/
			   	(1 - 0.9*Axes[1:3,3]/observer),
			   length=0,
			   col=col.axes)
		text(Axes[1:3,1]/(1 - Axes[1:3,3]/observer),
			 Axes[1:3,2]/(1 - Axes[1:3,3]/observer),
			 names.axes,
			 cex=1/(1 - Axes[1:3,3]/observer),
			 col=col.axes)
		if (nraddaxes > 0){
			jj <- 4 + (1:nraddaxes)
			arrows(x0=rep(0,nraddaxes),
				   x1=0.9*Axes[jj,1]/
				   	(1 - 0.9*Axes[jj,3]/observer),
				   y0=rep(0,nraddaxes),
				   y1=0.9*Axes[jj,2]/
				   	(1 - 0.9*Axes[jj,3]/observer),
				   length=0,
				   col=col.addaxes)
			text(Axes[jj,1]/(1 - Axes[jj,3]/observer),
				 Axes[jj,2]/(1 - Axes[jj,3]/observer),
				 names.addaxes,
				 cex=1/(1 - Axes[jj,3]/observer),
				 col=col.addaxes)
		}
	}
	Xproj <- rmax*X
	return(Xproj)
}


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
