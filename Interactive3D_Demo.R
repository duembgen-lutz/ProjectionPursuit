source("Interactive3D.R")

# To illustrate the use of Interactive3D(),
# generate a cube in dimension 3:
X3 <- rbind(c(1,1,1),
			c(1,1,-1),
			c(1,-1,1),
			c(-1,1,1),
			c(1,-1,-1),
			c(-1,1,-1),
			c(-1,-1,1),
			c(-1,-1,-1))
Edges <- rbind(c(1,2),
			   c(1,3),
			   c(1,4),
			   c(2,5),
			   c(2,6),
			   c(3,5),
			   c(3,7),
			   c(4,6),
			   c(4,7),
			   c(5,8),
			   c(6,8),
			   c(7,8))
# Look at the cube from different angles.
# In the end, try to restore the original
# view, with first axis pointing to the right,
# the second axis pointing up, the third axis
# pointing to you.
# By clicking on the plot, you can turn the data set,
# and clicking in the four corners, you can zoom in (ZI),
# zoom out (ZO), look at the points from direction (1,1,1),
# or terminate the procedure (End).
Interactive3D(X3,Edges=Edges,center.X=FALSE)
# The output "Xproj" is the point cloud (vertices)
# you produced...

# The same with mirrored main axes added:
Interactive3D(X3,Edges=Edges,center.X=FALSE,
			  AddAxes0 = -diag(3),names.addaxes = rep('',3))


# If a data analysis produced a particularly interesting
# view of the data, one can recover the corresponding
# linear transformation:

n <- 200
X <- matrix(rnorm(n*2),n,2)
Y <- X[,1]*X[,2]
XY <- cbind(X,Y)

XYb <- Interactive3D(XY)
trans <- RecoverTransformation(XY,XYb)
trans$a
trans$B
# XYb[i,] = trans$a + t(trans$B)%*%XY[i,]
Error <- XYb - rep(1,n)%*%t(trans$a) - XY%*%trans$B
sqrt(sum(Error^2))

