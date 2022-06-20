The file "ProjectionPursuit.R" contains procedures for the combination of independent coordinate selection with projection pursuit, as described in
   L. Dümbgen, K. Gysel and F. Perler:
   Refining Invariant Coordinate Selection via Local Projection Pursuit.
   https://arxiv.org/abs/2112.11998.

The goal is to find interesting d-dimensional projections of a q-dimensional data set, where q > d, and d is typically 2 or 3, although other values are possible as well.

The R script "ProjectionPursuit_Demo.R" illustrates the use of the procedures.

For the case of 3-dimensional projections, a procedure Interactive3D() is used, see "Interactive3D.R". Its use is explained in the R scipt file "Interactive3D_Demo.R".

Lutz Duembgen, June 20, 2022.