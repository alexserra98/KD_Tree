# Parallel KD_Tree
This repository contains two implementations of a parallel kd-tree using MPI and OpenMP in C++11. The work was written as an assignment of FHPC course 2021-2022 @ DSSC, Units.
## Introduction
A K-dimensional tree is a data structure widely used for partitioning and organizing
points in a k-dimensional space, they’re involved in many different applications
such as searches involving a multidimensional search key (e.g. range
searches and nearest neighbor searches). The task of the assignement was to write a program which takes as input a 2d dataset and compute over it the associated kd-tree.
The algorithm conist essentialy in these three steps
- picking round-robin an axis
- select the median point over that axis
- reiterate on the two halves at the left and right of the selected point
The algorithm stops when the region of the spaces that the program is parsing contains only one point. 
## Parallelization
### MPI
The parallelization of the serial algorithm through the MPI interface consists
in distributing the recursive calls among the processes provided. The communication among processess is organized in the following way: 
starting from the master process - rank 0 in our implementation - each worker
keep half of the data set received and the send the other half to his son, once
all the processors has received their input data they proceed the computation
serially. At the end of this section each processor, with exception of the master,
starting from the bottom of the tree, start sending to his parent the tree computed
on the received portion of the data set. The parent, on his behalf, will
merge the received buffer with his tree so that eventually the master process
will have the complete tree.
### OMP 
The OpenMP implementation is similar to the MPI one although it follows a
much simpler scheme. The first call of the function is made by master thread
and then each recursive call is assigned to a new task.
## Compilation
The src folder contains a Makefile; you can specify the parallalelization framework in the following way:    
```
    make omp/mpi
```   

In addition you can specify the data type of the point of the dataset with the following flags:
| Parameter | Value           | Default   | 
|-----------|-----------------|-----------|          
| TYPE      | float/double    | float     |
| MPITY     | MPI_FLOAT/DOUBLE| MPI_FLOAT |

## Run
