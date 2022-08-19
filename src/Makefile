CXX=g++
CC=mpic++
CXXFLAGS = -DTYPE=float
CXXMPIFLAGS = -DMPITYPE=MPI_FLOAT
CXXOMPFLAGS = -fopenmp

all:  mpi omp
omp: main_openmp.cpp funct_openmp.hpp
	$(CXX) $(CXXFLAGS) -fopenmp main_openmp.cpp funct_openmp.hpp -o main_mpi.x
mpi: main_mpi.cpp funct_mpi.hpp
	$(CC) $(CXXFLAGS) $(CXXMPIFLAGS)  main_mpi.cpp funct_mpi.hpp -o main_mpi.x
clean:
	rm -f omp mpi
