CXX=g++
CC=mpic++
TYPE = float
MPITY = MPI_FLOAT
OUTPUT = 0
CXXTYPE = -DTYPE=$(TYPE)
CXXMPITYPE = -DMPITYPE=$(MPITY)
CXXOMPFLAGS = -fopenmp
CXXOUTPUT = -DOUTPUT=$(OUTPUT)

all:  mpi omp
omp: main_openmp.cpp funct_openmp.hpp
	$(CXX) $(CXXTYPE) $(CXXOMPFLAGS) $(CXXOUTPUT) main_openmp.cpp funct_openmp.hpp -o main_omp.x
mpi: main_mpi.cpp funct_mpi.hpp
	$(CC) $(CXXTYPE) $(CXXMPITYPE) $(CXXOUTPUT)  main_mpi.cpp funct_mpi.hpp -o main_mpi.x
clean:
	rm -f omp mpi
