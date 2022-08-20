#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <memory.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <chrono>
#include "funct_openmp.hpp"
#include "omp.h"
#define COUNT 10
#define xdim 10
#define ydim 10
#define ndim 2




int main(int argc, char* argv[]){

    int dim_tot, n;
    dim_tot =  atoi(argv[1]);
    n = atoi(argv[2]);
    auto tbegin = std::chrono::high_resolution_clock::now();
    std::vector<point<TYPE>> mydata;
    std::vector<node<TYPE>> parallel_tree;
    mydata.resize(dim_tot);
    parallel_tree.reserve(dim_tot+1);
    auto begin1 = std::chrono::high_resolution_clock::now();
    dataset_creation(mydata, dim_tot);
    auto end1 = std::chrono::high_resolution_clock::now();
    auto elapsed1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1);
    std::size_t start{0};
    std::size_t end{mydata.size()-1};
    int provax=0;
    
    #ifdef DEBUG
    int index = kthSmallest(mydata,0,dim_tot-1, 0,1);
    mydata[index].print_node();
    std::cout<<index<<std::endl;
    for(std::size_t i{0}; i<mydata.size(); ++i){
        mydata[i].print_node();
    }

    #endif
    #ifdef _OPENMP

        #pragma omp parallel num_threads(n)
        {
            #pragma omp master
            {
                compute_median_p(parallel_tree,mydata,2,0,start,end);
            }   
            #pragma omp barrier
            #pragma omp master
            {
                #if OUTPUT == 1
                print2D(0, parallel_tree);
		    #endif
            }
        } 
       #else
        compute_median_s(parallel_tree,mydata,2,0,start,end);
        #ifdef OUTPUT
	if(dim_tot<100)   print2D(0, parallel_tree);
	    #endif
    #endif   
    auto tend = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tbegin);
    auto tmp=(elapsed.count()- elapsed1.count())*1e-9;
    std::cout<<n<<"," <<tmp<<std::endl;

    return 0;
}
