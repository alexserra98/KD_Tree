#include <iostream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <memory.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <chrono>
#include <random>
#include <fstream>
#include </home/alexserra98/uni/HPC/KD_Tree/src/funct_mpi.hpp>
int main(int argc, char* argv[]){
    int mydim = std::stoi(argv[1]);
    std::fstream myfile;
    myfile.open("data",std::ios::out);
    if(!myfile){
        std::cout<<"File couldn't be opened"<<std::endl;
    }
    else{
        std::default_random_engine generator;
        std::uniform_real_distribution<TYPE> distribution(0.0,100.0);
        for(std::size_t i{0}; i<mydim; ++i){
            myfile << TYPE(distribution(generator))<<"\t"<<TYPE(distribution(generator))<<std::endl;
        }
        myfile.close();
    }
    

}
