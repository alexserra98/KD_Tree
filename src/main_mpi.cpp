#include <iostream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <memory.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include "funct_mpi.hpp"
#include <chrono>
#include <mpi.h>
#include <exception>
#include <random>

#define COUNT 10
#define xdim 10
#define ydim 10
#define ndim 2
#define MAX_LENGTH 30
//#define DEBUG 0
//#define DEBUG_b 0
//#define DEBUG_Container 0	   
//#define DEBUG_r0 0
//#define DEBUG_2
//#define DEBUG_cm


int main(int argc, char* argv[]){
   
  double begin1, end1, tbegin,tend, elapsed1,telapsed; 
  int dim_tot;
  //dim_tot = std::stoi(argv[1]); 
  std::vector<point<TYPE>> mydata;



  MPI_Datatype MPI_point; // custom datatype
  MPI_Datatype MPI_node; // custom datatype
  
  std::vector<int> list_rank; // order of ranks for updating indexes
  list_rank.push_back(0);   
  MPI_Init ( NULL, NULL );
    /* dataset creation */
  tbegin =MPI_Wtime();
  dataset_creation(mydata);
  tend = MPI_Wtime();
  telapsed = tend - tbegin;
  dim_tot = mydata.size();
  /* MPI Vars */
  MPI_Status status;
  int rank;
  int size;
  
  

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  
  if(rank==0)   begin1 = MPI_Wtime();

  /* check size */
  try{
    check_size_proc(size);
  }
  catch(const std::exception& e){
    std::cerr << e.what() << '\n';
    return -2;
  }
  
  try{
    check_size_dataset(size,dim_tot);
  }
  catch(const std::exception& e){
    std::cerr << e.what() << '\n';
    return -3;
  }


  const int lengths[3] = { 2,1,1 };

  /* Building custom types */
  MPI_Type_contiguous(2, MPITYPE,&MPI_point);
  MPI_Type_commit(&MPI_point);
  MPI_Aint displacements[3];
  node<float> dummy;
  MPI_Aint base_address;
  MPI_Get_address(&dummy, &base_address);    
  MPI_Get_address(&dummy.p, &displacements[0]);
  MPI_Get_address(&dummy.right, &displacements[1]);
  MPI_Get_address(&dummy.left, &displacements[2]);
  displacements[0] = MPI_Aint_diff(displacements[0], base_address);
  displacements[1] = MPI_Aint_diff(displacements[1], base_address);
  displacements[2] = MPI_Aint_diff(displacements[2], base_address);
  MPI_Datatype types[3] ={MPI_point,MPI_INT,MPI_INT};
  MPI_Type_create_struct(1, lengths, displacements,types , &MPI_node);
  MPI_Type_commit(&MPI_node);
  

  MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
  MPI_Request request;
  
  /* vectors init */
  std::vector<point<TYPE>> mypile; // list of nodes of each process
  std::vector<node<TYPE>> serial_tree; // first levels of the tree computed in parallel
  std::vector<node<TYPE>> parallel_tree; // serial tree calculated by each tree
  std::vector<std::pair<int,int>> dim_tree; // tree which map send and receive dimensions
  
  int mystride[size];
  int boundaries[size];
  
  mydata.reserve(dim_tot+1);   

  /* position in the tree */
  int level,axis,mydim, levels;
  position_in_the_tree(levels,level,axis,mydim,dim_tot,rank,size);

  /*init stride */
  for(auto i{0};i<size;++i) mystride[i]=0;


  /* filling dim_tree */
  dim_tree.resize(pow(2,levels+1)-1);
  dimtree_fill(dim_tree,mystride,boundaries, 0,dim_tot,levels,-1,0,-1,list_rank);  
  #ifdef DEBUG_dimtree
  std::cout<<"size: " <<dim_tree.size()<<std::endl;
  if(rank==0){
  for(auto i: dim_tree){
       std::cout<<"#dimtree "<<i.first<<" , "<<i.second<<std::endl;
  }
  }
  #endif

  #ifdef DEBUG
  std::cout<<"level = "<<level<<" axis = "<<axis<<" rank = "<<rank<<std::endl;
  #endif

  /* vars */

  int dim_c; //
  int dim_n; // dimension of serial tree
  int myindex = -1; // upperlimit of portion on mydata in which the process operates, default = -1 aka last element in the vector
  int tmpind;
  int counter = 0; // counter used in while section of Building tree
  int stride; // offset in the recv vector in the Sending Back section
  int parent_index; // index of the node to which attach the recv tree
  int dim_r; // recv dimension 
  int dim_s[2]; // send dimension
  int flag=-1; // flag for while of dim search; 
 
  


  /******************************************
   * Building tree                          *
   ******************************************/
  
  if(rank==0){
    mydim=dim_tot;
    mydata.resize(mydim);
    parallel_tree.reserve(mydim+1);
    /* dataset creation */
    /*
    tbegin =MPI_Wtime();
    dataset_creation(mydata,mydim);
    tend = MPI_Wtime();
    telapsed = tend - tbegin;
    */
  }

  else{ 
    auto tmp_dim=levels-floor(log2(rank));
    parallel_tree.reserve(dim_tot+1);
  }
  
  
  
  while(counter<=levels){
    if(rank<pow(2,counter)){ //only nodes occurring in that level can work
      int flag = pow(2,counter)-1;
      while(flag!=-1){ // searching for dimension in dim_tree
        if(dim_tree[flag].first==rank){
          dim_r=dim_tree[flag].second;
          
	  int tmp_pos=pow(2,counter);
          dim_s[0]=dim_tree[2*flag+1].second; // 2^(level+1)+2*(pos_in_the_level) = 2*flag (+1 beacause is 0-indexed)
          dim_s[1]=dim_tree[2*flag+2].second;
	    #ifdef DEBUG
	  std::cout<<"rank "<<rank<<" mydata size "<<mydata.size()<<" I'm sending msg 1 of size "<<dim_s[0]<<" msg 2 of size "<<dim_s[1]<<" I'm recv msg of size "<<dim_r<<std::endl;
          #endif
	  flag=-2;
        }
	++flag;
      }      

      #ifdef DEBUG
      std::cout<<"rank "<<rank<<" myindex "<<myindex<<" mydim "<<mydim<<std::endl;
       #endif
      int son = rank + pow(2,counter);
      if(son<=size-1){ // if son rank don't exceed the number of process at disposal
        compute_median(parallel_tree,mydata,axis,dim_r-1);
        myindex = dim_s[0];
        #ifdef DEBUG
        std::cout<<"rank "<<rank<<" is sending to "<<son<<std::endl;
        #endif
        MPI_Isend(&mydata[myindex+1], dim_s[1],MPI_point,son,son,MPI_COMM_WORLD,&request);
      }
      else{ // proceed serially 
        dim_n = dim_r;
        myindex = (myindex-1)/2;

	      #ifdef DEBUG_Container
	      std::cout<<" rank: "<<rank<<" dim_r : "<<dim_r<<std::endl;
        #endif
	auto size_p=parallel_tree.size();
	if(size_p!=0)  parallel_tree[size_p-1].left = size_p;
        compute_median_s(parallel_tree,mydata,ndim,axis,0,dim_r-1);
      }

      #ifdef DEBUG
      std::cout<<"rank "<<rank<<" son "<<son<<" myindex "<<myindex<< " counter"<<counter<<std::endl;
      #endif

    }


      #ifdef DEBUG
          std::cout<<"rank "<<rank<<" mydata size "<<mydata.size()<<" I'm sending msg 1 of size "<<dim_s[0]<<" msg 2 of size "<<dim_s[1]<<" I'm recv msg of size "<<dim_r<<std::endl;
      #endif
	  if(rank>=pow(2,counter)&&rank<pow(2,counter+1)){ // only the right sons are recv
      int flag1 = pow(2,counter+1)-1;
      while(flag1!=-1){ // searching for dimension in dim_tree
        if(dim_tree[flag1].first==rank){
          dim_r=dim_tree[flag1].second;
          
	  int tmp_pos=pow(2,counter);
          dim_s[0]=dim_tree[2*flag1+1].second; // 2^(level+1)+2*(pos_in_the_level) = 2*flag (+1 beacause is 0-indexed)
          dim_s[1]=dim_tree[2*flag1+2].second;
	    #ifdef DEBUG
	  std::cout<<"rank "<<rank<<" mydata size "<<mydata.size()<<" I'm sending msg 1 of size "<<dim_s[0]<<" msg 2 of size "<<dim_s[1]<<" I'm recv msg of size "<<dim_r<<std::endl;
          #endif
	  flag1=-2;
        }
	++flag1;
      }      

      	int exp=int(log2(rank));
        int parent= rank-pow(2,exp);

        #ifdef DEBUG
        std::cout<<"rank "<<rank<<" parent "<<parent<<" counter"<<counter<<std::endl;
        #endif

        MPI_Recv(mydata.data(), dim_r+1,MPI_point,parent,rank,MPI_COMM_WORLD, &status);          
        dim_c=dim_r;
     

      }
     
    MPI_Barrier(MPI_COMM_WORLD);
    ++counter;
  }

  /*merging parallel_tree with serial*/
  
  auto size_p=parallel_tree.size();
  

  #ifdef DEBUG_Container
  //for(std::size_t i{0}; i<mypile.size();++i) std::cout<<"rank "<<rank<<" elem "<<i<<": "<<mypile[i].data[0]<<","<<mypile[i].data[1]<<std::endl;
  std::cout<<" rank: "<<rank<<" mysize: "<<parallel_tree.size()<<std::endl;
  //print2D(0, serial_tree);
  //print2D(0, parallel_tree);
  for(std::size_t i{0}; i<dim_tot; ++i) std::cout<<"rank: "<<rank<<" , "<<"i: "<<i<<" , "<<parallel_tree[i].p.data[0]<<" , "<<parallel_tree[i].p.data[1]<<std::endl;
  #endif  
  

  /*************************************
   * Sending it back                   *
   ************************************/
#ifdef LASTBUG
  std::cout<<"here i am"<<std::endl;
#endif  
  counter=levels; // the while start from the leaves up to the root

  if(rank==0){stride = dim_n + levels ; parent_index=levels-1;}
  else {stride = dim_n + levels-1 - floor(log2(rank)); parent_index=levels -2- floor(log2(rank));}
 
  dim_r = dim_n;
  int dim_ss = dim_n;
  int ierr;
  #ifdef DEBUG_b
  //  std::cout<<
  #endif
  flag=-1;
  int flag1=-1;
  while(counter>0){
      int son = rank + pow(2,counter-1);
     
	 	  if(rank>=pow(2,counter-1)&&rank<pow(2,counter)){ // rank has to send

	 flag = pow(2,counter)-1;
	 while(flag!=-1){ // searching for dimension in dim_tree
	    if(dim_tree[flag].first==rank){
	      dim_ss=dim_tree[flag].second;
	      int tmp_pos=pow(2,counter);
	      dim_r=dim_tree[2*flag+2].second; // 2^(level+1)+2*(pos_in_the_level) = 2*flag (+1 beacause is 0-indexed)
	      flag=-2;
	      }
	    ++flag;
	 }
 #ifdef DEBUG_b
      std::cout<<"counter: "<<counter <<" rank: "<<rank<<" stride: "<<stride<<" dim_r: "<<dim_r<<" dim_ss: "<<dim_ss<<std::endl;
      #endif

          int exp=int(log2(rank));
          int parent= rank-pow(2,exp); 
	  
	  ierr=MPI_Send(parallel_tree.data(),dim_ss, MPI_node, parent, parent, MPI_COMM_WORLD);
	 
          #ifdef DEBUG_b
          std::cout<<"IERR  "<<ierr<<std::endl;
	  std::cout<<"rank "<<rank<< " is sending to  parent "<<parent<<std::endl;
          #endif

          } 


	 if(rank<pow(2,counter-1)){  
	   // if(rank==0||rank==1||rank==2||rank==3){
	   	 flag1 = pow(2,counter-1)-1;
	 while(flag1!=-1){ // searching for dimension in dim_tree
	    if(dim_tree[flag1].first==rank){
	      dim_ss=dim_tree[flag1].second;
	      int tmp_pos=pow(2,counter);
	      dim_r=dim_tree[2*flag1+2].second; // 2^(level+1)+2*(pos_in_the_level) = 2*flag (+1 beacause is 0-indexed)
	      flag1=-2;
	      }
	    ++flag1;
	 } 
 #ifdef DEBUG_b
      std::cout<<"counter: "<<counter <<" rank: "<<rank<<" stride: "<<stride<<" dim_r: "<<dim_r<<" dim_ss: "<<dim_ss<<std::endl;
      #endif

	   MPI_Recv(&parallel_tree[stride], dim_r+1, MPI_node, son, rank, MPI_COMM_WORLD,&status);
            #ifdef DEBUG_b
            std::cout<<"rank "<<rank<<" is receiving from  son "<<son<<"a msg of size "<<dim_r<<" on counter "<<counter <<std::endl;
            #endif

            parallel_tree[parent_index].right = stride; // linking left branch to main chain
            --parent_index;

	     stride+=dim_r;

          }
    --counter;
  
   } 


  if(rank==0){
        #ifdef DEBUG_b
    for(auto h{0}; h<8;++h){
      std::cout<<"b: "<<boundaries[h]<<std::endl;
    }
    #endif

    int r=0;
    int stride1=0;
    for(int i=0; i<size-1; ++i){
      int tmp_rank=list_rank[i];
      stride1+=mystride[tmp_rank];
      r=mystride[list_rank[i+1]];
      for(int k=stride1; k<(r+stride1); ++k){
	if(parallel_tree[k].right!=-1) parallel_tree[k].right+=stride1;
        if(parallel_tree[k].left!=-1) parallel_tree[k].left+=stride1;
      }
      
      }
    
    #if OUTPUT==1
    //for(std::size_t i{0}; i<dim_tot; ++i) std::cout<<"{data: "<<parallel_tree[i].p.data[0]<<" , "<<parallel_tree[i].p.data[1]<<"; right: "<<parallel_tree[i].right<<" left: "<<parallel_tree[i].left<<"}" <<std::endl;
     if(dim_tot<100) print2D(0, parallel_tree);
    #endif
     end1 = MPI_Wtime();
     elapsed1 = end1 - begin1;
     auto time_elapsed = elapsed1-telapsed;
    std::cout<<size<<","<<time_elapsed<<std::endl;
  }

  
  MPI_Finalize();

  return 0;
}

