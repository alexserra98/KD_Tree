#pragma once
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
#define DEBUG2 1
#define COUNT 10
#define xdim 10
#define ydim 10
#define ndim 2
using  uint = unsigned int;

/******* child ranks for MPI call *******/
void child_nodes(uint n, uint kids[2]){
    uint level= floor(log2(n)); // level in which the node is located
    uint pos = n-pow(2,level); // position in the level
    kids[0]= pow(2,level+1)+2*pos; // sum of the nodes in previous levels +1 +shift due to the position
    kids[1]= pow(2,level+1)+1+2*pos;
}
/******* point implementation *******/
template<typename T>
struct point{
  T data[ndim];
  point() = default;
  point(const T& a, const T& b): data{a,b} {} 
  point(const T mydata[ndim])  {for(std::size_t i{0}; i<ndim; ++i) data[i] = mydata[i];} 
  point& operator=(const point<T>& mypoint){
    for(std::size_t i{0}; i<ndim; ++i) data[i] = mypoint.data[i]; 
    return *this; 
  }                                                                                                                                                   
  point(const point<T>& mypoint) noexcept  { // delegating ctor                                                                                                              
    *this=mypoint;                                                                                                                             
  }
  point(const point&& mypoint) noexcept  { // delegating ctor                                                                                                           
    *this=std::move(mypoint); 
  }

};


/******* node implementation *******/
template<typename T>
struct node {
  point<T> p;
  int right;
  int left;

  node() = default;
  node(const point<T>& mypoint): p(mypoint), right{-1}, left{-1} {}
  node(const point<T>& mypoint, const int r,const int l): p(mypoint), right{r}, left{l} {}
  
  node& operator=(const node<T>& mynode){
    p=mynode.p;
    right=mynode.right;
    left=mynode.left;
    return *this;
  }                                                                                                       
  node(const node<T>& mynode)  noexcept {
      *this=mynode;
  }
  node(const node<T>&& mynode) noexcept {
      *this=std::move(mynode);
  }
};
 
template <typename T>
int compute_median(std::vector<node<T>>& allc,std::vector<point<T>>& data,const int axis, int dim){
  int r;  
  r= dim;
  int l= 0;                                                                                          
  // find median                                                                                                                       
  int tmp=r/2;
  #ifdef DEBUG_cm
  std::cout<<"size"<<v.size()<<" r: "<<r<<" l: "<<l<<" v: "<<data[r].data[0]<<", "<<data[r].data[1]<<std::endl;
  std::cout<<"***********"<<std::endl;
  for(std::size_t i{0}; i<=r; ++i) std::cout<<"r: "<<r<<" dim: "<<dim<<"v: "<<data[i].data[0]<<", "<<v[i].data[1]<<std::endl;
  std::cout<<"***********"<<std::endl;
  #endif

  std::nth_element(data.begin(),data.begin()+(r)/2, data.begin()+r, [axis](const point<TYPE>& lhs, const point<TYPE>& rhs) {
      return lhs.data[axis] < rhs.data[axis];
    });

  
  int ind = allc.size();
  if(ind!=0) allc[ind-1].left = ind;
  allc.push_back(node<T>(data[tmp]));
  return tmp;
}
   
template<typename T>
int compute_median_s(std::vector<node<T>>& allc, std::vector<point<T>>& data, unsigned int nd, int axis,const int l,const int r){
  if(l>r){
    return -1;
  }
  else if(l==r){
    int ind = allc.size();
    allc.push_back(node<T>(data[l]));
    return ind;
  }
  //find median                                                                                                                                                                                                   
  int tmp=l+ (r-l)/2;
  std::nth_element(data.begin()+l,data.begin()+l+(r-l)/2, data.begin()+r, [axis](const point<TYPE>& lhs, const point<TYPE>& rhs) {
      return lhs.data[axis] < rhs.data[axis];
    });


  int ind = allc.size();
  allc.push_back(node<T>(data[tmp])); // filling all the data from the left

  axis= (axis+1)%nd; //                                                                                                                                                                     

    
  allc[ind].left = compute_median_s(allc,data,nd,axis,l, tmp-1);


  allc[ind].right = compute_median_s(allc,data,nd,axis,tmp+1, r);
  
  return ind;
}

/*******Print tree*******/
template<typename T>
void print2DUtil(int root, int space, const std::vector<node<T>>& allc){
  // Base case
  if (root == -1)
    return;
 
  // Increase distance between levels
  space += COUNT;
 
  // Process right child first
  print2DUtil(allc[root].right, space,allc);
 
  // Print current node after space
  // count
  std::cout<<std::endl;
  for (int i = COUNT; i < space; ++i) std::cout<<" ";
    
  std::cout<<allc[root].p.data[0]<<","<<allc[root].p.data[1]<<"\n";
 
  // Process left child
  print2DUtil(allc[root].left, space,allc);
}
 
// Wrapper over print2DUtil()
template<typename T>
void print2D(int root,const std::vector<node<T>>& allc){
  // Pass initial space count as 0
  print2DUtil(root, 0,allc);
}
 
/***********************************/

/********** Aux functions **********/
void position_in_the_tree(int& levels, int& level, int& axis, int& mydim,const int dim_tot,const int rank,const int size) noexcept {
  level=int(floor(log2(rank+1))); 
  axis=level%ndim;
  mydim=dim_tot;
  levels = int(log2(size));

}
template <typename T>
void dataset_creation(std::vector<point<T>>& mydata, const int mydim) noexcept {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,100.0);
  for(std::size_t i{0}; i<mydim; ++i){
    auto tmpx= TYPE(distribution(generator));
    auto tmpy = TYPE(distribution(generator));
    mydata[i]=point<TYPE>(tmpx,tmpy);
  }
}
template <typename T>
void merge_trees(std::vector<node<T>>& pt, std::vector<node<T>>& st, const int size_p){
  if(size_p!=0)  pt[size_p-1].left = size_p;
  for(std::size_t i{0}; i<st.size(); ++i){
    st[i].right = (st[i].right!=-1)*(st[i].right+size_p)+(st[i].right==-1)*(st[i].right);
    st[i].left = (st[i].left!=-1)*(st[i].left+size_p)+(st[i].left==-1)*(st[i].left);
    pt.push_back(std::move(st[i]));
  }

}

void check_size_proc(const int size){
  float pow=log2(size);
  if (pow!=int(pow))
    throw std::invalid_argument("Number of process must be a power of 2");
  else return;
}

void check_size_dataset(const int sizeproc, const int size_dataset){
  if (2*sizeproc > size_dataset)
    throw std::invalid_argument("2*number of process must be > of size of dataset");
  else return;
}
/*******************/

/******* Build tree with msg dimension *******/

void child_nodes(const int n, int kids[2]){
  int level= floor(log2(n)); // level in which the node is located
  int pos = n-pow(2,level); // position in the level
  kids[0]= pow(2,level+1)+2*pos; // sum of the nodes in previous levels +1 +shift due to the position
  kids[1]= pow(2,level+1)+1+2*pos;
}

void dimtree_fill(std::vector<std::pair<int,int>>& v,int *stride , int *boundaries ,int pos, int val,const int levels,const int  level, const int parent, int flag, std::vector<int>& list_rank){
  int stop=pow(2,levels+1)-1;
  int tag;
  tag = (flag!=-1)*(parent + pow(2,level))+ (flag==-1)*parent;
  if(pos>=stop){
    stride[parent]+=val;
    return;
  }
  if(parent!=tag){
    list_rank.push_back(tag);
    boundaries[tag]=val-1; 
  }
  int kid[2];
  child_nodes(pos+1,kid);
  --kid[0];
  --kid[1];

  ++stride[tag];//updating stride

  auto tmp = (val-1)/2;
  auto tmp1 = (val-1)/2+1*(((val-1))%2);
  
  v[pos]=std::pair<int,int>(tag,val);

  
  dimtree_fill(v,stride,boundaries, kid[0],tmp,levels,level+1,tag,-1, list_rank);
  dimtree_fill(v,stride,boundaries, kid[1], tmp1,levels,level+1,tag,0, list_rank);
} 

void dimtree_search(const std::vector<std::pair<int,int>> v, int flag, int& dim_r, int dim_s[2],int rank){
      while(flag!=-1){ // searching for dimension in dim_tree
        if(v[flag].first==rank){
          dim_r=v[flag].second;
          dim_s[0]=v[2*flag+1].second; // 2^(level+1)+2*(pos_in_the_level) = 2*flag (+1 beacause is 0-indexed)
          dim_s[1]=v[2*flag+2].second;
	        flag=-2;
        }
	      ++flag;

      }
}
