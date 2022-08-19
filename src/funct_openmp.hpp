#include <iostream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <memory.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <chrono>
#define DEBUG2 1
#define COUNT 10
#define xdim 10
#define ydim 10
#define ndim 2
using  uint = unsigned int;

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
                                                                                                                                                                                                              
  void print_node() noexcept {                                                                                                                                                                               
    std::cout<<"The content of the node is: "<<std::endl;                                                                                                                                                  
    for(std::size_t i{0}; i<2; ++i)                                                                                                                                                                     
      std::cout<<i<<"th element is = "<<data[i]<<std::endl;                                                                                                                                                  
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
int compute_median_p(std::vector<node<T>>& allc, std::vector<point<T>>& data, unsigned int nd, int axis,const int l,const int r){
  if(r-l<20){
    return compute_median_s(allc,  data,  nd,  axis,l,r);
  }
  else{

    if(l>r){
      return -1;
    }
    else if(l==r){
      int ind = allc.size();
      allc.push_back(node<T>(data[l]));
      return ind;
    }
    //find median                                                                                                                                                                                                    
  int tmp=l+(r-l)/2;
  std::nth_element(data.begin()+l,data.begin()+l+(r-l)/2, data.begin()+r, [axis](const point<TYPE>& lhs, const point<TYPE>& rhs) {
        return lhs.data[axis] < rhs.data[axis];
   });
    int ind = allc.size();
    allc.push_back(std::move(node<T>(data[tmp]))); // filling all the data from the left

    axis= (axis+1)%nd;                                                                                                                                                                      

#pragma omp task shared(allc,data) firstprivate(l,axis,tmp,ind) untied
    allc[ind].left = compute_median_p(allc,data,nd,axis,l, tmp-1);

#pragma omp task shared(allc,data) firstprivate(r,axis,tmp,ind) untied
    allc[ind].right = compute_median_p(allc,data,nd,axis,tmp+1, r);
    
    return ind;
  }
  
}
   
template<typename T>
int compute_median_s(std::vector<node<T>>& allc, std::vector<point<T>>& data, unsigned int nd, int axis,const int l,const int r){
  if(l>r){
    return -1;
  }
  else if(l==r){
    int ind = allc.size();
    allc.push_back(std::move(node<T>(data[l])));
    return ind;
  }
  //find median                                                                                                                                                                
  int tmp=l+(r-l)/2;
  std::nth_element(data.begin()+l,data.begin()+l+(r-l)/2, data.begin()+r, [axis](const point<TYPE>& lhs, const point<TYPE>& rhs) {
        return lhs.data[axis] < rhs.data[axis];
   });

  int ind = allc.size();
  allc.push_back(std::move(node<T>(data[tmp]))); // filling all the data from the left

  axis= (axis+1)%nd;                                                                                                                                                                      

    
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
