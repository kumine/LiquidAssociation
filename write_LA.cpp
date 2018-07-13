#include<Rcpp.h> 
#include<vector>     
#include<string>     
#include<fstream>
#include <sstream>


using namespace std; 
using namespace Rcpp;

//output LA result matrix
// res:  LA result.
// x:   x  gene name
//colid: Column name
//rowid: row name


// [[Rcpp::export]]
int write_LA(NumericMatrix res,string x,double cut,CharacterVector colid,CharacterVector rowid,std::string file){
  
  ofstream myfile;
  myfile.open (file.c_str(),ios::app);
  
  for(int i=0;i<res.nrow();i++){
    for(int j=i+1; j<res.ncol();j++){
      if(abs(res(i,j))>=cut){
        myfile <<x<<"\t"<<rowid(i)<<"\t"<<colid(j)<<"\t"<<res(i,j)<<endl;
      }
    }
  }
  
  return 0;
}


// x*Y

// [[Rcpp::export]]
NumericMatrix DXY(NumericVector X,NumericMatrix Y ){
  NumericMatrix XY(Y.nrow(),Y.ncol());
  for(int i=0;i<Y.nrow();i++){
    for(int j=0; j<Y.ncol();j++){
       XY(i,j)=Y(i,j)*X(i);
    }
  }
  return XY;
  
}



