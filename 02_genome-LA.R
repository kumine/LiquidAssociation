library("Rcpp")
sourceCpp("write_LA.cpp")


#==================================================
# gene expression data which had qqnormed. 368(inds) * 24907 (genes)
load("exp_qqnorm.Rdata")

cutoff=0.33                    # LA Threshold.
ofile="genome_LA_result.txt"   #result output file


#data=exp           # realitily we used all data
data=exp[,1:100]   # smaller dataset used to test.

gnm=colnames(data)

start=1
end=ncol(data)-2       

for(i in start:end){
  cat("start run:",i,"\n")
  xnm=gnm[i]
  x=data[,i]
  #XY=apply(data[,(i+1):ncol(data)],2,function(y){x*y})
  XY=DXY(x,data[,(i+1):ncol(data),drop=F])    # much faster than apply.
  colnames(XY)=colnames(data)[(i+1):ncol(data)]
  LA=crossprod(XY,data[,(i+1):ncol(data),drop=F])/nrow(data)
  write_LA(LA,xnm,cut=cutoff,colid=colnames(LA),rowid=rownames(LA),file=ofile)
}


#================================================
# verification
i=4
xnm=gnm[i]
x=data[,i]

XY=DXY(x,data[,(i+1):ncol(data),drop=F])    
colnames(XY)=colnames(data)[(i+1):ncol(data)]
LA=crossprod(XY,data[,(i+1):ncol(data),drop=F])/nrow(data)

yi=8
zi=10

ynm=gnm[yi]
znm=gnm[zi]

LA[ynm,znm]

mean(x*exp[,ynm]*exp[,znm])

