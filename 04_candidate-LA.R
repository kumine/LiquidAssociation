library("Rcpp")
library(data.table)
sourceCpp("write_LA.cpp")


#==================================================
# gene expression data which had qqnormed. 368(inds) * 24907 (genes)
load("exp_qqnorm.Rdata")

cutoff=0.33                    # LA Threshold.
ofile="candidate_LA_result.txt"   #result output file

infile="candidate_gene_list.txt"  #candiate gene list

#======================================================

xlist=read.table(infile,head=F,as.is=T)[,1]

data=exp

X=data[,xlist]
Y=X
gnm=colnames(X)

start=1
end=ncol(X)

for(i in start:end){
  cat("start run:",i,"\n")
  xnm=gnm[i]
  x=X[,i]
  XY=DXY(x,Y)    # much faster than apply.
  colnames(XY)=colnames(Y)
  LA=crossprod(XY,data)/nrow(data)
  write_LA(LA,xnm,cut=cutoff,colid=colnames(LA),rowid=rownames(LA),file=ofile)
}

#==========================================================
#remove any x==y, x==z, y==z

d=fread(ofile,head=F,data.table = F)
index=sapply(1:nrow(d),function(i){
  !any(d[i,1]==d[i,2],d[i,1]==d[i,3],d[i,2]==d[i,3])
})
d=d[index,]
write.table(d,file=ofile,quote=F,col.names = F,row.names = F,sep="\t")

#================================================
#just  verification  

LA=sapply(1:10, function(i){
  x=data[,d[i,1]]
  y=data[,d[i,2]]
  z=data[,d[i,3]]
  mean(x*y*z)
})

LA





