#===================================================================================
# Definition of LA-scouting genes in genome-wide significant LA triplets by MLA
library(data.table)

load("exp_qqnorm.Rdata")
infile="genome_LA_result.txt"
ofile="genome_LA_result_MLA.txt"


chuck=6   # normally the LA_result is very big, we can read it by chuck. chuck=6 is just for test
start=1


br=c(-3,-0.429,0.429,3) 
index=1:nrow(exp)
expPos=apply(exp,2,function(z){
  split(index,cut(z, breaks = br, labels = c("lower", "middle", "high")))
})

MLA=function(index){
  x=exp[,index[1]]
  y=exp[,index[2]]
  pos=expPos[[index[3]]] 
  (-1.072588*cor(x[pos[[1]]],y[pos[[1]]])+1.072588*cor(x[pos[[3]]],y[pos[[3]]]))/3
}

cat("gene1\tgene2\tgene3\tMLA(Z=gene1)\tMLA(Z=gene2)\tMLA(Z=gene3)\tLA\n",file=ofile,sep="")
while(T){
  data=fread(infile,sep="\t",head=F,skip=start-1,nrow=chuck,colClasses = rep("numeric",4),data.table=F)
  data=as.matrix(data)
  res=apply(data,1,function(x){
    mla=c(
      MLA(x[c(3,2,1)]), #X
      MLA(x[c(1,3,2)]), #Y
      MLA(x[1:3])            #Z
    )
    index=order(abs(mla))
    mla=as.character(round(mla[index],3))
    c(x[index],mla,x[4])
  })
  res=t(res)
  write.table(res,file=ofile,append=T,col.names=F,row.names = F,sep="\t",quote = F)
  if(nrow(data)<chuck){break}
  start=start+chuck
}














