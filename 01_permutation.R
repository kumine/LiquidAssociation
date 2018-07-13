#=====================================================================
#Genome-wide LA significance assessment by simulation 

#==================================================
# gene expression data which had qqnormed. 368(inds) * 24907 (genes)

load("exp_qqnorm.Rdata")

#==========================================================
#consider all value as a pool. sample 368 values as  gene random expression
sample_data=function(){
  sapply(1:368,function(i){
    nr=sample(1:nrow(exp),1)
    nc=sample(1:ncol(exp),1)
    exp[nr,nc]
  })
}

perm_LA=function(...){
  X=sample_data()
  Y=sample_data()
  XY=X*Y
#consider all real gene expression as Z gene
# LA=mean(x*y*z)=sum(x*y*z)/length(x)
# we can use matrix cross-product to accelerate. 
  LA=crossprod(XY,exp)[1,]/length(XY)
  range(LA)
}

#=======================

N=1000 # the simulation times. here just use 1000. In fact, we use 1 million.

ref_LA=sapply(1:N,perm_LA)


#=================================
#the random genome-wide LA

GW_LA=function(...){
  X=exp[,sample(1:ncol(exp),1)]
  Y=exp[,sample(1:ncol(exp),1)]
  XY=X*Y
  LA=crossprod(XY,exp)[1,]/length(XY)
  range(LA)
}


g_LA=sapply(1:N,GW_LA)

#=========================================================
#Q-Q plot

plot(sort(ref_LA[1,]),sort(g_LA[1,]),
     ylab="Genome-wide LA score",
     xlab="Random LA scores",
     main="negative"
     )

plot(sort(ref_LA[2,]),sort(g_LA[2,]),
     ylab="Genome-wide LA score",
     xlab="Random LA scores",
     main="negative"
)

