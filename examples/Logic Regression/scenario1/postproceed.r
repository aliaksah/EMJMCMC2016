library(hash)
library(stringi)
#setwd("/mn/sarpanitu/ansatte-u2/aliaksah/abeldata/logic-g-prior/logscenario1")

#experiment i
temp = list.files(pattern="post3etaG_*")
myfiles = lapply(FUN = read.csv,X = temp)

X<- as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
rhash<-hash()

scenairo=6.5
if(scenairo == 1|2)
{
  tps<-array(data = 0, dim =c(1000,3))
  tps[,1]<-((1-X$V1)*X$V4)
  tps[,2]<-(X$V8*X$V11)
  tps[,3]<-(X$V5*X$V9)
}

if(scenairo==3)
{
  tps<-array(data = 0, dim =c(1000,3))
  tps[,1]<- X$V4*X$V17*X$V30*X$V10
  tps[,2]<- X$V7*X$V20*X$V12
  tps[,3]<- X$V2*X$V9
}

if(scenairo == 4)
{
  tps<-array(data = 0, dim =c(1000,3))
  tps[,1]<-(X$V1*X$V4)
  tps[,2]<-(X$V8*X$V11)
  tps[,3]<-(X$V5*X$V9)
}

if(scenairo == 5)
{
  tps<-array(data = 0, dim =c(1000,4))
  tps[,1]<- X$V4*X$V17*X$V30*X$V10
  tps[,2]<- X$V7*X$V20*X$V12
  tps[,3]<-(X$V2*X$V9)
  tps[,4]<-(X$V37)
}

if(scenairo==6)
{
  tps<-array(data = 0, dim =c(1000,8))
  tps[,1]<- X$V4*X$V17*X$V30*X$V10
  tps[,2]<- X$V37*X$V20*X$V12
  tps[,3]<- X$V2*X$V9
  tps[,4]<- X$V1*X$V27*X$V3
  tps[,5]<- as.integer((X$V50*X$V19+X$V13*X$V11)>0)
  tps[,6]<- X$V21*X$V18
  tps[,7]<- X$V7
  tps[,8]<- X$V8
}

if(scenairo==6.5)
{
  tps<-array(data = 0, dim =c(1000,11))
  tps[,1]<- X$V4*X$V17*X$V30*X$V10
  tps[,2]<- X$V37*X$V20*X$V12
  tps[,3]<- X$V2*X$V9
  tps[,4]<- X$V1*X$V27*X$V3
  tps[,5]<- as.integer((X$V50*X$V19+X$V13*X$V11)>0)
  tps[,6]<- X$V21*X$V18
  tps[,7]<- X$V7
  tps[,8]<- X$V8
  tps[,9]<- X$V50*X$V19*X$V13*X$V11
  tps[,10]<- X$V50*X$V19
  tps[,11]<- X$V13*X$V11
}

N<-length(myfiles)
alpha<-0.5
clear(rhash)
fdr<-0
for(i in 1:min(100,N))
{
  tpi<-0
  fpi<-0
  for(j in 1:length(myfiles[[i]]$posterior))
  {
    if(myfiles[[i]]$posterior[j]>=alpha)
    {
      expr<-as.character(myfiles[[i]]$tree[j])
      print(expr)
      res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
      res[,1]<-res[,1]-res[,2]
      TP<-0
      for(u in 1:dim(tps)[2])
      {
        if(sum(tps[,u]-res[,2])==0)
        {
          TP<-1
          break
        }
      }
      if(TP)
      {
        tpi<-tpi+1
      }else{
        fpi<-fpi+1
      }
      
      ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),1,expr)
      if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
        rhash[[ress[1]]]<-ress
      else
      {
        if(ress[1] %in% keys(rhash))
        {
          rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3])) + as.numeric(1)
          if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
            rhash[[ress[1]]][4]<-expr
        }
        else
        {
          rhash[[ress[2]]][3]<- (as.numeric(rhash[[ress[2]]][3])) + as.numeric(1)
          if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
            rhash[[ress[2]]][4]<-expr
        }
      }
    }
  }
  print(fpi)
  print(tpi)
  if((fpi+tpi)>0)
  fdr<-fdr+(fpi/(fpi+tpi))
}
print("FDRFDRFDRFDR!!!!:")
print(fdr/100)



res<-as.data.frame(t(values(rhash)[c(4,3),]))
res$V1<-as.numeric(as.character(res$V1))
res<-res[order(res$V1, decreasing = T),]
colnames(res)<-c("posterior","tree")
write.csv(x = t(values(rhash))[,c(3,4)],file = "expn3.csv",row.names = F,col.names = F)



simplifyposteriors<-function(X,posteriors,th=0.0001,thf=0.5)
{
  posteriors<-posteriors[-which(posteriors[,2]<th),]
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    res[,1]<-res[,1]-res[,2]
    ress<-c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
      rhash[[ress[1]]]<-ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]<-expr
      }
      else
      {
        rhash[[ress[2]]][3]<- (as.numeric(rhash[[ress[2]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
          rhash[[ress[2]]][4]<-expr
      }
    }
    
  }
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1<-as.numeric(as.character(res$V1))
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  return(res)
}



