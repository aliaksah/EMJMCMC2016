#inference
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


#***********************IMPORTANT******************************************************
# if a multithreaded backend openBLAS for matrix multiplications
# is installed on your machine, please force it to use 1 thread explicitly
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)
#***********************IMPORTANT******************************************************


#read in the train and test data sets
test = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/breast%20cancer/test.csv",header = T,sep=",")[,-1]
train = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/breast%20cancer/train.csv",header = T,sep=",")[,-1]

#transform the train data set to a data.example data.frame that EMJMCMC class will internally use
data.example = as.data.frame(train)

#perfrom garbage collection
gc()
pinferunemjmcmc = function(s.seed = 0, n.cores = 4, mcgmj = mcgmjpar, report.level =  0.5, simplify = F, num.mod.best = 1000, predict = F,  test.data = 1, link.function = function(z)z, runemjmcmc.params)
{
  
  if(predict)
  {
    runemjmcmc.params$save.beta = T
    
    if(length(test.data)==0)
    {
      warning("Test data is not provided. No predictions will be made!")
    }
  }
  
  params = list(runemjmcmc.params)[rep(1,n.cores)]
  for(i in 1:n.cores)
  {
    params[[i]]$test = test.data
    params[[i]]$link = link.function
    params[[i]]$predict = predict
    params[[i]]$NM = num.mod.best
    params[[i]]$cpu=i + s.seed*n.cores
  }
  M = n.cores
  #results = runpar.infer(params[[1]])
  results=mcgmj(X = params,FUN = runpar.infer,mc.cores = n.cores)
  #print(results)
  #clean up
  gc()
  #prepare the data structures for final analysis of the runs
  compmax = runemjmcmc.params$interact.param$Nvars.max + 1
  resa=array(data = 0,dim = c(compmax,M*3))
  post.popul = array(0,M)
  max.popul = array(0,M)
  nulls=NULL
  not.null=1
  #check which threads had non-zero exit status
  for(k in 1:M)
  {
    if(length(results[[k]])<=1||length(results[[k]]$cterm)==0||length(results[[k]]$p.post)!=runemjmcmc.params$interact.param$Nvars.max)
    {
      nulls=c(nulls,k)
      #warning(paste0("Thread ",k,"did not converge or was killed by OS!"))
      next
    }
    else
    {
      not.null = k
    }
    
  }
  
  if(length(nulls) == M)
  {
    warning("All threads did not converge or gave an error! Returning stats from the threads only!")
    return(list(feat.stat = NULL,predictions = NULL,allposteriors = NULL, threads.stats = results))
  }
  
  
  #for all of the successful runs collect the results into the corresponding data structures
  for(k in 1:M)
  {
    if(k %in% nulls)
    {
      results[[k]]=results[[not.null]]
    }
    max.popul[k]=results[[k]]$cterm
    post.popul[k]=results[[k]]$post.populi
    resa[,k*3-2]=c(results[[k]]$fparam,"Post.Gen.Max")
    resa[,k*3-1]=c(results[[k]]$p.post,results[[k]]$cterm)
    resa[,k*3]=rep(post.popul[k],length(results[[k]]$p.post)+1)
    
  }
  #renormalize estimates of the marginal inclusion probabilities
  #based on all of the runs
  ml.max=max(max.popul)
  post.popul=post.popul*exp(-ml.max+max.popul)
  p.gen.post=post.popul/sum(post.popul)
  
  #perform BMA of the redictions across the runs
  pred = NULL
  if(predict){
    pred = results[[1]]$preds*p.gen.post[1]
    if(M > 1) {
      
      for(i in 2:M)
      {
        
        pred=pred+results[[i]]$preds*p.gen.post[i]
        
      }
    }
  }
  hfinal=hash()
  for(ii in 1:M)
  {
    resa[,ii*3]=p.gen.post[ii]*as.numeric(resa[,ii*3-1])
    resa[length(resa[,ii*3]),ii*3]=p.gen.post[ii]
    if(p.gen.post[ii]>0)
    {
      for(jj in 1:(length(resa[,ii*3])-1))
      {
        if(resa[jj,ii*3]>0)
        {
          if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
            hfinal[[resa[jj,ii*3-2]]]=as.numeric(resa[jj,ii*3])
          else
            hfinal[[resa[jj,ii*3-2]]]=hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
        }
        
      }
    }
  }
  
  posteriors=values(hfinal)
  clear(hfinal)
  #delete the unused further variables
  rm(hfinal)
  rm(resa)
  rm(post.popul)
  rm(max.popul)
  #simplify the found trees and their posteriors
  posteriors=as.data.frame(posteriors)
  posteriors=data.frame(X=row.names(posteriors),x=posteriors$posteriors)
  posteriors$X=as.character(posteriors$X)
  res1 = NULL
  if(simplify){
    
    
    res1=simplifyposteriors.infer(X = runemjmcmc.params$data,posteriors = posteriors, thf = report.level,resp = as.character(runemjmcmc.params$formula)[2])
    rownames(res1) = NULL
    res1$feature = as.character(res1$feature)
  }
  posteriors=posteriors[order(posteriors$x, decreasing = T),]
  colnames(posteriors)=c("feature","posterior")
  rm(params)
  gc()
  return(list(feat.stat = cbind(res1$feature,res1$posterior),predictions = pred,allposteriors = posteriors, threads.stats = results))
  
}




#specify the initial formula
formula1 = as.formula(paste(colnames(test)[31],"~ 1 +",paste0(colnames(test)[-31],collapse = "+")))

#a set of nonlinearities that will be used in the DBRM model
p0 = function(x) log(abs(x)+0.00001) #2
pm1 = function(x) (x+0.00001)^(-1)
pm2 = function(x) (x+0.00001)^(-2)
pm05 = function(x) (abs(x)+0.00001)^(-0.5)
p05 = function(x) (abs(x)+0.00001)^(0.5)
p2 = function(x) x^(2)
p3 = function(x) x^(3)

p0p0 = function(x) p0(x)*p0(x) #4
p0pm1 = function(x) p0(x)*(x+0.00001)^(-1)
p0pm2 = function(x) p0(x)*(x+0.00001)^(-2)
p0pm05 = function(x)p0(x)*(abs(x)+0.00001)^(-0.5)
p0p05 = function(x) p0(x)*(abs(x)+0.00001)^(0.5)
p0p1 = function(x) p0(x)*x
p0p2 = function(x) p0(x)*x^(2)
p0p3 = function(x) p0(x)*x^(3)
#relu=function(x)max(0,x)

#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients
#prediction
estimate.bas.glm.cpen = function(formula, data, link, distribution, family, prior = "BIC", logn,r = 0.1,yid=31,relat1=c("p0","p2","p3","p05","pm05","pm1","pm2"),relat2=c("p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2"))
{
  capture.output({out = glm(family = family,formula = formula,data = data)})
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
  
  
  sj=(stri_count_fixed(str = fparam, pattern = "*"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "+"))
  for(rel in relat1)
    sj=sj+log(2)*(stri_count_fixed(str = fparam, pattern = rel))
  for(rel in relat2)
    sj=sj+log(4)*(stri_count_fixed(str = fparam, pattern = rel))
  
  mlik = ((-out$deviance +2*log(r)*sum(sj)))/2
  return(list(mlik = mlik,waic = -(out$deviance + 2*out$rank) , dic =  -(out$deviance + logn*out$rank),summary.fixed =list(mean = coefficients(out))))
}


#define the number or cpus
M = 32
#define the size of the simulated samples
NM= 1000
#define \k_{max} + 1 from the paper
compmax = 46
#define treshold for preinclusion of the tree into the analysis
th=(10)^(-5)
#define a final treshold on the posterior marginal probability for reporting a tree
thf=0.05
#specify tuning parameters of the algorithm for exploring DBRM of interest
#notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
#
#specify the link function that will be used in the prediction phase
g=function(x)
{
  return((x = 1/(1+exp(-x))))
}
results=array(0,dim = c(2,100,5))

for(j in 1:100)
{
  
  #specify the initial formula
  set.seed(j)
  
  res1 = pinferunemjmcmc(s.seed = j, n.cores = M, report.level =  0.2 , num.mod.best = NM,simplify = T, predict = T,
                         test.data = as.data.frame(test),link.function = g,
                         runemjmcmc.params =list(formula = formula1,data = data.example,
                                                 estimator = estimate.bas.glm.cpen,
                                                 estimator.args =  list(data = data.example, yid=1, family = binomial(), logn = log(dim(data.example)[1]),r=1/(dim(data.example)[1])),
                                                 recalc_margin = 249, save.beta = T,interact = T,gen.prob = c(5,3,1,0,0),relations=c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2"),
                                                 relations.prob = rep(0.1,16),interact.param=list(allow_offsprings=3,mutation_rate = 100,
                                                                                                  last.mutation=10000, max.tree.size = 4, Nvars.max = 45,p.allow.replace=0.5,p.allow.tree=0.3,p.nor=0.3,p.and = 0.9),
                                                 n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,
                                                 burn.in = 100,print.freq = 1000,advanced.param = list(
                                                   max.N.glob=as.integer(10),
                                                   min.N.glob=as.integer(5),
                                                   max.N=as.integer(3),
                                                   min.N=as.integer(1),
                                                   printable = F)))
  
  print(res1$feat.stat)
  out = as.integer(res1$threads.stats[[1]]$preds>=0.5)
  
  ps=which(test$X==1)


  ns=which(test$X==0)
  
  results[1,j,1]=  (1-sum(abs(out-test$X[1:length(out)]))/length(out))
  results[1,j,2]=   sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))
  results[1,j,3] =   sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))
  
  out = as.integer(res1$predictions>=0.5)
  
  
  results[2,j,1]= (1-sum(abs(out-test$X[1:length(out)]))/length(out))
  results[2,j,2]=   sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))
  results[2,j,3] =  sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))
  
  write.csv(x =res1$feat.stat,row.names = F,file = paste0("posteriorbreastcancerFPBIC_",j,"_interactions.csv"))
  print(paste0("end simulation ",j))
  
  #print the run's metrics and clean the results
  write.csv(file =paste0("resultsrunbreastcancerBIC_",j,"_interactions.csv"),x= results[,j,])
  #rm(results)
  
  
  print(sqrt(mean((res1$predictions - test$Age)^2)))
  
  gc()
}


for(j in 1:100)
{
  
  tmp = read.csv(paste0("resultsrunbreastcancerBIC_",j,"_interactions.csv"))
  
  results[1,j,1]=  tmp[1,2]
  results[1,j,2]=  tmp[1,3]
  results[1,j,3] =   tmp[1,4]
  
  if(tmp[2,2] == 0)
    print(j)
  
  results[2,j,1]=  tmp[2,2]
  results[2,j,2]=   tmp[2,3]
  results[2,j,3] =   tmp[2,4]
  
}

#make the joint summary of the runs, including min, max and medians of the performance metrics
summary.results=array(data = NA,dim = c(2,15))

for(i in 1:2){
  for(j in 1:5)
  {
    summary.results[i,(j-1)*3+1]=min(results[i,,j])
    summary.results[i,(j-1)*3+2]=median(results[i,,j])
    summary.results[i,(j-1)*3+3]=max(results[i,,j])
  }
}
summary.results=as.data.frame(summary.results)

View(summary.results[,c(2,1,3,5,4,6,8,7,9)])

#featgmj = hash()

simplifyposteriors<-function(X,post,th=0.0001,thf=0.1,y = "X")
{
  posteriors = (cbind(as.character(post[,1]),as.numeric(as.character(post[,2]))))
  
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    
    
    expr<-posteriors[i,1]
    print(expr)
    res<-model.matrix(data=X,object = as.formula(paste0(y,"~",expr)))
    ress<-c(stri_flatten(round(sum(res[,2]),digits = 4),collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!((ress[1] %in% keys(rhash))))
      rhash[[ress[1]]]<-ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {new
        rhash[[ress[1]]][3]<- (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]<-expr
      }
    }
    
  }
  res<-as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1<-as.numeric(as.character(res$V1))
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  #res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  
  row.names(res) = 1:length(res$posterior)
  return(res)
}




featgmj = hash()

for(j in 1:100)
{
  print(j)
  tmpf = read.csv(paste0("posteriorbreastcancerFPBIC_",j,"_interactions.csv"))
  #tmp = simplifyposteriors(X = data.example,post =tmpf,y = "Age")
  for(feat in as.character(tmpf$V1))
  {
    if(!has.key(hash = featgmj,key =  feat ))
    {
      featgmj[[feat]] = as.numeric(1)
    } else{
      
      featgmj[[feat]] =as.numeric(featgmj[[feat]]) + 1
    }
  }
}

tmp = simplifyposteriors(th=0.0001,thf=0.01, X = data.example,post =as.data.frame(cbind(keys(featgmj),as.numeric(values(featgmj)))),y = "X")

write.csv(x =cbind(tmp[1:15,c(2,1)],tmp[16:30,c(2,1)]),row.names = F,file = "breastcancerfeatFP1_interactions.csv",sep = "&")
#print(paste0("end simulation ",j))


