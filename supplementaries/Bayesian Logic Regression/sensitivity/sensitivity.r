#read the appropriate version of the package
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package2.r")
library(bindata)
#define the function for perfroming parallel computations
parall.gmj = mclapply
#define some function for cleaning up after parallel computations
library(inline)
includes = '#include <sys/wait.h>'
code = 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait = cfunction(body=code, includes=includes, convention='.C')

#define the function estimating parameters of a given Gaussian logic regression with robust g prior
estimate.logic.lm.tCCH = function(formula = NULL, data, n=1000, m=50, r = 1, p.a = 1, p.b = 2, p.r = 1.5, p.s = 0, p.v=-1, p.k = 1,k.max=21)
{
  if(is.na(formula)||is.null(formula))
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  if(fmla.proc[2]=="-1")
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  out = lm(formula = formula,data = data)
  p = out$rank
  if(p>k.max)
  {
    return(list(mlik = -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = 0)))
  }
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj=(stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  p.v = (n+1)/(p+1)
  R.2 = summary(out)$r.squared
  
  mlik = (-0.5*p*log(p.v) -0.5*(n-1)*log(1-(1-1/p.v)*R.2) + log(beta((p.a+p)/2,p.b/2)) - log(beta(p.a/2,p.b/2)) + log(phi1(p.b/2,(n-1)/2,(p.a+p.b+p)/2,p.s/2/p.v,R.2/(p.v-(p.v-1)*R.2))) - hypergeometric1F1(p.b/2,(p.a+p.b)/2,p.s/2/p.v,log = T)+log(Jprior) + p*log(r)+n)
  if(mlik==-Inf||is.na(mlik)||is.nan(mlik))
    mlik = -10000
  return(list(mlik = mlik,waic = AIC(out)-n , dic =  BIC(out)-n,summary.fixed =list(mean = coef(out))))
}


#define the function estimating parameters of a given Gaussian logic regression with Jeffrey's prior
estimate.logic.lm = function(formula= NULL, data, n, m, r = 1,k.max=21)
{
  if(is.na(formula)||is.null(formula))
    return(list(mlik =  -10000,waic =10000 , dic =  10000,summary.fixed =list(mean = 1)))
  out = lm(formula = formula,data = data)
  p = out$rank
  if(p>k.max)
  {
    return(list(mlik = -10000,waic = 10000 , dic =  10000,summary.fixed =list(mean = 0)))
  }
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]]
  sj=(stri_count_fixed(str = fparam, pattern = "&"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "|"))
  sj=sj+1
  Jprior = prod(factorial(sj)/((m^sj)*2^(2*sj-2)))
  mlik = (-BIC(out)+2*log(Jprior) + 2*p*log(r)+n)/2
  if(mlik==-Inf)
    mlik = -10000
  return(list(mlik = mlik,waic = AIC(out)-n , dic =  BIC(out)-n,summary.fixed =list(mean = coef(out))))
}

#define the function simplifying logical expressions at the end of the search
simplifyposteriors=function(X,posteriors,th=0.0001,thf=0.5)
{
  posteriors=posteriors[-which(posteriors[,2]<th),]
  rhash=hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr=posteriors[i,1]
    print(expr)
    res=model.matrix(data=X,object = as.formula(paste0("V1~",expr)))
    res[,1]=res[,1]-res[,2]
    ress=c(stri_flatten(res[,1],collapse = ""),stri_flatten(res[,2],collapse = ""),posteriors[i,2],expr)
    if(!(ress[1] %in% values(rhash)||(ress[2] %in% values(rhash))))
      rhash[[ress[1]]]=ress
    else
    {
      if(ress[1] %in% keys(rhash))
      {
        rhash[[ress[1]]][3]= (as.numeric(rhash[[ress[1]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[1]]][4])>stri_length(expr))
          rhash[[ress[1]]][4]=expr
      }
      else
      {
        rhash[[ress[2]]][3]= (as.numeric(rhash[[ress[2]]][3]) + as.numeric(ress[3]))
        if(stri_length(rhash[[ress[2]]][4])>stri_length(expr))
          rhash[[ress[2]]][4]=expr
      }
    }
    
  }
  res=as.data.frame(t(values(rhash)[c(3,4),]))
  res$V1=as.numeric(as.character(res$V1))
  res=res[which(res$V1>thf),]
  res=res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]=1
  colnames(res)=c("posterior","tree")
  return(res)
}



#define number of simulations
MM = 800
#define number of threads to be used
M = 32
#define the size of the simulated samples
NM= 1000
#define \k_{max} + 1 from the paper
compmax = 21
#define treshold for preinclusion of the tree into the analysis
th=(10)^(-5)
#define a final treshold on the posterior marginal probability for reporting a tree  
thf=0.05

#define a function performing the map step for a given thread
runpar=function(vect)
{
  
  set.seed(as.integer(vect[22]))
  do.call(runemjmcmc, vect[1:21])
  vals=values(hashStat)
  fparam=mySearch$fparam
  cterm=max(vals[1,],na.rm = T)
  ppp=mySearch$post_proceed_results_hash(hashStat = hashStat)
  post.populi=sum(exp(values(hashStat)[1,][1:NM]-cterm),na.rm = T)
  clear(hashStat)
  rm(hashStat)
  rm(vals)
  return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam))
}

#define type of the sensitivity analysis
type1=F
type2=F
type3=F
type4=F
#define whether to address Jeffrey's prior (T) or robust G prior (F)
typeJ=F

#start the sensitivity analysis
for(j in MM:1)
{
  print(paste0("begin iteration ",j))
  tryCatch({
    
    j.id = j%%80
    if(j.id == 0)
      j.id = 80
    jj=j
    compmax = 21
    if(j.id%in%1:10)
    {
      #vary slopes in the true signal, run with Jeffrey's prior
      type1=T
      type2=F
      type3=F
      type4=F
      typeJ=T
      jj=j.id
      id.sym = "T1J"
    }else if(j.id%in%11:20)
    {
      #vary slopes in the true signal, run with robust g prior
      type1=T
      type2=F
      type3=F
      type4=F
      typeJ=F
      jj=j.id-10
      id.sym = "T1G"
    }else if(j.id%in%21:30)
    {
      #vary sample size, run with with Jeffrey's prior
      type1=F
      type2=T
      type3=F
      type4=F
      typeJ=T
      jj=j.id-20
      id.sym = "T2J"
    }else if(j.id%in%31:40)
    {
      #vary sample size, run with with robust g prior
      type1=F
      type2=T
      type3=F
      type4=F
      typeJ=F
      jj=j.id-30
      id.sym = "T2G"
    }else if(j.id%in%41:50)
    {
      #vary population size, run with with Jeffrey's prior
      type1=F
      type2=F
      type3=T
      type4=F
      typeJ=T
      jj=j.id-40
      compmax = jj*15+1
      id.sym = "T3J"
    }else if(j.id%in%51:60)
    {
      #vary population size, run with with robust g prior
      type1=F
      type2=F
      type3=T
      type4=F
      typeJ=F
      jj=j.id-50
      compmax = jj*15+1
      id.sym = "T3G"
    }else if(j.id%in%61:70)
    {
      #vary misspecification correlation, run with with Jeffrey's prior
      type1=F
      type2=F
      type3=F
      type4=T
      typeJ=T
      jj=j.id-60
      id.sym = "T4J"
    }else if(j.id%in%c(71:80,0))
    {
      #vary misspecification correlation, run with with robust g prior
      type1=F
      type2=F
      type3=F
      type4=T
      typeJ=F
      jj=j.id-70
      id.sym = "T4G"
    }
    
    #prepare the data for a current simulation
    set.seed(jj)
    X2= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
    Y2=rnorm(n = 1000,mean = 1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*X2$V37,sd = 1)
    X2$Y2=Y2
    
    #specify the initial formula
    formula1 = as.formula(paste(colnames(X2)[51],"~ 1 +",paste0(colnames(X2)[-c(51)],collapse = "+")))
    data.example = as.data.frame(X2)
    
    #specify default tuning parameters of the algorithm for exploring the space of Bayesian logic regressions of interest
    #notice that allow_offsprings=1 corresponds to the GMJMCMC algorithm for Bayesian logic regression
    vect=list(formula = formula1,data = X2,estimator = estimate.logic.lm.tCCH,estimator.args =  list(data = data.example,n = 1000, m = 50),recalc_margin = 250, save.beta = F,interact = T,relations = c("","lgx2","cos","sigmoid","tanh","atan","erf"),relations.prob =c(0.4,0.0,0.0,0.0,0.0,0.0,0.0),interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =(compmax-1),p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9),n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,outgraphs=F,print.freq = 1000,advanced.param = list(
      max.N.glob=as.integer(10),
      min.N.glob=as.integer(5),
      max.N=as.integer(3),
      min.N=as.integer(1),
      printable = F))
    
    params = list(vect)[rep(1,M)]
    
    #adjust simulation and tuning paramters w.r.t. the type of sensitivity analysis and the model prior
    for(i in 1:M)
    {
      if(type1)
      {
        Y2=rnorm(n = 1000,mean = 1+jj/10*((7*X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*X2$V37),sd = 1)
        X2$Y2=Y2
        data.example = as.data.frame(X2)
        params[[i]]$data = data.example
        params[[i]]$estimator.args =  list(data = data.example,n = 1000, m = 50)
        params[[i]]$formula = formula1
        if(typeJ)
        {
          params[[i]]$estimator = estimate.logic.lm
        }
      }else if(type2)
      {
        NN = jj*100
        X2= as.data.frame(array(data = rbinom(n = 50*NN,size = 1,prob = runif(n = 50*NN,0,1)),dim = c(NN,50)))
        Y2=rnorm(n = NN,mean = 1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*X2$V37,sd = 1)
        X2$Y2=Y2
        data.example = as.data.frame(X2)
        params[[i]]$data = data.example
        params[[i]]$estimator.args =  list(data = data.example,n = NN, m = 50)
        params[[i]]$formula = formula1
        if(typeJ)
        {
          params[[i]]$estimator = estimate.logic.lm
        }
      }else if(type3)
      {
    
        data.example = as.data.frame(X2)
        params[[i]]$interact.param=list(allow_offsprings=1,mutation_rate = 300,last.mutation = 5000, max.tree.size = 4, Nvars.max =jj*15,p.allow.replace=0.9,p.allow.tree=0.2,p.nor=0,p.and = 0.9)
        if(typeJ)
        {
          params[[i]]$estimator = estimate.logic.lm
        }
        
      }else if(type4)
      {
        
        ## Construct a binary correlation matrix
        rho = 0.1*jj
        m = matrix(c(1,rho,rho,1), ncol=2)   
        ## Simulate 1000 pairs
        ## correlation structure
        xtmp = rmvbin(1000, margprob = c(0.5, 0.5), bincorr = m)
        X2$V4 = xtmp[,1]
        Y2=rnorm(n = 1000,mean = 1+((7*X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*X2$V37),sd = 1)
        X2$V4 = xtmp[,2]
        X2$Y2=Y2
        data.example = as.data.frame(X2)
        params[[i]]$data = data.example
        params[[i]]$estimator.args =  list(data = data.example,n = 1000, m = 50)
        params[[i]]$formula = formula1
        if(typeJ)
        {
          params[[i]]$estimator = estimate.logic.lm
        }
        
      }
      
      params[[i]]$cpu=i
      params[[i]]$simul="scenario_Sensitivity5_"
      params[[i]]$simid=jj
    }
    
    #perform garbage collection
    gc()
    #explore Bayesian logic regression on M threads in parallel using the GMJMCMC algorithm
    results=parall.gmj(X = params,FUN = runpar,mc.preschedule = T, mc.cores = M)
    #clean up
    gc()
    wait()
    
    #prepare the data structures for final analysis of the runs
    resa=array(data = 0,dim = c(compmax,M*3))
    post.popul = array(0,M)
    max.popul = array(0,M)
    nulls=NULL
    not.null=1
    #check which threads had non-zero exit status
    for(k in 1:M)
    {
      if(length(results[[k]])<=1||length(results[[k]]$cterm)==0)
      {
        nulls=c(nulls,k)
        next
      }
      else
      {
        not.null = k
      }
      
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
    
    
    
    #delete the unused further variables and perfrom garbage collection
    rm(results)
    gc()
    #renormalize estimates of the marginal inclusion probabilities
    #based on all of the runs
    ml.max=max(max.popul)
    post.popul=post.popul*exp(-ml.max+max.popul)
    p.gen.post=post.popul/sum(post.popul)
    hfinal=hash()
    for(ii in 1:M)
    {
      resa[,ii*3]=p.gen.post[ii]*as.numeric(resa[,ii*3-1])
      resa[length(resa[,ii*3]),ii*3]=p.gen.post[ii]
      if(p.gen.post[ii]>0)
      {
        for(jjj in 1:(length(resa[,ii*3])-1))
        {
          if(resa[jjj,ii*3]>0)
          {
            #print(paste0(ii,"  and ",jj))
            if(as.integer(has.key(hash = hfinal,key =resa[jjj,ii*3-2]))==0)
              hfinal[[resa[jjj,ii*3-2]]]=as.numeric(resa[jjj,ii*3])
            else
              hfinal[[resa[jjj,ii*3-2]]]=hfinal[[resa[jjj,ii*3-2]]]+as.numeric(resa[jjj,ii*3])
          }
          
        }
      }
    }
    
    posteriors=values(hfinal)
    #delete the unused further variables
    clear(hfinal)
    rm(hfinal)
    rm(resa)
    rm(post.popul)
    rm(max.popul)
    #simplify the found trees and their posteriors
    posteriors=as.data.frame(posteriors)
    posteriors=data.frame(X=row.names(posteriors),x=posteriors$posteriors)
    posteriors$X=as.character(posteriors$X)
    tryCatch({
      res1=simplifyposteriors(X = X2,posteriors = posteriors, th,thf)
      write.csv(x =res1,row.names = F,file = paste0("post2eta_",id.sym,"_",j,".csv"))
    },error = function(err){
      print("error")
      print(err)
      write.csv(x =posteriors,row.names = F,file = paste0("posteriors2eta_",id.sym,"_",j,".csv"))
    },finally = {
      
      print(paste0("end simulation ",j))
      
    })
    #delete the unused further variables and perfrom garbage collection
    rm(X2)
    rm(data.example)
    rm(vect)
    rm(params)
    gc()
    print(paste0("end simulation ",j))
  },error = function(err){
    print("error")
    print(err)
    print(paste0("repeat  simulation ",j))
  },finally = {
    
    print(paste0("end simulation ",j))
    
  })
  
  
}


