#####################################################################################################
##################################### Modell 2 #####################################################
####################################################################################################
#read in the package most recent version
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

#dir.create("fractional polynomials", showWarnings = FALSE)
setwd("/mn/sarpanitu/ansatte-u2/aliaksah/RprojectsRstudio/EMJMCMC2016_JAIR/examples/fractional polynomials/")
options(bitmapType="cairo")
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

X = read.csv2(file = "art.csv",sep = ",",dec = ".")[,c(16,1:3,5:8,10:14)]

summary(X)
#number of simulations per scenario
SN =20
#number of observations in the data
n = 250 
#number of covariates
p = 14  

#0.00001
#epsilon = + 1
#x penalty 1
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

g =  function(x) x 

set.seed(040590)

#a parallelizarion function (just ignore it for now)
pinferunemjmcmc = function(n.cores = 4, startseed = 0, mcgmj = mcgmjpar, report.level =  0.5, simplify = F, num.mod.best = 1000, predict = F,  test.data = 1, link.function = function(z)z, runemjmcmc.params)
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
    params[[i]]$cpu=startseed+i
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

#an estimator function that returns the marginal likelihood estimates, parameters and two additional model selection criteria
#it is important as here we can play around with priors 
estimate.gamma.cpen = function(formula,null.mod = 0, data,r = 1.0/exp(2),logn=2,relat1=c("p0","p2","p3","p05","pm05","pm1","pm2"),relat2=c("p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2"))
{
  fparam=NULL
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
  
  tryCatch(capture.output({
    out = glm(formula = formula,data = data, family = gaussian)
    
    
    mlik = (-(BIC(out) -2*log(r)*sum(sj)))/2 - null.mod
    waic = (out$deviance + 2*out$rank) + null.mod
    dic =  (out$deviance + logn*out$rank) + null.mod
    summary.fixed =list(mean = coefficients(out))
    
  }, error = function(err) {
    print(err)
    mlik = -2*null.mod
    waic = 2*null.mod
    dic =  2*null.mod
    summary.fixed =list(mean = array(0,dim=length(fparam)))
  }))
  return(list(mlik = mlik,waic = waic , dic = dic,summary.fixed =summary.fixed))
  
}

dir.create(file.path('resultstmpe100_new__2e45'), showWarnings = FALSE)

#run the simulations
mu = 0.1 + p05(X$x1) + X$x1 + pm05(X$x3) + p0pm05(X$x3) + X$x4a + pm1(X$x5) + p0(X$x6) + X$x8 + X$x10

for(sds in  c(50, 25, 10, 5,1*(0.1)^{0:10}))
{
  
  print(sds) 
  set.seed(040590)
  X[[paste0("y_",sds,collapse = "")]] = rnorm(n =n, mean = mu,sd = sds)
  
}  
  formula1 = as.formula(paste(colnames(X)[1],"~ 1 +",paste0(colnames(X)[-c(1)],collapse = "+")))
  data.example = data.frame(X)
  
  mod.null = logLik(lm(y~1,X))
  
  for(i in 1:SN)
  {
    print(i)
    set.seed(i)
    res = pinferunemjmcmc(predict = F, startseed = i*320,  test.data = 1, n.cores = 32, report.level =  0.2, num.mod.best = 1000 ,simplify = T,link.function = g,
                          runemjmcmc.params = list(formula = formula1,outgraphs=F,data = X,estimator = estimate.gamma.cpen,  
                                                   estimator.args =  list(null.mod = mod.null,data = X, r = 1/n, logn = log(n)),recalc_margin = 249, 
                                                   save.beta = F,interact = T,gen.prob = c(5,0,1,0,0),relations=c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2"),
                                                   relations.prob = rep(0.1,16),
                                                   interact.param=list(deep.method=1,allow_offsprings=3,mutation_rate = 250,last.mutation = 15000, max.tree.size = 2, 
                                                                       Nvars.max =20,p.allow.replace=0.9,p.allow.tree=0.3,p.nor=0.1,p.and = 0.9),n.models =20000,unique = F,
                                                   max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,
                                                   advanced.param = list(
                                                     max.N.glob=as.integer(10),
                                                     min.N.glob=as.integer(5),
                                                     max.N=as.integer(3),
                                                     min.N=as.integer(1),
                                                     printable = F)))
    
    saveRDS(res,file = paste0("resultstmpe100_new__2e45/simulationnorm_",i,"_sd=",sds,".Rdata"))
    print(res$feat.stat)
  }
  
}

#first, analyze on a covariate related level
#tps = c("I(x1)","I(x3)","I(x4a)","I(x5)","I(x6)","I(x8)","I(x10)")

#tps = c("I(p05(I(x1)))","I(x1)","I(pm05(I(x3)))","I(p0pm05(I(x3)))","I(x4a)","I(pm1(I(x5)))","I(p0(I(x6)))","I(x8)","I(x10)")

rhash=hash()
#detection treshold

mmlik = rep(-100000,15)
MLIKtrue =  rep(-100000,15)

rhash=hash()
#detection treshold
alpha = 0.5
mmlik = array(-100000,c(100,15))
MLIKtrue =  array(-100000,c(100,15))

FDR = array(0,c(100,15))
POWER =array(0,c(100,15))
FPS = array(0,c(100,15))
k=0
SN=20
IPOWER = array(-100000,c(100,length(tps),15))
for(l in 1:100)
{
  k=0
  #rhash=hash()
  for(sds in  c(50, 25, 10, 5,1*(0.1)^{0:10}))
  {
    #IPOWER[[l]]=NULL
    k = k + 1
    clear(rhash)
    for(tp in tps)
      rhash[[tp]] = 0
    fdr=0
    ttp=0
    tfp=0
    
    tpcountot = array(0,length(tps))
    set.seed(040590)
    X$y = rnorm(n =n, mean = mu,sd = sds)
    
    set.seed(k*l)
    for(i in sample.int(ifelse(k==9,20,20),SN,replace = T))
    {
      tpcount = array(0,length(tps))
      SNC = SN
      tpi=0
      fpi=0
      
      
      #if(sds == 1 || sds == 10 || sds = 0.1) 
      #  res = readRDS(file = paste0("resultstmpe100_new__2e45/simulation_",i,"_sd=",sds,".Rdata"))
      #else 
      
 
      res = readRDS(file = paste0("resultstmpe100_new__2e45/simulation_",i,"_sd=",sds,".Rdata"))
      
      posterior = as.numeric(res$feat.stat[,2])
      
      
      for(j in 1:length(res$threads.stats))
      {
        if(res$threads.stats[[j]]$cterm>mmlik[l,k])
        {
          mmlik[l,k] = res$threads.stats[[j]]$cterm
        }
      }      
      if(length(res$feat.stat)==0)
      {
        print(paste0("SD: ",sds," run ",i, "did not converge!"))
        SNC = SN-1
        next
      }
      for(j in 1:length(res$feat.stat[,2]))
      {
        
        if(posterior[j]>=alpha)
        {
          TP.which = which(stri_detect_fixed(res$feat.stat[j,1],tps))
          TP.feat = tps[TP.which]
          TP = length( TP.feat > 0)
          
          if(TP)
          { 
            tpcount[TP.which] = tpcount[TP.which] + 1
            if(tpcount[TP.which]>1)
              next
            tpcountot[TP.which] = tpcountot[TP.which] + 1 
          }else
          {
            print(res$feat.stat[j,1])
          }
          
          if(TP)
          {
            if(length(as.numeric(rhash[[TP.feat]]))>0)
              tmpc = as.numeric(rhash[[TP.feat]])
            else
            { 
              tmpc = 0 
              rhash[[TP.feat]]=tmpc+1
            }
            rhash[[TP.feat]]=tmpc+1
            
            tpi=tpi+1
            ttp= ttp + 1
          }else{
            fpi=fpi+1
            tfp= tfp + 1
          }
          
        }
      }
      
      if((fpi+tpi)>0)
        fdr=fdr+(fpi/(fpi+tpi))
      #print(fpi)
      #print(tpi)
    }
    
    IPOWER[l,1:length(tps),k] = tpcountot/SNC
    if(l==1)
      MLIKtrue[l,k] = estimate.gamma.cpen(data = X, r = 1/n, logn = log(n), 
                                      formula = y ~ 1+ I(p05(I(x1))) + I(x1) + I(pm05(I(x3))) + I(p0pm05(I(x3))) + I(x4a) + I(pm1(I(x5))) + I(p0(I(x6))) + I(x8) + I(x10) )$mlik
    
    #print the FDR out
    print("Results for:")
    print(sds)
    print("MLIKs:")
    print(c(mmlik[l,k],MLIKtrue[l,k]))
    print("FDR:")
    FDR[l,k] = fdr/SNC
    print(FDR[l,k])
    print("Power:")
    POWER[l,k] = ttp/(SNC*length(tps))
    print(POWER[l,k])
    print("FPS:")
    FPS[l,k] = tfp/(SNC)
    print(FPS[l,k])
    print(rhash)
    
  }
}

L = 100
sd = c(100, 50, 25, 10, 5,1*(0.1)^{0:10})
k = k + 1

mm0 = c(0,unlist(lapply(X = 1:15,FUN = function(k)median(POWER[1:100,k]))))
ll0 = c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(x=POWER[1:100,k],probs = 0.025))))
uu10 = c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(x=POWER[1:100,k],probs = 0.975))))

write.csv(file = "powerssoft.csv",x = rbind(mm0,ll0,uu10))

pow.mfp = read.csv("Power_MFP.csv")
pow.bfp = read.csv("Power_BFP_U.csv")
pow.bfps = read.csv("Power_BFP_S.csv")
pow.bfpd = read.csv("Power_BFP_D.csv")
pow.ideal = read.csv("Power_IDEAL_strict.csv")[,-1]

plot(y = mm0,x = 1:k,type = "l",ylab = "",ylim = c(min(ll0),1),main="TPR (strict)",xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
#lines(y =ll0 ,x = 1:k,type = "l",lty = "dashed")
#lines(y =uu10 ,x = 1:k,type = "l",lty = "dashed")
lines(y = pow.ideal[1,],x = 1:k,"l",col = 7)
#lines(y = pow.ideal[2,],x = 1:k,"l",col = 7,lty = "dashed")
#lines(y = pow.ideal[3,],x = 1:k,"l",col = 7,lty = "dashed")
lines(y = pow.mfp[4,],x = 1:k,"l",col = 2)
#lines(y = pow.mfp[5,],x = 1:k,"l",col = 2,lty = "dashed")
#lines(y = pow.mfp[6,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = pow.bfp[4,],x = 1:k,"l",col = 3)
#lines(y = pow.bfp[5,],x = 1:k,"l",col = 3,lty = "dashed")
#lines(y = pow.bfp[6,],x = 1:k,"l",col = 3,lty = "dashed")
lines(y = pow.bfps[4,],x = 1:k,"l",col = 4)
#lines(y = pow.bfps[5,],x = 1:k,"l",col = 4,lty = "dashed")
#lines(y = pow.bfps[6,],x = 1:k,"l",col = 4,lty = "dashed")
lines(y = pow.bfpd[4,],x = 1:k,"l",col = 6)
#lines(y = pow.bfpd[5,],x = 1:k,"l",col = 6,lty = "dashed")
#lines(y = pow.bfpd[6,],x = 1:k,"l",col = 6,lty = "dashed")
legend(x = 11.22, y = 0.320, cex = 0.8, c("BGNLM_FP_IDEAL","BGNLM_FP", "MFP", "BFP_U","BFP_S","BFP_D"), col = c(7,1,2,3,4,6),
       text.col = "black", lty = c(1,1,1,1,1), 
       merge = TRUE, bg = "white", border = "white",box.col = "white")

pow.ideal = read.csv("Power_IDEAL_soft.csv")[,-1]

plot(y = mm0,x = 1:k,type = "l",ylab = "",ylim = c(min(ll0),max(uu10)),main="TPR (soft)",xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
#lines(y =ll0 ,x = 1:k,type = "l",lty = "dashed")
#lines(y =uu10 ,x = 1:k,type = "l",lty = "dashed")
lines(y = pow.mfp[1,],x = 1:k,"l",col = 2)
#lines(y = pow.mfp[2,],x = 1:k,"l",col = 2,lty = "dashed")
#lines(y = pow.mfp[3,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = pow.ideal[1,],x = 1:k,"l",col = 7)
#lines(y = pow.ideal[2,],x = 1:k,"l",col = 7,lty = "dashed")
#lines(y = pow.ideal[3,],x = 1:k,"l",col = 7,lty = "dashed")
lines(y = pow.bfp[1,],x = 1:k,"l",col = 3)
#lines(y = pow.bfp[2,],x = 1:k,"l",col = 3,lty = "dashed")
#lines(y = pow.bfp[3,],x = 1:k,"l",col = 3,lty = "dashed")
lines(y = pow.bfps[1,],x = 1:k,"l",col = 4)
#lines(y = pow.bfps[2,],x = 1:k,"l",col = 4,lty = "dashed")
#lines(y = pow.bfps[3,],x = 1:k,"l",col = 4,lty = "dashed")
lines(y = pow.bfpd[1,],x = 1:k,"l",col = 6)
#lines(y = pow.bfpd[3,],x = 1:k,"l",col = 6,lty = "dashed")
#lines(y = pow.bfpd[3,],x = 1:k,"l",col = 6,lty = "dashed")
legend(x = 11.22, y = 0.370, c("BGNLM_BFP_PRL", "FFP", "BFP_U","BFP_S","BFP_D"), col = c(1,2,3,4,6),
       text.col = "black", lty = c(1,1,1,1,1), 
       merge = TRUE, bg = "white", border = "white",box.col = "white")




fdr.mfp = read.csv("fdr_MFP.csv")
fdr.bfp = read.csv("fdr_BFP_U.csv")

fdr.bfps = read.csv("fdr_BFP_S.csv")
fdr.bfpd = read.csv("fdr_BFP_D.csv")

fdr.ideal = read.csv("FDR_IDEAL_strict.csv")[,-1]



mm1 = c(0,unlist(lapply(X = 1:15,FUN = function(k)median(FDR[1:100,k]))))
ll1 = c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(x=FDR[1:100,k],probs = 0.025))))
uu11 = c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(x=FDR[1:100,k],probs = 0.975))))

write.csv(file = "fdrsoft.csv",x = rbind(mm1,ll1,uu11))



plot(y =mm1,x = 1:k,type = "l",ylab = "",main= "FDR (strict)",ylim = c(min(ll1),max(fdr.bfp)),xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
#lines(y =ll1 ,x = 1:k,type = "l",lty = "dashed")
#lines(y =uu11 ,x = 1:k,type = "l",lty = "dashed")
lines(y = fdr.ideal[1,],x = 1:k,"l",col = 7)
#lines(y = fdr.ideal[2,],x = 1:k,"l",col = 7,lty = "dashed")
#lines(y = fdr.ideal[3,],x = 1:k,"l",col = 7,lty = "dashed")
lines(y = fdr.mfp[4,],x = 1:k,"l",col = 2)
#lines(y = fdr.mfp[5,],x = 1:k,"l",col = 2,lty = "dashed")
#lines(y = fdr.mfp[6,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = fdr.bfp[4,],x = 1:k,"l",col = 3)
#lines(y = fdr.bfp[5,],x = 1:k,"l",col = 3,lty = "dashed")
#lines(y = fdr.bfp[6,],x = 1:k,"l",col = 3,lty = "dashed")
lines(y = fdr.bfps[4,],x = 1:k,"l",col = 4)
#lines(y = fdr.bfps[5,],x = 1:k,"l",col = 4,lty = "dashed")
#lines(y = fdr.bfps[6,],x = 1:k,"l",col = 4,lty = "dashed")
lines(y = fdr.bfpd[4,],x = 1:k,"l",col = 6)
#lines(y = fdr.bfpd[5,],x = 1:k,"l",col = 6,lty = "dashed")
#lines(y = fdr.bfpd[6,],x = 1:k,"l",col = 6,lty = "dashed")
legend(x = 11.22, y = 0.350, c("BGNLM_BFP_PRL", "FFP", "BFP_U","BFP_S","BFP_D"), col = c(1,2,3,4,6),
       text.col = "black", lty = c(1,1,1,1,1), 
       merge = TRUE, bg = "white", border = "white",box.col = "white")

fdr.ideal = read.csv("FDR_IDEAL_soft.csv")[,-1]


plot(y =mm1,x = 1:k,type = "l",ylab = "",main= "FDR (soft)",ylim = c(min(ll1),max(fdr.bfp)),xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
#lines(y =ll1 ,x = 1:k,type = "l",lty = "dashed")
#lines(y =uu11 ,x = 1:k,type = "l",lty = "dashed")
lines(y = fdr.ideal[1,],x = 1:k,"l",col = 7)
#lines(y = fdr.ideal[2,],x = 1:k,"l",col = 7,lty = "dashed")
#lines(y = fdr.ideal[3,],x = 1:k,"l",col = 7,lty = "dashed")
lines(y = fdr.mfp[1,],x = 1:k,"l",col = 2)
#lines(y = fdr.mfp[2,],x = 1:k,"l",col = 2,lty = "dashed")
#lines(y = fdr.mfp[3,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = fdr.bfp[1,],x = 1:k,"l",col = 3)
#lines(y = fdr.bfp[2,],x = 1:k,"l",col = 3,lty = "dashed")
#lines(y = fdr.bfp[3,],x = 1:k,"l",col = 3,lty = "dashed")
lines(y = fdr.bfps[1,],x = 1:k,"l",col = 4)
#lines(y = fdr.bfps[2,],x = 1:k,"l",col = 4,lty = "dashed")
#lines(y = fdr.bfps[3,],x = 1:k,"l",col = 4,lty = "dashed")
lines(y = fdr.bfpd[1,],x = 1:k,"l",col = 6)
#lines(y = fdr.bfpd[2,],x = 1:k,"l",col = 6,lty = "dashed")
#lines(y = fdr.bfpd[3,],x = 1:k,"l",col = 6,lty = "dashed")
legend(x = 11.22, y = 0.950, c("BGNLM_BFP_PRL", "FFP", "BFP_U","BFP_S","BFP_D"), col = c(1,2,3,4,6),
       text.col = "black", lty = c(1,1,1,1,1), 
       merge = TRUE, bg = "white", border = "white",box.col = "white")




mm2 = unlist(lapply(X = 1:15,FUN = function(k)median(FPS[1:100,k])))
ll2 = unlist(lapply(X = 1:15,FUN = function(k)quantile(x=FPS[1:100,k],probs = 0.025)))
uu12 = unlist(lapply(X = 1:15,FUN = function(k)quantile(x=FPS[1:100,k],probs = 0.975)))

plot(y =mm2,x = 1:k,type = "l",ylab = "",main= "FPS",ylim = c(min(ll2),max(uu12)),xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
lines(y =ll2 ,x = 1:k,type = "l",lty = "dashed")
lines(y =uu12 ,x = 1:k,type = "l",lty = "dashed")



mmlik = rep(-100000,16)
MLIKtrue =  rep(-100000,16)

rhash=hash()
#detection treshold
alpha = 0.5
mmlik = array(-100000,c(100,16))
MLIKtrue =  array(-100000,c(100,16))
RSQRtrue =  array(-100000,c(100,16))

X = read.csv2(file = "art.csv",sep = ",",dec = ".")[,c(16,1:3,5:8,10:14)]
for(l in 1:100)
{
  k=0
  #rhash=hash()
  for(sds in  c(100,50, 25, 10, 5,1*(0.1)^{0:10}))
  {
    #IPOWER[[l]]=NULL
    k = k + 1
   
    set.seed(040590)
    X$y = rnorm(n =n, mean = mu,sd = sds)
    
    set.seed(k*l)
    if(k>1)
    {for(i in sample.int(ifelse(k==9,20,20),SN,replace = T))
      {
      
        res = readRDS(file = paste0("resultstmpe100_new__2e45/simulation_",i,"_sd=",sds,".Rdata"))
        posterior = as.numeric(res$feat.stat[,2])
        
        
        for(j in 1:length(res$threads.stats))
        {
          if(res$threads.stats[[j]]$cterm>mmlik[l,k])
          {
            mmlik[l,k] = res$threads.stats[[j]]$cterm
          }
        }      
      }
    }
    else
    {
      mmlik[l,k] = estimate.gamma.cpen(data = X, null.mod = -500,r = 1/n, logn = log(n), 
                                     formula = y ~ 1)$mlik
      #ltrue = summary(lm(data = X, formula = y ~ 1+ I(p05(I(x1))) + I(x1) + I(pm05(I(x3))) + I(p0pm05(I(x3))) + I(x4a) + I(pm1(I(x5))) + I(p0(I(x6))) + I(x8) + I(x10)))
      #RSQRtrue[l,k] = ltrue$r.squared
    }
    if(l==1)
    {  
      MLIKtrue[l,k] = estimate.gamma.cpen(data = X,null.mod = -500, r = 1/n, logn = log(n), 
                                          formula = y ~ 1+ I(p05(I(x1))) + I(x1) + I(pm05(I(x3))) + I(p0pm05(I(x3))) + I(x4a) + I(pm1(I(x5))) + I(p0(I(x6))) + I(x8) + I(x10) )$mlik
      ltrue = summary(lm(data = X, formula = y ~ 1+ I(p05(I(x1))) + I(x1) + I(pm05(I(x3))) + I(p0pm05(I(x3))) + I(x4a) + I(pm1(I(x5))) + I(p0(I(x6))) + I(x8) + I(x10)))
      RSQRtrue[l,k] = ltrue$r.squared
    }
    #print the FDR out
    print("Results for:")
    print(sds)
    print("MLIKs:")
    print(c(mmlik[l,k],MLIKtrue[l,k]))

    
  }
}

k=16
mm3 = unlist(lapply(X = 1:16,FUN = function(k)median(mmlik[1:100,k])))
ll3 = unlist(lapply(X = 1:16,FUN = function(k)quantile(x=mmlik[1:100,k],probs = 0.025)))
uu3 = unlist(lapply(X = 1:16,FUN = function(k)quantile(x=mmlik[1:100,k],probs = 0.975)))

plot(y =mm3,x = 1:k,type = "l",ylab = "",main= "Marginal log likelihood",ylim = c(min(ll3),max(c(uu3,max(MLIKtrue[1,k:1])))+2000),xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
axis(1, 1:k, as.character(round(RSQRtrue[1,1:k],4)),cex.axis=1.5,las = 2,line = -15.8)
#lines(y =ll3 ,x = 1:k,type = "l",lty = "dashed")
#lines(y =uu13 ,x = 1:k,type = "l",lty = "dashed")
lines(MLIKtrue[1,1:k],col=2)
#lines(y =ll31 ,x = 1:k,type = "l",col= 3,lty = "dashed")
#lines(y =uu131 ,x = 1:k,type = "l",col= 3,lty = "dashed")
#lines(mm31,col= 3)


IPOWER1=IPOWER
colnames(IPOWER1) = tps
for(i in 1:7)
{  
  plot(unlist(lapply(X = 9:1,FUN = function(k)median(IPOWER1[,i,k]))),x = 1:k,ylab = "", ylim = c(0,1), type = "l",main = tps[i],xlab="",xaxt = "n",cex.axis=1.5)
  axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
  lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER1[,i,k],probs = 0.025))),x = 1:k,type = "l",lty = "dashed")
  lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER1[,i,k],probs = 0.975))),x = 1:k,type = "l",lty = "dashed")
}


pow.mfp.x3 = read.csv("x3_MFP.csv")
pow.bfp.x3 = read.csv("x3_BFP_U.csv")
pow.bfps.x3 = read.csv("x3_BFP_S.csv")
pow.bfpd.x3 = read.csv("x3_BFP_D.csv")

pow.ideal.x3 = read.csv("Power x3 strict.csv")[,-1]

pow.ideal.x3 = read.csv("Power x3 soft.csv")[,-1]

IPOWER1=IPOWER
colnames(IPOWER1) = tps
i = 3
for(i in c(3))
{  
  plot(c(0,unlist(lapply(X = 1:15,FUN = function(k)median(IPOWER1[,i,k])))),x = 1:k,ylab = "", ylim = c(0,1), type = "l",main = tps[i],xlab="",xaxt = "n",cex.axis=1.5)
  axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
  #lines(c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(IPOWER1[,i,k],probs = 0.025)))),x = 1:k,type = "l",lty = "dashed")
  #lines(c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(IPOWER1[,i,k],probs = 0.975)))),x = 1:k,type = "l",lty = "dashed")
  lines(y = pow.ideal.x3[1,],x = 1:k,"l",col = 7)
  #lines(y = pow.ideal.x3[2,],x = 1:k,"l",col = 7,lty = "dashed")
  #lines(y = pow.ideal.x3[3,],x = 1:k,"l",col = 7,lty = "dashed")
  lines(y = pow.mfp.x3[4,],x = 1:k,"l",col = 2)
  #lines(y = pow.mfp.x3[5,],x = 1:k,"l",col = 2,lty = "dashed")
  #lines(y = pow.mfp.x3[6,],x = 1:k,"l",col = 2,lty = "dashed")
  lines(y = pow.bfp.x3[4,],x = 1:k,"l",col = 3)
  #lines(y = pow.bfp.x3[5,],x = 1:k,"l",col = 3,lty = "dashed")
  #lines(y = pow.bfp.x3[6,],x = 1:k,"l",col = 3,lty = "dashed")
  lines(y = pow.bfps.x3[4,],x = 1:k,"l",col = 4)
  #lines(y = pow.bfps.x3[5,],x = 1:k,"l",col = 4,lty = "dashed")
  #lines(y = pow.bfps.x3[6,],x = 1:k,"l",col = 4,lty = "dashed")
  lines(y = pow.bfpd.x3[4,],x = 1:k,"l",col = 6)
  #lines(y = pow.bfpd.x3[5,],x = 1:k,"l",col = 6,lty = "dashed")
  #lines(y = pow.bfpd.x3[6,],x = 1:k,"l",col = 6,lty = "dashed")
  legend(x = 11.22, y = 0.320, cex = 0.8, c("BGNLM_FP_IDEAL","BGNLM_FP", "MFP", "BFP_U","BFP_S","BFP_D"), col = c(7,1,2,3,4,6),
         text.col = "black", lty = c(1,1,1,1,1), 
         merge = TRUE, bg = "white", border = "white",box.col = "white")
  i=2
  plot(c(0,unlist(lapply(X = 1:15,FUN = function(k)median(IPOWER1[,i,k])))),x = 1:k,ylab = "", ylim = c(0,1), type = "l",main = tps[i],xlab="",xaxt = "n",cex.axis=1.5)
  axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
  #lines(c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(IPOWER1[,i,k],probs = 0.025)))),x = 1:k,type = "l",lty = "dashed")
  #lines(c(0,unlist(lapply(X = 1:15,FUN = function(k)quantile(IPOWER1[,i,k],probs = 0.975)))),x = 1:k,type = "l",lty = "dashed")
  lines(y = pow.ideal.x3[1,],x = 1:k,"l",col = 7)
  #lines(y = pow.ideal.x3[2,],x = 1:k,"l",col = 7,lty = "dashed")
  #lines(y = pow.ideal.x3[3,],x = 1:k,"l",col = 7,lty = "dashed")
  lines(y = pow.mfp.x3[1,],x = 1:k,"l",col = 2)
  #lines(y = pow.mfp.x3[2,],x = 1:k,"l",col = 2,lty = "dashed")
  #lines(y = pow.mfp.x3[3,],x = 1:k,"l",col = 2,lty = "dashed")
  lines(y = pow.bfp.x3[1,],x = 1:k,"l",col = 3)
  #lines(y = pow.bfp.x3[2,],x = 1:k,"l",col = 3,lty = "dashed")
  #lines(y = pow.bfp.x3[3,],x = 1:k,"l",col = 3,lty = "dashed")
  lines(y = pow.bfps.x3[1,],x = 1:k,"l",col = 4)
  #lines(y = pow.bfps.x3[2,],x = 1:k,"l",col = 4,lty = "dashed")
  #lines(y = pow.bfps.x3[3,],x = 1:k,"l",col = 4,lty = "dashed")
  lines(y = pow.bfpd.x3[1,],x = 1:k,"l",col = 6)
  #lines(y = pow.bfpd.x3[2,],x = 1:k,"l",col = 6,lty = "dashed")
  #lines(y = pow.bfpd.x3[3,],x = 1:k,"l",col = 6,lty = "dashed")
}



tps = c("I(p05(I(x1)))","I(x1)","I(pm05(I(x3)))","I(p0pm05(I(x3)))","I(x4a)","I(pm1(I(x5)))","I(p0(I(x6)))","I(x8)","I(x10)")

rhash=hash()
#detection treshold

mmlik = rep(-100000,11)
MLIKtrue =  rep(-100000,11)
RSQRtrue =  array(-100000,c(100,16))
rhash=hash()
#detection treshold

mmlik = array(-100000,c(100,11))
MLIKtrue =  array(-100000,c(100,11))

FDR = array(0,c(100,11))
POWER =array(0,c(100,11))
FPS = array(0,c(100,11))
k=0
SN=20
IPOWER = array(-100000,c(100,length(tps),11))
for(l in 1:100)
{
  k=0
  #rhash=hash()
  for(sds in 1*(0.1)^{10:0})
  {
    #IPOWER[[l]]=NULL
    k = k + 1
    clear(rhash)
    for(tp in tps)
      rhash[[tp]] = 0
    fdr=0
    ttp=0
    tfp=0
    
    tpcountot = array(0,length(tps))
    set.seed(040590)
    X$y = rnorm(n =n, mean = mu,sd = sds)
    
    set.seed(k*l)
    for(i in sample.int(ifelse(k==9,20,20),SN,replace = T))
    {
      tpcount = array(0,length(tps))
      SNC = SN
      tpi=0
      fpi=0
      
      
      res = readRDS(file = paste0("resultstmpe100_new__2e45/simulation_",i,"_sd=",sds,".Rdata"))
      posterior = as.numeric(res$feat.stat[,2])
      
      
      for(j in 1:length(res$threads.stats))
      {
        if(res$threads.stats[[j]]$cterm>mmlik[l,k])
        {
          mmlik[l,k] = res$threads.stats[[j]]$cterm
        }
      }      
      if(length(res$feat.stat)==0)
      {
        print(paste0("SD: ",sds," run ",i, "did not converge!"))
        SNC = SN-1
        next
      }
      for(j in 1:length(res$feat.stat[,2]))
      {
        
        if(posterior[j]>=alpha)
        {
          TP.which = which(stri_detect_fixed(res$feat.stat[j,1],tps))
          TP.feat = tps[TP.which]
          TP = length( TP.feat > 0)
          
          if(TP)
          { 
            tpcount[TP.which] = tpcount[TP.which] + 1
            if(tpcount[TP.which]>1)
              next
            tpcountot[TP.which] = tpcountot[TP.which] + 1 
          }else
          {
            print(res$feat.stat[j,1])
          }
          
          if(TP)
          {
            if(length(as.numeric(rhash[[TP.feat]]))>0)
              tmpc = as.numeric(rhash[[TP.feat]])
            else
            { 
              tmpc = 0 
              rhash[[TP.feat]]=tmpc+1
            }
            rhash[[TP.feat]]=tmpc+1
            
            tpi=tpi+1
            ttp= ttp + 1
          }else{
            fpi=fpi+1
            tfp= tfp + 1
          }
          
        }
      }
      
      if((fpi+tpi)>0)
        fdr=fdr+(fpi/(fpi+tpi))
      #print(fpi)
      #print(tpi)
    }
    
    IPOWER[l,1:length(tps),k] = tpcountot/SNC
    if(l==1)
      MLIKtrue[l,k] = estimate.gamma.cpen(data = X, r = 1/n, logn = log(n), 
                                          formula = y ~ 1+ I(p05(I(x1))) + I(x1) + I(pm05(I(x3))) + I(p0pm05(I(x3))) + I(x4a) + I(pm1(I(x5))) + I(p0(I(x6))) + I(x8) + I(x10) )$mlik
    
    #print the FDR out
    print("Results for:")
    print(sds)
    print("MLIKs:")
    print(c(mmlik[l,k],MLIKtrue[l,k]))
    print("FDR:")
    FDR[l,k] = fdr/SNC
    print(FDR[l,k])
    print("Power:")
    POWER[l,k] = ttp/(SNC*length(tps))
    print(POWER[l,k])
    print("FPS:")
    FPS[l,k] = tfp/(SNC)
    print(FPS[l,k])
    print(rhash)
    
  }
}

sd = 1*(0.1)^{0:10}

mm10 = unlist(lapply(X = 11:1,FUN = function(k)median(POWER[1:100,k])))
ll10 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=POWER[1:100,k],probs = 0.025)))
uu10 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=POWER[1:100,k],probs = 0.975)))

plot(y =mm10,x = 1:k,type = "l",ylab = "",ylim = c(min(pow.bfp[5,]),max(uu10)),main="Power (strict)",xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
lines(y =ll10 ,x = 1:k,type = "l",lty = "dashed")
lines(y =uu10 ,x = 1:k,type = "l",lty = "dashed")
lines(y = pow.mfp[4,],x = 1:k,"l",col = 2)
lines(y = pow.mfp[5,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = pow.mfp[6,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = pow.bfp[4,],x = 1:k,"l",col = 3)
lines(y = pow.bfp[5,],x = 1:k,"l",col = 3,lty = "dashed")
lines(y = pow.bfp[6,],x = 1:k,"l",col = 3,lty = "dashed")
legend(x = 7.22, y = 0.675, c("BGNLM", "MFP", "BFP","","95% CI",""), col = c(1, 2, 3,1,2,3),
       text.col = "black", lty = c(1, 1, 1,2,2,2), 
       merge = TRUE, bg = "gray90")

lines(y =mm0 ,x = 1:k,type = "l",col=2)
lines(y =ll0 ,x = 1:k,type = "l",lty = "dashed",col=2)
lines(y =uu0 ,x = 1:k,type = "l",lty = "dashed",col=2)


mm11 = unlist(lapply(X = 11:1,FUN = function(k)median(FDR[1:100,k])))
ll11 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=FDR[1:100,k],probs = 0.025)))
uu11 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=FDR[1:100,k],probs = 0.975)))

plot(y =mm11,x = 1:k,type = "l",ylab = "",main= "FDR (strict)",ylim = c(min(ll1),max(fdr.bfp)),xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
lines(y =ll11 ,x = 1:k,type = "l",lty = "dashed")
lines(y =uu11 ,x = 1:k,type = "l",lty = "dashed")
lines(y = fdr.mfp[4,],x = 1:k,"l",col = 2)
lines(y = fdr.mfp[5,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = fdr.mfp[6,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = fdr.bfp[4,],x = 1:k,"l",col = 3)
lines(y = fdr.bfp[5,],x = 1:k,"l",col = 3,lty = "dashed")
lines(y = fdr.bfp[6,],x = 1:k,"l",col = 3,lty = "dashed")
legend(x = 7.22, y = 0.333, c("BGNLM", "MFP", "BFP","","95% CI",""), col = c(1, 2, 3,1,2,3),
       text.col = "black", lty = c(1, 1, 1,2,2,2), 
       merge = TRUE, bg = "gray90")


#lines(y =mm1 ,x = 1:k,type = "l",col=2)
#lines(y =ll1 ,x = 1:k,type = "l",lty = "dashed",col=2)
#lines(y =uu1 ,x = 1:k,type = "l",lty = "dashed",col=2)

mm12 = unlist(lapply(X = 11:1,FUN = function(k)median(FPS[1:100,k])))
ll12 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=FPS[1:100,k],probs = 0.025)))
uu12 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=FPS[1:100,k],probs = 0.975)))

plot(y =mm12,x = 1:k,type = "l",ylab = "",main= "FPS",ylim = c(min(ll2),max(uu12)),xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
lines(y =ll12 ,x = 1:k,type = "l",lty = "dashed")
lines(y =uu12 ,x = 1:k,type = "l",lty = "dashed")
lines(y =mm2 ,x = 1:k,type = "l",col=2)
lines(y =ll2 ,x = 1:k,type = "l",lty = "dashed",col=2)
lines(y =uu2 ,x = 1:k,type = "l",lty = "dashed",col=2)

mm13 = unlist(lapply(X = 11:1,FUN = function(k)median(mmlik[1:100,k])))
ll13 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=mmlik[1:100,k],probs = 0.025)))
uu13 = unlist(lapply(X = 11:1,FUN = function(k)quantile(x=mmlik[1:100,k],probs = 0.975)))

plot(y =mm13,x = 1:k,type = "l",ylab = "",main= "Log marginal posterior",ylim = c(min(ll13),max(c(uu13,max(MLIKtrue[1,k:1])))),xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
lines(y =ll13 ,x = 1:k,type = "l",lty = "dashed")
lines(y =uu13 ,x = 1:k,type = "l",lty = "dashed")
lines(MLIKtrue[1,k:1],col=2)
legend(x = 0.84, y = 4000, c("Found","95% CI","True"), col = c(1, 1, 4),
       text.col = "black", lty = c(1, 2, 1), 
       merge = TRUE, bg = "gray90")

IPOWER2= IPOWER
colnames(IPOWER2) = tps
for(i in 1:9)
{  
  plot(unlist(lapply(X = 9:1,FUN = function(k)median(IPOWER2[,i,k]))),x = 1:k,ylab = "", ylim = c(0,1), type = "l",main = tps[i],xlab="",xaxt = "n",cex.axis=1.5)
  axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
  lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER2[,i,k],probs = 0.025))),x = 1:k,type = "l",lty = "dashed")
  lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER2[,i,k],probs = 0.975))),x = 1:k,type = "l",lty = "dashed")
  
}

plot(unlist(lapply(X = 9:1,FUN = function(k)median(IPOWER2[,3,k]))),x = 1:k,ylab = "", ylim = c(0,1), type = "l",main = "Power x3^(-0.5) (strict)",xlab="",xaxt = "n",cex.axis=1.5)
axis(1, at=1:k, labels=sd,cex.axis=1.5,las = 2)
lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER2[,3,k],probs = 0.025))),x = 1:k,type = "l",lty = "dashed")
lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER2[,3,k],probs = 0.975))),x = 1:k,type = "l",lty = "dashed")
#lines(unlist(lapply(X = 9:1,FUN = function(k)median(IPOWER2[,4,k]))),x = 1:k,type = "l",col=4)
#lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER2[,4,k],probs = 0.025))),x = 1:k,type = "l",lty = "dashed",col=4)
#lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER2[,4,k],probs = 0.975))),x = 1:k,type = "l",lty = "dashed",col=4)
lines(y = pow.mfp.x3[4,],x = 1:k,"l",col = 2)
lines(y = pow.mfp.x3[5,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = pow.mfp.x3[6,],x = 1:k,"l",col = 2,lty = "dashed")
lines(y = pow.bfp.x3[4,],x = 1:k,"l",col = 3)
lines(y = pow.bfp.x3[5,],x = 1:k,"l",col = 3,lty = "dashed")
lines(y = pow.bfp.x3[6,],x = 1:k,"l",col = 3,lty = "dashed")
legend(x = 7.22, y = 0.5, c("BGNLM", "MFP", "BFP","","95% CI",""), col = c(1, 2, 3,1,2,3),
       text.col = "black", lty = c(1, 1, 1,2,2,2), 
       merge = TRUE, bg = "gray90")

lines(unlist(lapply(X = 9:1,FUN = function(k)median(IPOWER1[,2,k]))),x = 1:k,type = "l",col=2)
lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER1[,2,k],probs = 0.025))),x = 1:k,type = "l",lty = "dashed",col=2)
lines(unlist(lapply(X = 9:1,FUN = function(k)quantile(IPOWER1[,2,k],probs = 0.975))),x = 1:k,type = "l",lty = "dashed",col=2)

cor(p05(X$x1),X$x1)
cor(pm05(X$x3),p0pm05(X$x3))
