library(BMA)
library(pec)
options(bitmapType="cairo")
# get data
download.file('https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/Dateien/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/gbsg_br_ca.zip',
              'gbsg_br_ca.zip')
data <- read.csv(unz('gbsg_br_ca.zip',
                     'gbsg_br_ca/gbsg_br_ca.csv'),
                 header = TRUE)
#system('rm whitehall1.zip')


head(data)
time <- data$rectime
cens <- data$censrec
X <- data[, c(2:4, 6:8, 10:12)]
data <- cbind(time, cens, X)

i = 1
for(nam in names(data))
{ 
  data[[nam]] = as.numeric(data[[nam]])
}

rm(X)

gc()

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

#dir.create("fractional polynomials", showWarnings = FALSE)
setwd("/mn/sarpanitu/ansatte-u2/aliaksah/RprojectsRstudio/EMJMCMC2016_JAIR/examples/fractional polynomials/")

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)


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



data.test = data[501:686,]
data.example = data[1:500,]

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
estimate.bic.surv.cpen = function(formula, data,r = 1.0/500.0,logn=log(500.0),relat1=c("p0","p2","p3","p05","pm05","pm1","pm2"),relat2=c("p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2"))
{
  #print(names(data))
  mlik = NULL
  toret = NULL
  fparam=1
  tryCatch({
    
    fmla.proc=as.character(formula)[2:3]
    fobserved = fmla.proc[1]
    fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
    fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
    fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
    
    if(max(table(fparam))>1)
      return(list(mlik = -10000,waic =  10000 , dic = 10000, summary.fixed =list(mean = array(NA,dim=length(fparam)))))
    
    sj=sum(stri_count_fixed(str = fparam, pattern = "*"))
    sp = sum(stri_count_fixed(str = fparam, pattern = "+"))
    sj=sj + sp
    for(rel in relat1)
      sj=sj+log(2)*(stri_count_fixed(str = fparam, pattern = rel))
    for(rel in relat2)
      sj=sj+log(4)*(stri_count_fixed(str = fparam, pattern = rel))
    
    #print(fparam)
    if(length(fparam)==1)
    {  
      #out = coxph(Surv(time,cens)~1, data = data)
      return(list(mlik = log(r)*sum(sj),waic = -log(r)*sum(sj) , dic =  -log(r)*sum(sj),summary.fixed =list(mean = 0))) 
    }
    formula1 = as.formula(paste0("Surv(time,cens)","~",fmla.proc[2]))
    
    
    
    
    
    if(length(fparam)==2)
    {
      out = coxph(formula1, data = data)
      out.null = coxph(Surv(time,cens)~1, data = data)
      mlik = out$loglik[1] - out.null$loglik[1]+ 2*log(r)*sum(sj)
      toret = list(mlik = mlik,waic = -mlik , dic = -mlik,summary.fixed =list(mean =  c(0,out$coefficients)))
    }else{
      
      #print(formula1)
      out = bic.surv(formula1, data = data,
                     factor.type = F,prior.param = 1,method = "breslow", iter.max = 30)
      
      mlik = ((-out$bic + log(r)*length(out$mle) +2*log(r)*sum(sj)))/2
      
      toret = list(mlik = mlik,waic = -mlik , dic =  -mlik,summary.fixed =list(mean =  c(0,out$mle))) 
    }
  }, error = function(err) {
    #print(err)
    toret = NULL
    mlik = NULL
  })
  
  if(is.null(toret)||is.null(mlik)||is.na(mlik))
    return(list(mlik = -10000,waic =  10000 , dic = 10000, summary.fixed =list(mean = array(0,dim=length(fparam)))))
  
  return(toret)
}


predict.surv.thread <- function(thread.id,data.train,data.test)
{
  
  transf = thread.id$fparam
  coeff = thread.id$betas
  model.weights = thread.id$mliks-min(thread.id$mliks)
  model.weights = exp(model.weights)/sum(exp(model.weights))
  # coeff: matrix with the regression coefficients (columns) for each model (row)
  # data: full data
  # train: indexes of the train data
  # model_weights: vector of dimension nrow(coeff) with the model weights
  # transf: list of FP, if null, assume that x.test has already the transformed variables
  
  # if the test data are given as variables and not FP, transform
  if(length(transf) > 0)
  {
    form.data <- as.formula(paste0("time~",paste0(transf,collapse = "+")))
    x.train <- model.matrix(form.data,data.train)
    x.test <-  model.matrix(form.data,data.test)
  }else{
    # for now, suppose x.test has the right transformations
    x.train <- as.matrix(data.train)
    x.test <- as.matrix(data.test)
  }
  
  # if the format of coeff has NA instead of 0, use 0
  if (sum(is.na(coeff)) > 0)  coeff[is.na(coeff)] <- 0
  
  # scores
  average.lin.pred.train <- apply(x.train%*%t(coeff), 1, mean, na.rm = TRUE,  w = model.weights)
  
  # average the risk scores (step b' page 441 )
  average.lin.pred.test <- apply(x.test%*%t(coeff), 1, mean, na.rm = TRUE, w = model.weights)
  
  return(list(train = average.lin.pred.train,test = average.lin.pred.test))
  
}


predict.surv.bma <- function(res.object,data.train,data.test)
{
  converged <- NULL
  for(i in 1:length(res.object$threads.stats))
  {
    if(length(res.object$threads.stats[[i]])==1)
    {
      warning(paste0("Thread ",i," did not converge or had an error. Removing thread if from the results..."))
    }else converged <- c(converged,i)
  }
  
  if(length(converged)==0)
  {
    warning("No threads converged. Returning null...")
    return(NULL)
  }
  
  all.threads <- mclapply(converged,FUN = function(i) predict.surv.thread(thread.id = res.object$threads.stats[[i]],data.train = data.train, data.test = data.test))
  threads.weights <- unlist(lapply(converged, function(x) res.object$threads.stats[[x]]$cterm))
  threads.weights <- exp(threads.weights-min(threads.weights))/sum(exp(threads.weights-min(threads.weights)))
  
  average.lin.pred.train <- array(0,dim(data.train)[1])
  average.lin.pred.test <- array(0,dim(data.test)[1])
  for(i in converged)
  {
    average.lin.pred.train <- average.lin.pred.train + all.threads[[i]]$train*threads.weights[i]
    average.lin.pred.test <- average.lin.pred.test + all.threads[[i]]$test*threads.weights[i]
  }
  
  data.train <- cbind(data.train[,1],data.train[,2], average.lin.pred.train)
  colnames(data.train)[1:3] <- c('time', 'cens','average.lin.pred')
  
  data.test <- cbind(data.test[,1], data.test[,2], average.lin.pred.test)
  colnames(data.test)[1:3] <- c('time', 'cens','average.lin.pred')
  
  mod <- coxph(Surv(time, cens) ~ average.lin.pred, data = as.data.frame(data.train),
               init = 1, control = coxph.control(iter.max = 0), x = TRUE)
  pred <- pec(mod, formula = mod$formula, data = as.data.frame(data.test))
  
  
  # compute the prediction measurements
  ibs <- ibs(pred, times = min(max(data.train[,1]), max(data.test[,1])))
  cindex <- cindex(mod, mod$formula, data = as.data.frame(data.test), cens.model = 'cox')$AppCindex
  
  c(as.vector(ibs), cindex)
  
  return(list(curves.pec = pred, ibs = as.vector(ibs),cindex = cindex))
}





formula1 = as.formula(paste("time ~ 1 +",paste0(colnames(data)[-c(1,2)],collapse = "+")))



#define the number or cpus
M = 32
#define the size of the simulated samples
NM= 1000
#define \k_{max} + 1 from the paper
compmax = 16
#define treshold for preinclusion of the tree into the analysis
th=(10)^(-5)
#define a final treshold on the posterior marginal probability for reporting a tree
thf=0.05
#specify tuning parameters of the algorithm for exploring DBRM of interest
#notice that allow_offsprings=3 corresponds to the GMJMCMC runs and
#
g = function(x) x
results=array(0,dim = c(2,100,5))

for(j in 1:100)
{
  
  #specify the initial formula
  set.seed(j)
  
  res1 = pinferunemjmcmc(predict = F,  test.data = NULL, n.cores = 32, report.level =  0.1, num.mod.best = 1000 ,simplify = T,link.function = g,
                         runemjmcmc.params = list(formula = formula1,outgraphs=F,data = data.example,estimator = estimate.bic.surv.cpen,  
                                                  estimator.args =  list(data = data.example),recalc_margin = 249, 
                                                  save.beta = T,interact = T,gen.prob = c(5,0,1,0,0),relations=c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2"),
                                                  relations.prob = rep(0.1,16),
                                                  interact.param=list(deep.method=1,allow_offsprings=3,mutation_rate = 250,last.mutation = 15000, max.tree.size = 2, 
                                                                      Nvars.max = 15,p.allow.replace=0.5,p.allow.tree=0.3,p.nor=0.3,p.and = 0.9),n.models =  20000,unique = F,
                                                  max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,
                                                  advanced.param = list(
                                                    max.N.glob=as.integer(10),
                                                    min.N.glob=as.integer(5),
                                                    max.N=as.integer(3),
                                                    min.N=as.integer(1),
                                                    printable = F)))
  print(res1$feat.stat)
  
  
  preds = predict.surv.bma(res.object = res1,data.train = data.example, data.test = data.test)
  
  preds.first =  predict.surv.bma(res.object = list(threads.stats = list(res1$threads.stats[[2]])),data.train = data.example, data.test = data.test)
  
  results[1,j,1]=  preds.first$ibs[1]
  results[1,j,2]=  preds.first$ibs[2]
  results[1,j,3] = preds.first$cindex$coxph
  
  
  results[2,j,1]=   preds$ibs[1]
  results[2,j,2]=   preds$ibs[2]
  results[2,j,3]=   preds$cindex$coxph  
  
  
  write.csv(x =res1$feat.stat,row.names = F,file = paste0("posteriorsurvFP_",j,".csv"))
  print(paste0("end simulation ",j))
  
  #print the run's metrics and clean the results
  write.csv(file =paste0("resultsurvrun_",j,".csv"),x= results[,j,])
  #rm(results)
  
  
  print(results[1:2,j,3])
  
  gc()
  
}

for(j in 1:100)
{
  print(file.info(paste0("resultsurvrun_",j,".csv"))$ctime)
  results[1:2,j,1:5] = as.matrix(read.csv(paste0("resultsurvrun_",j,".csv"))[,-1])
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

round(summary.results[,c(8,7,9,2,1,3,5,4,6)],4)

#featgmj = hash()

simplifyposteriors<-function(X,post,th=0.0001,thf=0.1,y = "time")
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
  tmpf = read.csv(paste0("posteriorsurvFP_",j,".csv"))
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

tmp = simplifyposteriors(X = data.example,post =as.data.frame(cbind(keys(featgmj),as.numeric(values(featgmj)))),y = "time")

write.csv(x =tmp,row.names = T,file = "survivalfeatFP.csv")
#print(paste0("end simulation ",j))


