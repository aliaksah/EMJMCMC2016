#library(dplyr)
library(plyr)
library(dplyr)
library(magrittr)
library(rootSolve)

#setwd("/home/michaelh/SIMULATION_paper/")

#options("exppression" = 500000)
##########################################
#: Calculate True h^2 for each Scenario :#
##########################################
simPars <- 
  read.table("scripts/simulationParameters_new.txt", 
             header = TRUE, 
             stringsAsFactors = FALSE) %>%
  select(SNPId = SNP_causal, Pos = BP_causal, MAF = MAF_causal, S1, S2, S3, S4)

causSNPsS1 <- with(simPars, SNPId[S1 != 0])
causSNPsS2 <- with(simPars, SNPId[S2 != 0]) 
causSNPsS3 <- with(simPars, SNPId[S3 != 0]) 
causSNPsS4 <- with(simPars, SNPId[S4 != 0]) 

genoData <- read.table("data/CHR1_NFBC.raw",
                       header = TRUE,
                       stringsAsFactors = FALSE)
names(genoData) <- gsub("_.$", "", names(genoData)) %>% gsub("\\.", "-", .)

genoData <- genoData[c("FID", "IID", simPars$SNPId)]

expectedData <- data.frame(as.matrix(genoData[-(1:2)]) %*% as.matrix(simPars[c('S1', 'S2', 'S3', 'S4')])) 

explainedVar <- apply(expectedData, 2, function(x) var(x)/(var(x)+1)) # True h^2


############################################
#: Calculate Neigbourhoods of causal SNPs :#
#: Definition: Neigbourhood 		  :#
#	       corr > .3 and +/- 1(2) MBP :#
############################################
distance_cutoff <- 1e6
corr_cutoff <- .3

physMap <- read.table("data/CHR1_NFBC.bim",
                      header = FALSE,
                      stringsAsFactors = FALSE)[c(2, 4)] %>%
  select(SNPid = V2, pos = V4) %>%
  filter(!grepl('^cnv', SNPid))	   

genoData <- read.table("data/CHR1_NFBC.raw", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
names(genoData) <- gsub("_.$", "", names(genoData)) %>% gsub("\\.", "-", .)	 

findNeigbour <- function(causSNP) {
  #causSNP <- "rs1498308"
  #print(causSNP)
  physMap %>%
    mutate(distance = abs(pos - pos[SNPid == causSNP])) %>%
    filter(distance <= distance_cutoff) %>%
    mutate(corr = sapply(SNPid, function(x) abs(cor(genoData[causSNP], genoData[x]))),
           causSNPid = causSNP) %>%
    filter(corr >= corr_cutoff)
}

#dataNeigbourhoodS1 <- plyr::ldply(causSNPsS1, findNeigbour)
dataNeigbourhoodS2 <- plyr::ldply(causSNPsS2, findNeigbour)
#dataNeigbourhoodS3 <- plyr::ldply(causSNPsS3, findNeigbour)
#dataNeigbourhoodS4 <- plyr::ldply(causSNPsS4, findNeigbour)

rm(genoData)
rm(expectedData)
rm(physMap)
rm(simPars)
rm(causSNPsS4)
rm(causSNPsS3)
rm(causSNPsS1)
rm(findNeigbour)
gc()

source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")


pheno<-read.csv(paste0("data_S2_nocausal_5402/pimass/data.recode.pheno_",1,".txt"),header = F)
geno<-t(read.csv("data_S2_nocausal_5402/pimass/data.recode.mean.geno.txt",header = F,stringsAsFactors = F))
names<-geno[1,]
geno<-as.data.frame(geno[-c(1,2,3),])
names(geno)<-names
geno$Y<-pheno$V1
geno <- as.data.frame(mclapply(geno, as.numeric))




estimate.lm.MAIC2 <- function(formula, data, n = 5402, m = 24592, c = 4,u=170)
{
  size<-stri_count_fixed(str = as.character(formula)[3],pattern = "+")
  
  if(size>u)
  {
    return(list(mlik = (-50000 + rnorm(1,0,1) - size*log(m*m*n/c) + 2*log(factorial(size+1))),waic =0, dic =  0,summary.fixed =list(mean = array(0,size+1))))
  }else{
    out <- lm(formula = formula,data = data,x = T)
    logmarglik <- (2*logLik(out) - 2*out$rank - 2*out$rank*log(m/c-1) + 2*log(factorial(out$rank)))/2
    #sss<-summary(out)
    #her  = summary(out)$r.squared
    #ser  = ser<-2*her*(1-her*her)*(n-out$rank-1)/sqrt((n*n-1)*(3+n))
    sss<-summary(out)
    if(length(which(is.na(out$coefficients))>0))
      cfs<-out$coefficients[-which(is.na(out$coefficients))]#[-1]
    else
      cfs<-out$coefficients#[-1]
    
    covs<-cov(out$x)
    covb<-sss$cov.unscaled
    #print(dim(covs)[1])
    #print(length(cfs))
    
    if(length(cfs)==dim(covs)[1])
    {
      her<-(t(cfs)%*%covs%*%(cfs)*n/(sss$sigma^2*sss$df[2]))[1,1]
      #ser<-sqrt(4*t(cfs)%*%covs%*%covb%*%(covs)%*%cfs*n*n/(sss$sigma^4*(sss$df[2]^2))) #the Gosha's version 2*her*(1-her*her)*(n-out$rank-1)/sqrt((n*n-1)*(3+n))
      ser<-sqrt(4*her/((n-out$rank-1)))
      #ser3<-2*her*(1-her*her)*(n-out$rank-1)/sqrt((n*n-1)*(3+n))
    }
    if(is.na(her))
    {
      her  = 0
      ser  = 0
    }
    
    if(her==0)
    {
      return(list(mlik = (-50000 + rnorm(1,0,1) - size*log(m*m*n/c) + 2*log(factorial(size+1))),waic =0, dic = 0,summary.fixed =list(mean = array(0,size+1))))
    }
  
    return(list(mlik = logmarglik,waic = her , dic =  ser,summary.fixed =list(mean = coef(out))))
  }
}

do.call.emjmcmc<-function(vect)
{
  
  set.seed(as.integer(vect$cpu))
  do.call(runemjmcmc, vect[1:vect$simlen])
  vals<-values(hashStat)
  fparam<-mySearch$fparam
  cterm<-max(vals[1,],na.rm = T)
  ckey<-keys(hashStat)[which(vals[1,]==cterm)]
  ckey<-fparam[stri_locate_all_fixed(str = ckey,pattern = "1")[[1]][,1]]
  ppp<-mySearch$post_proceed_results_hash(hashStat = hashStat)
  post.populi<-sum(exp(values(hashStat)[1,][1:vect$NM]-cterm),na.rm = T)
  herac = sum(ppp$m.post*(values(hashStat)[2,]),na.rm = T)
  hers  = values(hashStat)[2,]
  sers = values(hashStat)[3,]
  f<-function(xu)
  {
    sum(pnorm(q = xu,mean = hers, sd = sers)*ppp$m.post)-1+0.025
  }
  hu<-uniroot(f = f,interval = c(-1000,1000), tol = 0.0001,extendInt="yes")$root
  rm(f)
  f<-function(xu)
  {
    sum(pnorm(q = xu,mean = hers, sd = sers)*ppp$m.post)-0.025
  }
  hl<-uniroot(f = f,interval = c(-1000,1000), tol = 0.0001,extendInt="yes")$root
  clear(hashStat)
  rm(hashStat)
  rm(vals)
  rm(f)
  gc()
  return(list(post.populi = post.populi, p.post =  ppp$p.post, cterm = cterm, fparam = fparam, best = ckey, herac = herac,CI = c(hl,hu)))
}

MM = 10
M = 31
size.init=1000
NM= 1000

compmax = 101
th<-(10)^(-5)
thf<-0.05

data.example<-geno

rm(geno)
gc()

simplifyposteriors<-function(posteriors,th=0.0001,thf=0.2, dataNeigbourhood = dataNeigbourhoodS2)
{
  tds<-which(posteriors[,2]<th)
  if(length(tds)>0)
    posteriors<-posteriors[-tds,]
  posteriors[,1]<-stri_replace(str =  posteriors[,1],fixed = "I(",replacement = "")
  posteriors[,1]<-stri_replace(str =  posteriors[,1],fixed = ")",replacement = "")
  rhash<-hash()
  for(i in 1:length(posteriors[,1]))
  {
    expr<-posteriors[i,1]
    #print(expr)
    id<-which(dataNeigbourhood$SNPid==expr)
    if(length(id)==0)
    {
      key<-expr
    }else
    {
      key<-dataNeigbourhood$causSNPid[id]
    }
    if(!((key %in% keys(rhash))))
      rhash[[key]]<-c(key,posteriors[i,2])
    else
    {
      if(key %in% keys(rhash))
      {
        rhash[[key]][2]<- (as.numeric(rhash[[key]][2]) + posteriors[i,2])
      }
    }
    
  }
  res<-as.data.frame(t(values(rhash)[c(2,1),]))
  res$V1<-as.numeric(as.character(res$V1))
  res<-res[which(res$V1>thf),]
  res<-res[order(res$V1, decreasing = T),]
  clear(rhash)
  rm(rhash)
  res[which(res[,1]>1),1]<-1
  colnames(res)<-c("posterior","tree")
  return(res)
}

j=1 

for(j in 31:100)
{
  tryCatch({
    
    set.seed(j)
    
    pheno<-read.csv(paste0("data_S2_nocausal_5402/pimass/data.recode.pheno_",j,".txt"),header = F)
    data.example$Y<-as.numeric(pheno$V1)
    rm(pheno)
    cors<-cor(data.example$Y,data.example[,1:24602])
    gc()
    cov.names<-names[which(abs(cors)>0.05)]
    names1<-names[which(abs(cors)>0.025)]
    print(length(names1))
    sum<-summary(lm(as.formula(paste0("Y~1+",paste(cov.names,collapse = "+"))),data = data.example))
    cov.names<-names(sum$coefficients[-1,4])
    detected<-cov.names
    detected<-stri_replace(str = detected,fixed = "I(",replacement = "")
    detected<-stri_replace(str = detected,fixed = ")",replacement = "")
    detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
    detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)
    detlen<-length(detect.true.unique)
    totlen<-length(detected)-length(detect.true)+length(detect.true.unique)
    print(detlen/20)
    print((totlen-detlen)/totlen)
    gc()
    formula1 <- as.formula(paste0("Y~1+",paste(cov.names,collapse = "+")))
    secondary <-names1[-which(names1 %in% cov.names)]
    
    # her = 0
    #for(i in 1:1000)
    #{
     #est = estimate.lm.MAIC2(data = data.example,formula = as.formula(paste0("Y~1+",paste(cov.names[sample.int(n = length(cov.names),size = runif(1,1,length(cov.names)))],collapse = "+"))))
     #print(est$dic)
     #print(c(est$waic-1.96*est$dic,est$waic,est$waic+1.96*est$dic))
    #}
    # print(her/1000)

    vect<-list(formula = formula1, locstop.nd = T, keep.origin = F,p.add = 0.1,max.time = 120, p.add.default = 0.1, pool.cor.prob = T,secondary <-names1[-which(names1 %in% cov.names)], outgraphs=F,data = data.example,estimator = estimate.lm.MAIC2,presearch=F, locstop =T,estimator.args =  list(data = data.example),recalc_margin = 999,gen.prob = c(1,0,0,0,0), save.beta = F,interact = T,relations=c("cos"),relations.prob =c(0.1),interact.param=list(allow_offsprings=3,mutation_rate = 1000, last.mutation = 15000, max.tree.size = 4, Nvars.max =(compmax-1),p.allow.replace=0.7,p.allow.tree=0.25,p.nor=0,p.and = 0.9),n.models = 25000,unique = T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 50,print.freq = 1000,advanced.param = list(
      max.N.glob=as.integer(130),
      min.N.glob=as.integer(10),
      max.N=as.integer(5),
      min.N=as.integer(1),
      printable = F))
    
    
    params <- list(vect)[rep(1,M)]
    
    for(i in 1:M)
    {
      #cov.names<-names[sample.int(n = length(names),size = size.init,prob = abs(cors))]
      #sum<-summary(lm(as.formula(paste0("Y~1+",paste(cov.names,collapse = "+"))),data = geno))
      #cov.names<-names(sum$coefficients[-1,4])
      #params[[i]]$
      params[[i]]$formula = formula1
      params[[i]]$cpu<-i*j
      params[[i]]$simul<-"scenario_JM_"
      params[[i]]$simid<-j
      params[[i]]$NM<-NM
      params[[i]]$simlen<-31
    }
    
    gc()
    
    #res<-do.call(runemjmcmc,args = params[[3]][1:27])
    #res$p.post
    #length(which(!is.na(res$m.post)))
    #detected<-mySearch$fparam[which(res$p.post>0.1)]
    
    #res =do.call(runemjmcmc, params[[3]][1:params[[3]]$simlen])
    
    gc()
    print(paste0("begin simulation ",j))
    results<-parall.gmj(X = params, M = M)
    #print(results)

    #wait()
    
    resa<-array(data = 0,dim = c(compmax,M*3))
    post.popul <- array(0,M)
    max.popul <- array(0,M)
    nulls<-NULL
    
    not.null<-1
    for(k in 1:M)
    {
      if(is.character(results[[k]]))
      {
        nulls<-c(nulls,k)
        next
      }
      if(length(results[[k]])==0)
      {
        nulls<-c(nulls,k)
        next
      }
      else
      {
        not.null <- k
      }
      
    }
    
    
    for(k in 1:M)
    {
      if(k %in% nulls)
      {
        results[[k]]<-results[[not.null]]
      }
      max.popul[k]<-results[[k]]$cterm
      post.popul[k]<-results[[k]]$post.populi
      if(length(resa[,k*3-2])==(length(results[[k]]$fparam)+1))
      {
        resa[,k*3-2]<-c(results[[k]]$fparam,"Post.Gen.Max")
        resa[,k*3-1]<-c(results[[k]]$p.post,results[[k]]$cterm)
        resa[,k*3]<-rep(post.popul[k],length(results[[k]]$p.post)+1)
      }else
      {
        #idsx<-order(results[[k]]$p.post,decreasing = T,na.last = T)
        resa[,k*3-2]<-rep(results[[k]]$fparam[1],length(resa[,k*3-2]))
        resa[,k*3-1]<-rep(0,length(resa[,k*3-1]))
        resa[,k*3]<-rep(-10^9,length(resa[,k*3]))
        max.popul[k]<- -10^9
        post.popul[k]<- -10^9
      }
      
    }
    
    
  
    
    ml.max<-max(max.popul)
    post.popul<-post.popul*exp(-ml.max+max.popul)
    p.gen.post<-post.popul/sum(post.popul)
    hfinal<-hash()
    for(ii in 1:M)
    {
      resa[,ii*3]<-p.gen.post[ii]*as.numeric(resa[,ii*3-1])
      resa[length(resa[,ii*3]),ii*3]<-p.gen.post[ii]
      if(p.gen.post[ii]>0)
      {
        for(jj in 1:(length(resa[,ii*3])-1))
        {
          if(resa[jj,ii*3]>0)
          {
            #print(paste0(ii,"  and ",jj))
            if(as.integer(has.key(hash = hfinal,key =resa[jj,ii*3-2]))==0)
              hfinal[[resa[jj,ii*3-2]]]<-as.numeric(resa[jj,ii*3])
            else
              hfinal[[resa[jj,ii*3-2]]]<-hfinal[[resa[jj,ii*3-2]]]+as.numeric(resa[jj,ii*3])
          }
          
        }
      }
    }
    
    her<-0
    heru<-0
    herl<-0
    for(i in 1:M)
    {
      her<-her+results[[i]]$herac*p.gen.post[[i]]
      heru<-heru+results[[i]]$CI[2]*p.gen.post[[i]]
      herl<-herl+results[[i]]$CI[1]*p.gen.post[[i]]
      #print(results[[i]]$herac)
    }
    print(c(herl,her,heru))
    
    write.csv(x =c(herl,her,heru),row.names = F,file = paste0("herestcmjmcaic2_",j,".csv"))
    
    KMK<-which(max.popul==max(max.popul,na.rm = T))
    
    write.csv(x =c(results[[KMK]]$cterm,results[[KMK]]$best),row.names = F,file = paste0("bestmodcmjmcaic2_",j,".csv"))
    
    gc()
    rm(results)
    
    posteriors<-values(hfinal)
    
    #print(posteriors)
    clear(hfinal)
    rm(hfinal)
    rm(resa)
    rm(post.popul)
    rm(max.popul)
    posteriors<-as.data.frame(posteriors)
    posteriors<-data.frame(X=row.names(posteriors),x=posteriors$posteriors)
    posteriors$X<-as.character(posteriors$X)
    tryCatch({
      res1<-simplifyposteriors(posteriors = posteriors, th,thf)
      row.names(res1)<-1:dim(res1)[1]
      
      detected<-res1$tree
      detected<-stri_replace(str = detected,fixed = "I(",replacement = "")
      detected<-stri_replace(str = detected,fixed = ")",replacement = "")
      
      detect.true.unique<-unique(dataNeigbourhoodS2$causSNPid[which(dataNeigbourhoodS2$SNPid %in% detected)])
      detect.true<-which(detected %in% dataNeigbourhoodS2$SNPid)
      
      detlen<-length(detect.true.unique)
      totlen<-length(detected)-length(detect.true)+length(detect.true.unique)
      
      print(detlen/20)
      print((totlen-detlen)/totlen)
      
      
      
      write.csv(x =res1,row.names = F,file = paste0("post1AGMJSIM_",j,".csv"))
    },error = function(err){
      print("error")
      write.csv(x =posteriors,row.names = F,file = paste0("post1EGMJSIM_",j,".csv"))
    },finally = {
      
      print(paste0("end simulation ",j))
      
    })
    rm(X)
    #rm(data.example)
    rm(vect)
    rm(params)
    gc()
    print(paste0("end simulation ",j))
  },error = function(err){
    print("error")
    j=j-1
    print(paste0("repeat  simulation ",j))
  },finally = {
    
    print(paste0("end simulation ",j))
    rm(vect)
    rm(params)
    gc()
  })
  
}




