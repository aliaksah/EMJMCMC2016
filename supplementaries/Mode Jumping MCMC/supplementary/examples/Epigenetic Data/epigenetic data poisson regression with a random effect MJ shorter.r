rm(list = ls(all = TRUE))
# install the required packges if needed
#install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/testing")
#install.packages("bigmemory")
#install.packages("snow")
#install.packages("Rmpi")
#install.packages("ade4")
#install.packages("sp")
#install.packages("BAS")
#install.packages("https://github.com/aliaksah/EMJMCMC2016/files/270429/EMJMCMC_1.2.tar.gz", repos = NULL, type="source")
#install.packages("RCurl")
#install.packages("hash")

library(hash)
library(RCurl)
library(EMJMCMC)
library(sp)
library(INLA)
library(parallel)
library(bigmemory)
library(snow)
library(MASS)
library(ade4)
library(compiler)
library(BAS)
require(stats)

#define the working directory

workdir<-""

# get the data
M<-5
size<-1

data.example <- read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/Epigenetic%20Data/epigen.txt"),sep = ",",header = T)[,2:30]


data.example<-data.example[,c(2,5,6,8:10,12:17,21,23,24,29)]
data.example$eg3000<-data.example$express_noisy>3000
data.example$eg10000<-data.example$express_noisy>10000
#data.example$eg295000<-data.example$express_noisy>295000
data.example$express_noisy<-NULL

fparam.example <- colnames(data.example )[-c(1,2,3,4)]
fobserved.example <- colnames(data.example)[2]


#create MySearch object with default parameters. N/B default estimator is INLA!
mySearch = EMJMCMC2016()
mySearch$parallelize = mclapply

SparseM::image(cor(data.example)[16:1,],axes = FALSE, col = grey(seq(1, 0, length = 256)))

# specify some INLA realted parameters (if the data is analyzed for the very first time)
mySearch$estimator = inla
args<-list(family = "poisson",data = data.example)
args$control.compute = list(dic = TRUE, waic = TRUE, mlik = TRUE)
mySearch$latent.formula  = "+offset(log(total_bases))+f(data.example$pos,model=\"ar1\")"
mySearch$estimator.args = args
mySearch$printable.opt = T

#example of the underlying model within INLA
formula2 <- as.formula(paste0("methylated_bases ~  1+",paste(fparam.example,collapse = "+"),"+offset(log(total_bases))+f(data.example$pos,model=\"ar1\")"))
fm4<-do.call(inla, c(args,formula = formula2))
summary(fm4)

#estimate.inla.ar1(formula = formula2,args)

#end defining the precalculated results

# create a big memory object to store the results

statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)

# carry out full enumeration to learn about the truth. This one MUST be completed before moving to the experiments in this example!
system.time(
  FFF<-mySearch$full_selection(list(statid=6, totalit =2^13+1, ub = 10,mlikcur=-Inf,waiccur =Inf))
)

write.big.matrix(statistics1, "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/submitted/aliaksah-EMJMCMC2016-04a333b/examples/Epigenetic Data/precalculated.csv", row.names = FALSE, col.names = FALSE,
                 sep = ",")

# use the precalculated results to save time (if available). This should not get addressed if the data set is analysed for the first time!
crits<-as.big.matrix(read.table("/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/submitted/aliaksah-EMJMCMC2016-04a333b/examples/Epigenetic Data/precalculated.csv",sep = ",")[,1:3,15])
Nvars<-mySearch$Nvars
bittodec<-mySearch$bittodec
# estimator based on precalculated and saved into crit data
esimator.inla.prec.constrained<-function(formula, crits)
{
  values <- strsplit(as.character(formula)[3],split = " + ",fixed = T)
  vec<-array(0,dim = Nvars)
  for(i in 2:(length(values[[1]])))
  {
    iid <- which(fparam.example == values[[1]][i])
    if(length(iid)>0)
      vec[iid]<-1
  }

  # incorporate constraints
  if(sum(vec[10:12])==3)
    return(list(mlik = -100000,waic = 100000 , dic = 100000))
  id<-bittodec(vec)+1
  return(list(mlik = crits[id,1]+1000,waic = crits[id,2] , dic =  crits[id,3]))
}

mySearch$estimator = esimator.inla.prec.constrained
mySearch$estimator.args = list(crits = crits)


statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
statistics <- describe(statistics1)
# try the estimator function based on precalculated values out
esimator.inla.prec.constrained(formula = formula2, crits = crits)

# carry out full enumeration to learn about the truth. This one MUST be completed before moving to the experiments in this example!
system.time(
  FFF<-mySearch$full_selection(list(statid=6, totalit =2^13+1, ub = 10,mlikcur=-Inf,waiccur =Inf))
)


# check that all models are enumerated during the full search procedure
idn<-which(!is.na(statistics1[,1]))
length(idn)

# draw the model space and get other graphical output
mySearch$visualize_results(statistics1, "test",1024, crit=list(mlik = T, waic = T, dic = T),draw_dist = TRUE)


# once full search is completed, get the truth for the experiment
ppp<-mySearch$post_proceed_results(statistics1 = statistics1)
truth = ppp$p.post # make sure it is equal to Truth column from the article
truth.m = ppp$m.post
truth.prob = ppp$s.mass
ordering = sort(ppp$p.post,index.return=T)
fake500 <- sum(exp(x = (sort(statistics1[,1],decreasing = T)[1:2])),na.rm = TRUE)/truth.prob
print("pi truth")
sprintf("%.10f",truth[ordering$ix])
sprintf(mySearch$fparam[ordering$ix])


# get the best performance results for a given number of iterations
iddx <- sort(statistics1[,1],decreasing = T,index.return=T,na.last = NA)$ix
statistics1[as.numeric(iddx[386:2^13]),1:15]<-NA

ppp.best<-mySearch$post_proceed_results(statistics1 = statistics1)
best = ppp.best$p.post # make sure it is equal to Truth column from the article
bset.m = ppp.best$m.post
best.prob = ppp.best$s.mass/truth.prob
print("pi best")
sprintf("%.10f",best[ordering$ix])

best.bias.m<-sqrt(mean((bset.m - truth.m)^2,na.rm = TRUE))*100
best.rmse.m<-sqrt(mean((bset.m - truth.m)^2,na.rm = TRUE))*100

best.bias<- best - truth
best.rmse<- abs(best - truth)

# view the "unbeatible" results
View((cbind(best.bias[ordering$ix],best.rmse[ordering$ix])*100))

# define parameters of the search

#mySearch$printable.opt=T
mySearch$printable.opt = F
mySearch$parallelize = lapply
#mySearch$printable.opt = TRUE
mySearch$max.cpu = as.integer(4)
mySearch$locstop.nd=FALSE
mySearch$max.cpu.glob = as.integer(4)
mySearch$max.N.glob=as.integer(6)
mySearch$min.N.glob=as.integer(4)
mySearch$max.N=as.integer(1)
mySearch$min.N=as.integer(1)
mySearch$recalc.margin = as.integer(2^13)
mySearch$max.cpu.hyper=as.integer(1)

distrib_of_proposals = c(10,10,10,10,1000)
#distrib_of_proposals = c(0,0,0,0,1000)

distrib_of_neighbourhoods=t(array(data = c(7.6651604,16.773326,14.541629,12.839445,2.964227,13.048343,7.165434,
                                           0.9936905,15.942490,11.040131,3.200394,15.349051,5.466632,14.676458,
                                           1.5184551,9.285762,6.125034,3.627547,13.343413,2.923767,15.318774,
                                           14.5295380,1.521960,11.804457,5.070282,6.934380,10.578945,12.455602,
                                           1,1,1,1,1,1,1),dim = c(7,5)))
distrib_of_neighbourhoods[7]=distrib_of_neighbourhoods[7]/50
# carry out the experiment (notice that result may vary depending on the part of genome addressed)

Niter <- 100
thining<-1
system.time({

  vect <-array(data = 0,dim = c(length(fparam.example),Niter))
  vect.mc <-array(data = 0,dim = c(length(fparam.example),Niter))
  inits <-array(data = 0,dim = Niter)
  freqs <-array(data = 100,dim = c(5,Niter))
  freqs.p <-array(data = 100,dim = c(5,7,Niter))
  masses <- array(data = 0,dim = Niter)
  biases.m <- array(data = 0,dim = 2 ^(length(fparam.example))+1)
  biases.m.mc <- array(data = 0,dim = 2 ^(length(fparam.example))+1)
  rmse.m <- array(data = 0,dim = Niter)
  rmse.m.mc <- array(data = 0,dim = Niter)
  iterats <- array(data = 0,dim = c(2,Niter))

  for(i in 1:Niter)
  {
    statistics1 <- big.matrix(nrow = 2 ^(length(fparam.example))+1, ncol = 15,init = NA, type = "double")
    statistics <- describe(statistics1)
    mySearch$g.results[4,1]<-0
    mySearch$g.results[4,2]<-0
    mySearch$p.add = array(data = 0.5,dim = length(fparam.example))
    print("BEGIN ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(i)
    set.seed(10*i)
    initsol=rbinom(n = length(fparam.example),size = 1,prob = 0.5)
    inits[i] <- mySearch$bittodec(initsol)
    freqs[,i]<- distrib_of_proposals
    resm<-mySearch$modejumping_mcmc(list(varcur=initsol,statid=-1, distrib_of_proposals =distrib_of_proposals,distrib_of_neighbourhoods=distrib_of_neighbourhoods, eps = 0.000000001, trit = 3150, trest = 375000, burnin = 3, max.time = 30, maxit = 20000000, print.freq = 500))
    vect[,i]<-resm$bayes.results$p.post
    vect.mc[,i]<-resm$p.post
    masses[i]<-resm$bayes.results$s.mass/truth.prob
    print(masses[i])
    freqs.p[,,i] <- distrib_of_neighbourhoods
    cur.p.post <- resm$bayes.results$m.post
    cur.p.post[(which(is.na(cur.p.post)))]<-0
    rmse.m[i]<-mean((cur.p.post - truth.m)^2,na.rm = TRUE)
    biases.m<-biases.m + (cur.p.post - truth.m)
    cur.p.post.mc <- resm$m.post
    cur.p.post.mc[(which(is.na(cur.p.post.mc)))]<-0
    rmse.m.mc[i]<-mean((cur.p.post.mc - truth.m)^2,na.rm = TRUE)
    biases.m.mc<-biases.m.mc + (cur.p.post.mc - truth.m)
    iterats[1,i]<-mySearch$g.results[4,1]
    iterats[2,i]<-mySearch$g.results[4,2]
    print("COMPLETE ITERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! with")
    print(iterats[2,i])
    if(i==1)
    {
      # draw the model space and get other graphical output
      mySearch$visualize_results(statistics1, "MJMCMC",1024, crit=list(mlik = T, waic = T, dic = T),draw_dist = TRUE)
      
    }
    remove(statistics1)
    remove(statistics)

  }
}
)


Nlim <- 1
order.deviat <- sort(masses,decreasing = TRUE,index.return=T)


print("model bias rm")
sqrt(mean((biases.m/Niter)^2,na.rm = TRUE))*100000
print("model rmse rm")
sqrt(mean(rmse.m))*100000

print("model bias mc")
sqrt(mean((biases.m.mc/Niter)^2,na.rm = TRUE))*100000
print("model rmse mc")
sqrt(mean(rmse.m.mc))*100000


print("model coverages")
mean(masses)
median(masses)
print("mean # of iterations")# even smaller on average than in BAS
mean(iterats[1,])
print("mean # of estimations")# even smaller on average than in BAS
mean(iterats[2,])

hist(masses)


# correlation between the MSE and the masses, obviously almost minus 1
cor(rmse.m,masses)
cor(rmse.m.mc,masses)
cor(iterats[2,],masses)

truth.buf <- array(data = 0,dim = c(length(fparam.example),Niter))
truth.buf[,1:Niter]<-truth
bias <- vect - truth.buf
bias.mc <- vect.mc - truth.buf
rmse <- (vect^2 +truth.buf^2 - 2*vect*truth.buf)
rmse.mc <- (vect.mc^2 +truth.buf^2 - 2*vect.mc*truth.buf)
bias.avg.rm<-rowMeans(bias)
rmse.avg.rm <-sqrt(rowMeans(rmse))
bias.avg.mc<-rowMeans(bias.mc)
rmse.avg.mc <-sqrt(rowMeans(rmse.mc))


print("pi biases rm")
sprintf("%.10f",bias.avg.rm[ordering$ix]*100)
print("pi rmse rm")
sprintf("%.10f",rmse.avg.rm[ordering$ix]*100)

print("pi biases mc")
sprintf("%.10f",bias.avg.mc[ordering$ix]*100)
print("pi rmse mc")
sprintf("%.10f",rmse.avg.mc[ordering$ix]*100)

# View the results
View((cbind(ordering$ix/100,truth[ordering$ix]/100,bias.avg.rm[ordering$ix],rmse.avg.rm[ordering$ix],bias.avg.mc[ordering$ix],rmse.avg.mc[ordering$ix])*100))


