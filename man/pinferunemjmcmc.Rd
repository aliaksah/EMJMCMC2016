\name{pinferunemjmcmc}
\alias{pinferunemjmcmc}
\title{A wrapper for running the GLMM, BLR, or DBRM based inference and predictions in an expert but rather easy to use way}
\usage{pinferunemjmcmc(n.cores = 4, mcgmj = mcgmjpar, report.level =  0.5, simplify = F, num.mod.best = 1000, predict = F,  test.data = 1, link.function = function(z)z, runemjmcmc.params)}
\arguments{
\item{runemjmcmc.params}{a vector of parameters of runemjmcmc function, see the help of runemjmcmc for details}
\item{mcgmj}{an mclapply like function for perfroming for perfroming parallel computing, do not change the default unless you are using Windows}
\item{simplify}{a logical value specifying in simplification of the features is to be done after the search is completed}
\item{predict}{a logical value specifying if predictions should be done by the run of pinferunemjmcmc}
\item{test.data}{covariates data.frame to be used for predictions}
\item{link.function}{the link functions to be used to make predictions}
\item{report.level}{a numeric value in (0,1) specifying the treshold for detections based on the marginal inclusion probabilities}
\item{n.cores}{the maximal number of cores (and (R)(G)MJMCMC threads) to be addressed in the analysis}
\item{num.mod.best}{the number of the best models in the thread to calculate marginal inclusion probabilities}
}
\value{a list of
\item{feat.stat}{detected features or logical expressions and their marginal inclusion probabilities}
\item{predictions}{predicted values if they are required, NULL otherwise}
\item{allposteriors}{all visited by (R)(G)MJMCMC features and logical expressions and their marginal inclusion probabilities}
\item{threads.stats}{a vector of detailed outputs of individual n.cores threads of (R)(G)MJMCMC run}
}
\seealso{EMJMCMC::runemjmcmc, EMJMCMC::LogrRegr, EMJMCMC::DeepRegr, EMJMCMC::LinRegr}
\examples{

#inference

X=read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/DBRM\%20supplementaries/kepler\%20and\%20mass/exa1.csv")
data.example = as.data.frame(X)

#specify the initial formula
formula1 = as.formula(paste(colnames(X)[5],"~ 1 +",paste0(colnames(X)[-5],collapse = "+")))

#a set of nonlinearities that will be used in the DBRM model
sini=function(x)sin(x/180*pi)
expi=function(x)exp(-abs(x))
logi =function(x)log(abs(x)+1)
troot=function(x)abs(x)^(1/3)
to23=function(x)abs(x)^(2.3)
to35=function(x)abs(x)^(3.5)


#specify the estimator function returning p(Y|m)p(m), model selection criteria and the vector of the modes for the beta coefficients
estimate.gamma.cpen = function(formula, data,r = 1.0/223.0,logn=log(223.0),relat=c("to23","expi","logi","to35","sini","troot","sigmoid"))
{
  fparam=NULL
  fmla.proc=as.character(formula)[2:3]
  fobserved = fmla.proc[1]
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
  fmla.proc[2]=stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
  fparam =stri_split_fixed(str = fmla.proc[2],pattern = "+I",omit_empty = F)[[1]]
  sj=(stri_count_fixed(str = fparam, pattern = "*"))
  sj=sj+(stri_count_fixed(str = fparam, pattern = "+"))
  for(rel in relat)
    sj=sj+(stri_count_fixed(str = fparam, pattern = rel))
  sj=sj+1
  tryCatch(capture.output({
    out = glm(formula = formula,data = data, family = gaussian)
    mlik = (-(out$deviance -2*log(r)*sum(sj)))/2
    waic = -(out$deviance + 2*out$rank)
    dic =  -(out$deviance + logn*out$rank)
    summary.fixed =list(mean = coefficients(out))

  }, error = function(err) {
    print(err)
    mlik = -10000
    waic = -10000
    dic =  -10000
    summary.fixed =list(mean = array(0,dim=length(fparam)))
  }))
  return(list(mlik = mlik,waic = waic , dic = dic,summary.fixed =summary.fixed))

}
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
#allow_offsprings=4 -to the RGMJMCMC runs
res1 = pinferunemjmcmc(n.cores = M, report.level =  0.5, num.mod.best = NM ,simplify = T,runemjmcmc.params = list(formula = formula1,data = data.example,estimator = estimate.gamma.cpen,estimator.args =  list(data = data.example),recalc_margin = 249, save.beta = F,interact = T,outgraphs=F,relations=c("to23","expi","logi","to35","sini","troot","sigmoid"),relations.prob =c(0.1,0.1,0.1,0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=3,mutation_rate = 250,last.mutation=10000, max.tree.size = 5, Nvars.max =15,p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.9,p.and = 0.9),n.models = 10000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F)))

print(res1$feat.stat)


#prediction

#load some of the possible nonlinearities of interest
troot=function(x)abs(x)^(1/3)
sini=function(x)sin(x/180*pi)
logi=function(x)log(abs(x+0.1))
gfquar=function(x)as.integer(x<quantile(x,probs = 0.25))
glquar=function(x)as.integer(x>quantile(x,probs = 0.75))
gmedi=function(x)as.integer(x>median(x))
cosi=function(x)cos(x/180*pi)
gmean=function(x)as.integer(x>mean(x))
gone=function(x)as.integer(x>0)
gthird=function(x)(abs(x)^(1/3))
gfifth=function(x)(abs(x)^(1/5))
grelu=function(x)(x*(x>0))
contrelu=function(x)log(1+exp(x))
gauss=function(x)exp(-x*x)
compmax = 21

#read in the train and test data sets
test = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/DBRM\%20supplementaries/breast\%20cancer/test.csv",header = T,sep=",")[,-1]
train = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/DBRM\%20supplementaries/breast\%20cancer/train.csv",header = T,sep=",")[,-1]

#transform the train data set to a data.example data.frame that EMJMCMC class will internally use
data.example = as.data.frame(train)

#specify the link function that will be used in the prediction phase
g=function(x)
{
  return((x = 1/(1+exp(-x))))
}

formula1 = as.formula(paste(colnames(data.example)[31],"~ 1 +",paste0(colnames(data.example)[-31],collapse = "+")))


res = pinferunemjmcmc(n.cores =30, report.level =  0.5 , num.mod.best = NM,simplify = T, predict = T,test.data = as.data.frame(test),link.function = g, runemjmcmc.params =list(formula = formula1,data = data.example,gen.prob = c(1,1,1,1,0),estimator =estimate.bas.glm.cpen,estimator.args =  list(data = data.example,prior = aic.prior(),family = binomial(),yid=31, logn = log(143),r=exp(-0.5)),recalc_margin = 95, save.beta = T,interact = T,relations = c("gauss","tanh","atan","sin"),relations.prob =c(0.1,0.1,0.1,0.1),interact.param=list(allow_offsprings=4,mutation_rate = 100,last.mutation=1000, max.tree.size = 6, Nvars.max = 20,p.allow.replace=0.5,p.allow.tree=0.4,p.nor=0.3,p.and = 0.9),n.models = 7000,unique =T,max.cpu = 4,max.cpu.glob = 4,create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,print.freq = 1000,advanced.param = list(
  max.N.glob=as.integer(10),
  min.N.glob=as.integer(5),
  max.N=as.integer(3),
  min.N=as.integer(1),
  printable = F)))

print(auc(prob = res$predictions,y = test$X))
for(i in 1:M){
  print(auc(prob = res$threads.stats[[i]]$preds,y = test$X))
  print( res$threads.stats[[i]]$post.populi)
  }

for(jjjj in 1:10)
{
  resw = as.integer(res$predictions>=0.1*jjjj)
  prec=(1-sum(abs(resw-test$X),na.rm = T)/length(resw))
  print(prec)
  #FNR
  ps=which(test$X==1)
  fnr=sum(abs(resw[ps]-test$X[ps]))/(sum(abs(resw[ps]-test$X[ps]))+length(ps))

  #FPR
  ns=which(test$X==0)
  fpr=sum(abs(resw[ns]-test$X[ns]))/(sum(abs(resw[ns]-test$X[ns]))+length(ns))

}

}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
