\name{LogicRegr}
\alias{LogicRegr}
\title{A wrapper for running the Bayesian logic regression based inference in a easy to use way}
\usage{LogicRegr(formula, data, family = "Gaussian",prior = "J",report.level = 0.5, d = 20, cmax = 5, kmax = 20, p.and = 0.9, p.not = 0.05, p.surv = 0.1,ncores = -1, n.mods = 1000 ,advanced)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{family}{a string taking values of either "Gaussian" or "Bernoulli" correspodning to the linear or logistic Bayesian logic regression contexts}
\item{prior}{character values "J" or "G" corresponing either to Jeffey's or robust g prior}
\item{report.level}{a numeric value in (0,1) specifying the treshold for detections based on the marginal inclusion probabilities}
\item{d}{population size for the GMJMCMC algorithm}
\item{cmax}{the maximal allowed depth of logical expressions to be considered}
\item{kmax}{the maximal number of logical expressions per model}
\item{p.and}{probability of AND parameter of GMJMCMC algorithm}
\item{p.not}{probability of applying logical NOT in GMJMCMC algorithm}
\item{p.surv}{minimal survival probabilities for the features to be allowed to enter the next population}
\item{ncores}{the maximal number of cores (and GMJMCMC threads) to be addressed in the analysis}
\item{n.mods}{the number of the best models in the thread to calculate marginal inclusion probabilities}
\item{advanced}{should only be adrresed by experienced users to tune advanced parameters of GMJMCMC, advanced corresponds to the vector of tuning parameters of runemjmcmc function}
}
\value{a list of
\item{feat.stat}{detected logical expressions and their marginal inclusion probabilities}
\item{predictions}{NULL currently, since LogrRegr function is not designed for predictions at the moment, which is still possible in its expert mother function pinferunemjmcmc}
\item{allposteriors}{all visited by GMJMCMC logical expressions and their marginal inclusion probabilities}
\item{threads.stats}{a vector of detailed outputs of individual ncores threads of GMJMCMC run}
}
\seealso{EMJMCMC::runemjmcmc, EMJMCMC::pinferunemjmcmc}
\examples{

set.seed(040590)
X1= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
Y1=rnorm(n = 1000,mean = 1+0.7*(X1$V1*X1$V4) + 0.8896846*(X1$V8*X1$V11)+1.434573*(X1$V5*X1$V9),sd = 1)
X1$Y1=Y1

#specify the initial formula
formula1 = as.formula(paste(colnames(X1)[51],"~ 1 +",paste0(colnames(X1)[-c(51)],collapse = "+")))
data.example = as.data.frame(X1)


#run the inference with robust g prior
res4G = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4G$feat.stat)
#run the inference with Jeffrey's prior
res4J = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "J",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
print(res4J$feat.stat)


}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
