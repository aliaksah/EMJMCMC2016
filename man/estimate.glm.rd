\name{estimate.glm}
\alias{estimate.glm}
\title{Obtaining Bayesian estimators of interest from a GLM model}
\usage{estimate.glm(formula, data, family, prior, n=1, g = 0)}
\arguments{
\item{formula}{a formula object for the model to be addressed}
\item{data}{a data frame object containing variables and observations corresponding to the formula used}
\item{family}{distribution family for the responces}
\item{prior}{integers 1,2 or 3 corresonding to AIC, BIC or Gelman's g-prior}
\item{n}{sample size}
\item{g}{g parameter of Gelman's g prior}
}
\value{a list of
\item{mlik}{marginal likelihood of the model}
\item{waic}{AIC model selection criterion}
\item{dic}{BIC model selection criterion}
\item{summary.fixed$mean}{a vector of posterior modes of the parameters}
}
\seealso{glm}
\examples{

X4= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = runif(n = 50*1000,0,1)),dim = c(1000,50)))
Y4=rnorm(n = 1000,mean = 1+7*(X4$V4*X4$V17*X4$V30*X4$V10)+7*(((X4$V50*X4$V19*X4$V13*X4$V11)>0)) + 9*(X4$V37*X4$V20*X4$V12)+ 7*(X4$V1*X4$V27*X4$V3)
             +3.5*(X4$V9*X4$V2) + 6.6*(X4$V21*X4$V18) + 1.5*X4$V7 + 1.5*X4$V8,sd = 1)
X4$Y4=Y4
data.example = as.data.frame(X4)
formula1 = as.formula(paste(colnames(X4)[51],"~ 1 +",paste0(colnames(X4)[-c(51)],collapse = "+")))

formula1 = as.formula(paste(colnames(data.example)[1],"~ 1 +",paste0(colnames(data.example)[-1],collapse = "+")))
EMJMCMC:::estimate.glm(formula = formula1, data = data.example,prior = 2, family = gaussian())
}
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{models}% __ONLY ONE__ keyword per line
