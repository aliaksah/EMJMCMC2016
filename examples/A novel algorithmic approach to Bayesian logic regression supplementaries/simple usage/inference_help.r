#read the most recent stable version of the package
source("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/R/the_mode_jumping_package4.r")

#make sure that you are using Mac Os or Linux (mclapply is currently not supported for Windows)

#simulate some Gaussian responses

set.seed(040590)
X1= as.data.frame(array(data = rbinom(n = 50*1000,size = 1,prob = 0.3),dim = c(1000,50)))
Y1=rnorm(1000,-0.7+1*(X1$V1)*(X1$V4) + 1*(X1$V8*X1$V11)+1*(X1$V5*X1$V9),1)
X1$Y1 = Y1


formula1 = as.formula(paste(colnames(X1)[51],"~ 1 +",paste0(colnames(X1)[-c(51)],collapse = "+")))
data.example = as.data.frame(X1)

#run the inference with robust g prior
res = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)


#run the inference with Jeffrey's prior
res = LogicRegr(formula = formula1,data = data.example,family = "Gaussian",prior = "J",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)



#change to Bernoulli responses
X1$Y1=round(1.0/(1.0+exp(-Y1)))


#run the inference with robust g prior
res = LogicRegr(formula = formula1,data = data.example,family = "Bernoulli",prior = "G",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)


#run the inference with Jeffrey's prior
res = LogicRegr(formula = formula1,data = data.example,family = "Bernoulli",prior = "J",report.level = 0.5,d = 15,cmax = 2,kmax = 15,p.and = 0.9,p.not = 0.01,p.surv = 0.2,ncores = 32)
