#check the tutorial

library(EMJMCMC)
library(clusterGeneration)
library(bindata)
library(ggplot2)

# set the seed
set.seed(040590)
# construct a correlation matrix for M = 50 variables
M = 50
m = clusterGeneration::rcorrmatrix(M,alphad=2.5) 
# simulate 1000 binary variables with this correlation matrix
X = bindata::rmvbin(1000, margprob = rep(0.5,M), bincorr = m)

# prepare the correlation matrix in the melted format
melted_cormat = reshape2::melt(cor(X))
# plot the heat-map of the correlations
ggplot2::ggplot(data = melted_cormat, 
                ggplot2::aes(x=Var1, y=Var2, fill=value)) + 
  ggplot2::geom_tile() +  
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.title.y =  ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank())

# simulate Gaussian responses from a model with two-way interactions
# which is considered in S.4 of the paper
df = data.frame(X)
df$Y=rnorm(n = 1000,mean = 1+1.43*(df$X5*df$X9)+
             0.89*(df$X8*df$X11)+0.7*(df$X1*df$X4),sd = 1)

help("LogicRegr")


# specify the initial formula
formula1 = as.formula(paste(colnames(df)[M+1],"~ 1 + ",
                            paste0(colnames(df)[-c(M+1)],collapse = "+")))
# Bayesian logic regression with the robust-g-prior
res4G = LogicRegr(formula = formula1, data = df,
                  family = "Gaussian", prior = "G", report.level = 0.5, 
                  d = 15,cmax = 2,kmax = 15, p.and = 0.9, p.not = 0.1, p.surv = 0.2, 
                  ncores = 32)
# Bayesian logic regression with the Jeffreys prior
res4J = LogicRegr(formula = formula1, data = df,
                  family = "Gaussian", prior = "J", report.level = 0.5, 
                  d = 15, cmax = 2,kmax = 15, p.and = 0.9, p.not = 0.1, p.surv = 0.2, 
                  ncores = 32)

# print the results for the robust g-prior
print(base::rbind(c("expressions","probabilities"),res4G$feat.stat))

#print the results for the Jeffreys prior
print(base::rbind(c("expressions","probabilities"),res4J$feat.stat))


# simulate Gaussian responses from a model with two-way interactions 
# and an age effect which is an extension of S.4 of the paper
Xp = data.frame(X)
Xp$age = rpois(1000,lambda = 34)
Xp$Y=rnorm(n = 1000,mean = 1+0.7*(Xp$X1*Xp$X4) + 
             0.89*(Xp$X8*Xp$X11)+1.43*(Xp$X5*Xp$X9) + 2*Xp$age, sd = 1)


teid  = sample.int(size =100,n = 1000,replace = F)
test  = Xp[teid,]
train = Xp[-teid,]



# specify the initial formula
formula1 = as.formula(paste("Y ~ 1 +",
                            paste0(colnames(test)[-c(51,52)],collapse = "+")))
# specify the link function
g = function(x) x
# specify the parameters of the custom estimator function
estimator.args = list(data = train, n = dim(train)[1], 
                      m =stri_count_fixed(as.character(formula1)[3],"+"),k.max = 15)
# specify the parameters of gmjmcmc algorithm
gmjmcmc.params = list(allow_offsprings=1,mutation_rate = 250,
                      last.mutation=10000, max.tree.size = 5, Nvars.max =15,
                      p.allow.replace=0.9,p.allow.tree=0.01,p.nor=0.01,p.and = 0.9)
# specify some advenced parameters of mjmcmc
mjmcmc.params = list(max.N.glob=10, min.N.glob=5, max.N=3, min.N=1,
                     printable = F)
# run the inference of BLR with a non-binary covariate and predicions
res.alt = pinferunemjmcmc(n.cores = 30, report.level =  0.2, 
                          num.mod.best = 100,simplify = T,predict = T,test.data = test,
                          link.function = g,
                          runemjmcmc.params = list(formula = formula1,latnames = c("I(age)"),
                                                   data = train,estimator = estimate.logic.lm.jef,
                                                   estimator.args =estimator.args,
                                                   recalc_margin = 249, save.beta = T,interact = T,outgraphs=F,
                                                   interact.param = gmjmcmc.params,
                                                   n.models = 10000,unique = T,max.cpu = 4,max.cpu.glob = 4,
                                                   create.table = F,create.hash = T,pseudo.paral = T,burn.in = 100,
                                                   print.freq = 1000,
                                                   advanced.param = mjmcmc.params))


print(base::rbind(c("expressions","probabilities"),res.alt$feat.stat))



print(sqrt(mean((res.alt$predictions-test$Y)^2)))
print(mean(abs((res.alt$predictions-test$Y))))


library(HDeconometrics)
ridge = ic.glmnet(x = train[,-51],y=train$Y,family = "gaussian",
                  alpha = 0)
predict.ridge = predict(ridge$glmnet,newx = as.matrix(test[,-51]),
                        type = "response")[,which(ridge$glmnet$lambda == ridge$lambda)]
print(sqrt(mean((predict.ridge-test$Y)^2)))
print(mean(abs((predict.ridge-test$Y))))


tmean = 1+2*test$age+0.7*(test$X1*test$X4) + 
  0.89*(test$X8*test$X11)+1.43*(test$X5*test$X9)
print(sqrt(mean((tmean -test$Y)^2)))
print(mean(abs((tmean -test$Y))))

