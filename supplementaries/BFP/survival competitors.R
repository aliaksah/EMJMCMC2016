rm(list=ls())
setwd("~/fractional_polynomials")

# get data
download.file('https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/Dateien/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/gbsg_br_ca.zip',
              'gbsg_br_ca.zip')
data <- read.csv(unz('gbsg_br_ca.zip',
                     'gbsg_br_ca/gbsg_br_ca.csv'),
                 header = TRUE)
#system('rm whitehall1.zip')

# load R libraries:
# - mfp: original fractional polynomial
# install.packages('mfp')
# note, bfp does not have the option to perform survival analysis

library(mfp)
#library(pec)
library(rms)

head(data)
time <- data$rectime
cens <- data$censrec
X <- data[, c(2:4, 6:8, 10:12)]
data <- cbind(time, cens, X)

trainAndTest <- function(id, data)
{
  set.seed(id)
  # split in training and test set
  train <- c(sample((1:nrow(data))[data$cens == 1], sum(data$cens)*2/3), # split separately events
             sample((1:nrow(data))[data$cens == 0], sum(!data$cens)*2/3)) # and censored observations
  # fit the model
  model <- mfp(formula = Surv(time, cens) ~ fp(age, df = 4, select = 0.05) + fp(size, df = 4, select = 0.05) +
                 fp(nodes, df = 4, select = 0.05) + fp(pgr, df = 4, select = 0.05) + fp(er, df = 4, select = 0.05) +
                 fp(meno, df = 4, select = 0.05) + fp(gradd1, df = 4, select = 0.05) + fp(gradd2, df = 4, select = 0.05) +
                 fp(hormon, df = 4, select = 0.05), data = data[train, ], family = 'cox')
  # change the variables otherwise pec does not work
  selected <- rownames(model$powers)[apply(!is.na(model$powers), 1, sum) > 0]
  # TBD: make this faster
  X.mod <- NULL
  for(i in colnames(X))
  {
    if(!is.na(model$powers[i, 1]))
    {
      x.rescaled <- ((X[, i] + model$scale[i, 1])/model$scale[i, 2])
      tmp <- x.rescaled^model$powers[i, 1]
      if(model$powers[i, 1] == 0) tmp <- log(x.rescaled)
      if(!is.na(model$powers[i, 2]))
      {
        if(model$powers[i, 2] == model$powers[i, 1])
        {
          tmp <- cbind(tmp, x.rescaled^model$powers[i, 2] * log(x.rescaled))
          if(model$powers[i, 2] == 0) tmp[, 2] <- log(x.rescaled)*log(x.rescaled)
        }
        if(model$powers[i, 2] != model$powers[i, 1])
        {
          tmp <- cbind(tmp, x.rescaled^model$powers[i, 2])
          if(model$powers[i, 2] == 0) tmp[, 2] <- log(x.rescaled)
        }
      }
      X.mod <- cbind(X.mod, tmp)
      colnames(X.mod)[ncol(X.mod)] <- paste0(i, '.1')
      if(!is.na(model$powers[i, 2])) colnames(X.mod)[(ncol(X.mod) - 1):ncol(X.mod)] <- paste0(i, c('.1', '.2'))
    }
  }
  data.mod <- as.data.frame(cbind(time, cens, X.mod))
  # to use pec, fit a Cox model with the edited columns
  model.cox <- coxph(Surv(time, cens) ~ ., data = data.mod[train, ], x = TRUE)
  pred_fp <- pec(model.cox, model.cox$formula, data = data.mod[-train, ], traindata = data.mod[train, ], cens.model = 'cox')
  ibs_fp <- ibs(pred_fp, times = min(max(data$time[train]), max(data$time[-train])))
  cindex_fp <- cindex(model.cox, model.cox$formula, data = data.mod[-train, ], cens.model = 'cox')$AppCindex
  # model linear effects
  mod.cox <- cph(Surv(time, cens) ~ age + size + nodes + pgr + er + meno + gradd1 + gradd2 + hormon, data = data[train,])
  selected <- fastbw(mod.cox, rule = 'p')$names.kept
  model.lin <- coxph(Surv(time, cens) ~ ., data = data[train, c('time', 'cens', selected)], x = TRUE)
  cindex_lin <- cindex(model.lin, model.lin$formula, data = data[-train, ], cens.model = 'cox')$AppCindex
  pred_lin <- pec(model.lin, model.lin$formula, data = data[-train, ], traindata = data[train, ], cens.model = 'cox')
  ibs_lin <- ibs(pred_lin, times = min(max(data$time[train]), max(data$time[-train])))
  
  print(id)
  rbind(model$power, as.vector(ibs_fp), as.vector(ibs_lin), c(cindex_fp, cindex_lin))
}

res <- sapply(c(1:906, 908:975, 977:1002), trainAndTest, data = data, simplify = 'array')
save(res, file = "results_gbsg.Rdata")

p <- ncol(X)
FP_form <- res[1:p, , ]
tmp <- apply(FP_form, c(1, 3), function(x) {x <- unlist(x); x[is.na(x)] <- 9; x[1] + 0.01 * x[2]})
recognise_FP <- function(x)
{
  if(is.na(x$power1)) return('not selected')
  if(is.na(x$power2))
  {
    if(x$power1 == 1) return('linear')
    if(x$power1 != 1) return(paste0('FP_1(',x$power1,')'))
  }
  return(paste0('FP_2(',x$power1,';',x$power2,')'))
}
tmp <- apply(FP_form, c(1, 3), recognise_FP)
FP <- apply(tmp, 1, table)

cindex_lin <- mean(unlist(res[p + 3, 1, ])) # res[p + 2, 2, ] should be always around 0.5
cindex_fp <- mean(unlist(res[p + 3, 2, ])) # res[p + 2, 2, ] should be always around 0.5
ibs_null <- mean(unlist(res[p + 1, 1, ], res[p + 1, 2, ])) #temporary, try to understand why I get two different values
ibs_mfp <- mean(unlist(res[p + 1, 2, ]))
ibs_lin <- mean(unlist(res[p + 2, 2, ]))

c(cindex_lin, cindex_fp, ibs_lin, ibs_mfp, ibs_null)
