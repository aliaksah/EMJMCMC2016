#load some required libraries
library(RCurl)
library(glmnet)
library(xgboost)
library(h2o)
library(BAS)
#define your working directory, where the data files are stored
workdir=""


#prepare the test set data
simx = read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NEAs.txt"),sep = ",",header = T,fill=TRUE)
simy =  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Recognize/NotNeas8%2B.txt"),sep = ",",header = T,fill=TRUE)
simx$neo=1
simy$neo=0
data.example1 = as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = T)
transform=colnames(data.example1)[-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25)]
nas=NULL
for(i in 1:length(transform))
{
  print(i)
  data.example1[[transform[i]]]=as.numeric(as.character(data.example1[[transform[i]]]))
  nas=c(nas,which(is.na(data.example1[[transform[i]]])))
}
data.example1=data.example1[-unique(nas),]



#prepare the training set data
simx = read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NeoPHA.txt"),sep = ",",header = T,fill=TRUE)
simy =  read.table(text = getURL("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/examples/asteroid%20data/Teach/NotNeo-Type7.txt"),sep = ",",header = T,fill=TRUE)
simx$neo=1
simy$neo=0
data.example = as.data.frame(t(cbind(t(simy),t(simx))),stringsAsFactors = T)
for(i in 1:length(transform))
{
  print(i)
  data.example[[transform[i]]]=as.numeric(as.character(data.example[[transform[i]]]))
}


gc()

#prepare the data structures for the final results
results=array(0,dim = c(11,100,5))

# h2o initialize
h2o.init(nthreads=-1, max_mem_size = "6G")
h2o.removeAll()

for(ii in 1:100)
{
  print(paste("iteration ",ii))
  capture.output({withRestarts(tryCatch(capture.output({

    set.seed(ii)
    #here we are no longer running DBRM, since DBRM algorithms are run via other scripts
    #for computational efficiency and speed



    # xgboost logLik gblinear
    t=system.time({
      param = list(objective = "binary:logistic",
                   eval_metric = "logloss",
                   booster = "gblinear",
                   eta = 0.05,
                   subsample = 0.86,
                   colsample_bytree = 0.92,
                   colsample_bylevel = 0.9,
                   min_child_weight = 0,
                   gamma = 0.005,
                   max_depth = 15)

      train=as.data.frame(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])
      test=as.data.frame(data.example1[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])

      dval=xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
      watchlist=list(dval=dval)

      m2 = xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA),
                     param, nrounds = 10000,
                     watchlist = watchlist,
                     print_every_n = 10)

    })
    results[3,ii,4]=t[3]
    # Predict
    t=system.time({
      dtest  = xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
    })

    t=system.time({
      out = predict(m2, dtest)
    })

    results[3,ii,5]=t[3]
    out=as.integer(out>=0.5)

    #compute and store the performance metrics
    print(results[3,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

    #FNR
    ps=which(data.example1$neo==1)
    results[3,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[3,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


    # xgboost logLik gbtree
    t=system.time({
      param = list(objective = "binary:logistic",
                   eval_metric = "logloss",
                   booster = "gbtree",
                   eta = 0.05,
                   subsample = 0.86,
                   colsample_bytree = 0.92,
                   colsample_bylevel = 0.9,
                   min_child_weight = 0,
                   gamma = 0.005,
                   max_depth = 15)

      train=as.data.frame(data.example[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])
      test=as.data.frame(data.example1[,-c(2,4,5,13,14,15,16,17,19,20,21,22,23,24,25,37,38)])

      dval=xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA)
      watchlist=list(dval=dval)


      m2 = xgb.train(data = xgb.DMatrix(data = data.matrix(train[,-1]), label = data.matrix(train[,1]),missing=NA),
                     param, nrounds = 10000,
                     watchlist = watchlist,
                     print_every_n = 10)

    })

    results[4,ii,4]=t[3]
    # Predict
    system.time({
      dtest  = xgb.DMatrix(data.matrix(test[,-1]),missing=NA)
    })

    t=system.time({
      out = predict(m2, dtest)
    })
    out=as.integer(out>=0.5)

    #compute and store the performance metrics
    print(results[4,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))
    results[4,ii,5]=t[3]
    #FNR
    ps=which(data.example1$neo==1)
    results[4,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[4,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))



    #GLMNET (elastic networks) # lasso a=1
    t=system.time({
      fit2 = glmnet(as.matrix(train)[,-1], train$neo, family="binomial")
    })
    results[5,ii,4]=t[3]

    mmm=as.matrix(test[,-1])
    mmm[which(is.na(mmm))]=0
    #Predict
    t=system.time({
      out = predict(fit2,mmm , type = "response")[,fit2$dim[2]]
    })
    results[5,ii,5]=t[3]

    out=as.integer(out>=0.5)

    #compute and store the performance metrics
    print(results[5,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

    #FNR
    ps=which(data.example1$neo==1)
    results[5,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[5,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

    # ridge a=0
    t=system.time({
      fit2 = glmnet(as.matrix(train)[,-1], train$neo, family="binomial",alpha=0)
    })
    results[6,ii,4]=t[3]

    mmm=as.matrix(test[,-1])
    mmm[which(is.na(mmm))]=0
    #Predict
    t=system.time({
      out = predict(fit2,mmm , type = "response")[,fit2$dim[2]]
    })

    results[6,ii,5]=t[3]

    out=as.integer(out>=0.5)

    #compute and store the performance metrics
    print(results[6,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

    #FNR
    ps=which(data.example1$neo==1)
    results[6,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[6,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

    gc()

    # h2o.random forest



    df = as.h2o(train)



    train1 = h2o.assign(df , "train1.hex")
    valid1 = h2o.assign(df , "valid1.hex")
    test1 = h2o.assign(as.h2o(test[,-1]), "test1.hex")

    train1[1:5,]

    features = names(train1)[-1]

    # in order to make the classification prediction
    train1$neo = as.factor(train1$neo)

    t=system.time({
      rf1 = h2o.randomForest( stopping_metric = "AUC",
                              training_frame = train1,
                              validation_frame = valid1,
                              x=features,
                              y="neo",
                              model_id = "rf1",
                              ntrees = 10000,
                              stopping_rounds = 3,
                              score_each_iteration = T,
                              ignore_const_cols = T,
                              seed = ii)
    })
    results[7,ii,4]=t[3]

    #Predict
    t=system.time({
      out=h2o.predict(rf1,as.h2o(test1))[,1]
    })
    results[7,ii,5]=t[3]
    out=as.data.frame(out)

    #compute and store the performance metrics
    out=as.integer(as.numeric(as.character(out$predict)))


    print(results[7,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

    #FNR
    ps=which(data.example1$neo==1)
    results[7,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[7,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

    #h2o deeplearning

    t=system.time({
      neo.dl = h2o.deeplearning(x = features, y = "neo",hidden=c(200,200,200,200,200,200),
                                distribution = "bernoulli",
                                training_frame = train1,
                                validation_frame = valid1,
                                seed = ii)
    })
    # now make a prediction

    results[8,ii,4]=t[3]

    #Predict
    t=system.time({
      out=h2o.predict(neo.dl,as.h2o(test1))[,1]
    })
    results[8,ii,5]=t[3]
    out=as.data.frame(out)

    out=as.integer(as.numeric(as.character(out$predict)))

    #compute and store the performance metrics
    print(results[8,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

    #FNR
    ps=which(data.example1$neo==1)
    results[8,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[8,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


    #h2o glm

    t=system.time({
      neo.glm = h2o.glm(x = features, y = "neo",
                        family = "binomial",
                        training_frame = train1,
                        validation_frame = valid1,
                        #lambda = 0,
                        #alpha = 0,
                        lambda_search = F,
                        seed = ii)
    })
    #Predict
    results[9,ii,4]=t[3]

    t=system.time({
      out=h2o.predict(neo.glm,as.h2o(test1))[,1]
    })
    results[9,ii,5]=t[3]
    out=as.data.frame(out)

    out=as.integer(as.numeric(as.character(out$predict)))

    #compute and store the performance metrics
    print(results[9,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

    #FNR
    ps=which(data.example1$neo==1)
    results[9,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[9,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))

    #h2o naive bayes
    t=system.time({
      neo.nb = h2o.naiveBayes(x = features, y = "neo",
                              training_frame = train1,
                              validation_frame = valid1,
                              seed = ii)
    })
    #Predict
    results[10,ii,4]=t[3]
    t=system.time({
      out=h2o.predict(neo.nb,as.h2o(test1))[,1]
    })
    results[10,ii,5]=t[3]
    out=as.data.frame(out)

    out=as.integer(as.numeric(as.character(out$predict)))

    #compute and store the performance metrics
    print(results[10,ii,1]=(1-sum(abs(out-test$neo[1:length(out)]))/length(out)))

    #FNR
    ps=which(data.example1$neo==1)
    results[10,ii,2]=sum(abs(out[ps]-data.example1$neo[ps]))/(sum(abs(out[ps]-data.example1$neo[ps]))+length(ps))

    #FPR
    ns=which(data.example1$neo==0)
    results[10,ii,3]=sum(abs(out[ns]-data.example1$neo[ns]))/(sum(abs(out[ns]-data.example1$neo[ns]))+length(ns))


  })), abort = function(){onerr=TRUE;out=NULL})})
  print(results[,ii,1])
}

ids=1:100

ress=results[,ids,]


#make the joint summary of the runs, including min, max and medians of the performance metrics
summary.results=array(data = NA,dim = c(15,15))
for(i in 1:11)
{
  for(j in 1:5)
  {
    summary.results[i,(j-1)*3+1]=min(ress[i,,j])
    summary.results[i,(j-1)*3+2]=median(ress[i,,j])
    summary.results[i,(j-1)*3+3]=max(ress[i,,j])
  }
}
summary.results=as.data.frame(summary.results)
names(summary.results)=c("min(prec)","median(prec)","max(prec)","min(fnr)","median(fnr)","max(fnr)","min(fpr)","median(fpr)","max(fpr)","min(ltime)","median(ltime)","max(ltime)","min(ptime)","median(ptime)","max(ptime)")
rownames(summary.results)[1:11]=c("GMJMCMC(AIC)","MJMCMC(AIC)","lXGBOOST(logLik)","tXGBOOST(logLik)","LASSO","RIDGE","RFOREST","DEEPNETS","LR","NAIVEBAYESS","KMEANS")



#write the final reults into the files
write.csv(x = round(summary.results,4),file = "/mn/sarpanitu/ansatte-u2/aliaksah/Desktop/package/EMJMCMC/examples/neo classification/asteroidsfinal.csv")


