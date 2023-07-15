#load some required libraries

library(mfp)
library(bfp)


#read in the train and test data sets
test = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/breast%20cancer/test.csv",header = T,sep=",")[,-1]
train = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/breast%20cancer/train.csv",header = T,sep=",")[,-1]

#transform the train data set to a data.example data.frame that EMJMCMC class will internally use
data.example = as.data.frame(train)
#perfrom garbage collection
gc()

#prepare the data structures for the final results
results=array(0,dim = c(2,100,5))

for(ii in 1:100)
{
  
  
  print(paste("iteration ",ii))
  
  t=system.time({
    formula1 = paste(paste(colnames(data.example)[31], ' ~ 1 + fp(',paste0(colnames(data.example)[-c(31)], collapse = ') + fp(')))
    formula1 = as.formula(paste(formula1, ')', sep = ''))
    mod.mfp = mfp(formula1, data = data.example, family = binomial, select = 0.05)
    
  })
  # Predict
  results[1,ii,4]=t[3]
  
  
  
  t=system.time({
    out =predict(test,object = mod.mfp)
  })
  results[1,ii,5]=t[3]
  
  out = as.integer(out>=0.5)
  
  ps=which(test$X==1)
  
  
  ns=which(test$X==0)
  
  #compute and store the performance metrics
  results[1,ii,1]= (1-sum(abs(out-test$X[1:length(out)]))/length(out))
  results[1,ii,2]= sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))
  results[1,ii,3] = sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))
  
  
  t=system.time({
    formula2 = paste(paste(colnames(data.example)[31], ' ~ 1 + bfp(',paste0(colnames(data.example)[-c(31)], collapse = ') + bfp(')))
    formula2 = as.formula(paste(formula2, ')', sep = ''))
    mod.mfp = BayesMfp(formula2, data = data.example, family = gaussian, priorSpecs = list(a = 4, modelPrior = "flat"),
                       method = 'sampling')
    
  })
  # Predict
  results[2,ii,4]=t[3]

  
  
  
  t=system.time({
    out =predict(test,object = mod.mfp)
  })
  results[2,ii,5]=t[3]
  
  
  out = as.integer(out>=0.5)
  
  #compute and store the performance metrics
  results[2,ii,1]= (1-sum(abs(out-test$X[1:length(out)]))/length(out))
  results[2,ii,2]= sum(abs(out[ps]-test$X[ps]))/(sum(abs(out[ps]-test$X[ps]))+length(ps))
  results[2,ii,3] = sum(abs(out[ns]-test$X[ns]))/(sum(abs(out[ns]-test$X[ns]))+length(ns))
  print( results[,ii,1])
  
  gc()
  #})), abort = function(){onerr=TRUE;out=NULL})})
}

ids=1:100

ress=results[,ids,]

#make the joint summary of the runs, including min, max and medians of the performance metrics
summary.results=array(data = NA,dim = c(15,15))
for(i in 1:2)
{
  for(j in 1:5)
  {
    summary.results[i,(j-1)*3+1]=min(ress[i,,j])
    summary.results[i,(j-1)*3+2]=median(ress[i,,j])
    summary.results[i,(j-1)*3+3]=max(ress[i,,j])
  }
}
summary.results=as.data.frame(summary.results)

View(summary.results[,c(2,1,3,5,4,6,8,7,9)])
write.csv(x = summary.results,file = "summarycompetebrcancer.csv")

