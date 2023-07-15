#load some required libraries
library(mfp)
library(bfp)



#define your working directory, where the data files are stored
workdir=""


data.example = read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/abalone%20age/abalone.data",header = F)
data.example$MS=as.integer(data.example$V1=="M")
data.example$FS=as.integer(data.example$V1=="F")
data.example$V1=data.example$V9
data.example$V9 = NULL

names(data.example) = c("Age","Length", "Diameter","Height","WholeWeight","ShuckedWeight","VisceraWeight","ShellWeight","Male","Femele")

set.seed(040590)
teid =  read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/supplementaries/BGNLM/abalone%20age/teid.csv",sep = ";")[,2]


test = data.example[teid,]
data.example = data.example[-teid,]
train = data.example


vect = NULL
vect$x = as.matrix(data.example[,-1])
vect$y = data.example$Age

sum(test$Age)


gc()

#prepare the data structures for the final results
results=array(0,dim = c(2,100,5))

features = names(train)[-1]
for(ii in 1:100)
{
  

  
  
  print(paste("iteration ",ii))
  
  set.seed(ii)
  t=system.time({
    
    formula1 = paste(paste(colnames(data.example)[1], ' ~ 1 + fp(',paste0(colnames(data.example)[-c(1)], collapse = ') + fp(')))
    formula1 = as.formula(paste(formula1, ')', sep = ''))
    mod.mfp = mfp(formula1, data = data.example, family = gaussian, select = 0.05)
    
    #formula2 = paste(paste(colnames(data.example)[1], ' ~ 1 + bfp(',paste0(colnames(data.example)[-c(1)], collapse = ') + bfp(')))
    #formula2 = as.formula(paste(formula2, ')', sep = ''))
    #mod.mfp = BayesMfp(formula2, data = data.example, family = gaussian, priorSpecs = list(a = 4, modelPrior = "dependent"),
    #                   method = 'sampling')
    
  })
  # Predict
  results[1,ii,4]=t[3]
  
  
  
  t=system.time({
    out =predict(test,object = mod.mfp)
  })
  results[1,ii,5]=t[3]
  
  
  
  #compute and store the performance metrics
  results[1,ii,1]= sqrt(mean((out - test$Age)^2))
  results[1,ii,2]=mean(abs(out - test$Age))
  results[1,ii,3] = cor(out,test$Age)
  
  set.seed(ii)
  
  t=system.time({
    formula2 = paste(paste(colnames(data.example)[1], ' ~ 1 + bfp(',paste0(colnames(data.example)[-c(1)], collapse = ') + bfp(')))
    formula2 = as.formula(paste(formula2, ')', sep = ''))
    mod.mfp = BayesMfp(formula2, data = data.example, family = gaussian, priorSpecs = list(a = 4, modelPrior = "sparse"),
                       method = 'sampling')
    
  })
  # Predict
  results[2,ii,4]=t[3]
  
  
  
  t=system.time({
    out =predict(test,object = mod.mfp)
  })
  results[2,ii,5]=t[3]
  
  
  
  #compute and store the performance metrics
  results[2,ii,1]= sqrt(mean((out - test$Age)^2))
  results[2,ii,2]=mean(abs(out - test$Age))
  results[2,ii,3] = cor(out,test$Age)
  
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

summary.results = summary.results[1:2,]

round(summary.results[,c(2,1,3,5,4,6,8,7,9)],4)

names(summary.results)=c("min(rmse)","median(rmse)","max(rmse)","min(mae)","median(mae)","max(mae)","min(corr)","median(corr)","max(cor)","min(ltime)","median(ltime)","max(ltime)","min(ptime)","median(ptime)","max(ptime)")
rownames(summary.results)=c("lXGBOOST(logLik)","tXGBOOST(logLik)","LASSO","RIDGE","RFOREST","DEEPNETS","GR","VARBAYESS")

write.csv(x = summary.results,file = "summarycompete.csv")
