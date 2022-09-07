
# Simulation: 
  # It is divided into three parts. 
  # The first part is data generation, 
  # the second part is logistic model comparison, 
  # and the third part is variable selection and classification.
# The same goes for real data.

# Some packages
library(gbm)
library('e1071')   #SVM
library("kernlab")
library(tree)
library(randomForest)
library(xgboost)
# library(Matrix)   
library(ROCR)     #  AUC 
library(EnvStats) #  rnormTrunc()
library(mvtnorm)  
library(glmnet)
library(mniw)     #  Inverse_Wishart
library(ggplot2)
library(reshape2)   


######### Simulated data generation and label generation；#############
m_len = 1
n_person = 100        # Individual number
n_variable = 200      # Number of variables
label_moni = matrix(0,m_len,n_person)  # label 
# set.seed(5)
x = sample2_xn(n_person,n_variable,0,0.9)
n = dim(x)[1]
m = dim(x)[2]
e = 0.1
n_cor = 5

aa <- 1
while(aa == 1){
  for(i in 1:m_len){
    ## Artificially constructed labels, two categories
    # x_01 = matrix(0,nrow = n_person,ncol = 10)
    x_gyh <- apply(x[,1:(e*m)],2,gyh1)
    cor_x_01 = cor(x[,1:(e*m)],method = "spearman")  # Correlation matrix
    label_train = sample(c(1,-1),size=n,replace = TRUE,prob = c(500,500))
    pre1 = label_train
    # par(mfrow=c(2,3))
    for(k in 1:4){
      label_train = ifelse(pre1<=0,-1,1)
      x_01 = Calculate_Bij(x[,1:(e*m)],label_train,n_cor)
      x_gen_label = as.data.frame(x_gyh*x_01)  
      p_train = ifelse(label_train==-1,0,1)
      # glm.fit = glm(p_train~.,family = binomial,x_gen_label,control=list(maxit=200))   # LR
      # result_z <- c(AIC(glm.fit),BIC(glm.fit),0,1)
      # glm.pre1 = predict(glm.fit,x_gen_label,type="response") # Threshold selection
      # # plot(glm.pre1,xlab="individuals",ylab="Pr(y=1|Emat)",cex.lab=1.5)   #plot(glm.pre1[pre1==1]) 
      # #plot(sort(glm.pre1))  
      # pre1 <- ifelse(glm.pre1 > median(glm.pre1),1,0)
      # label_moni[i,] = pre1
      
      label_train = ifelse(pre1<=0,-1,1)
      x_01 = Calculate_Bij(x[,1:(e*m)],label_train,n_cor)
      X1 <- as.matrix(as.data.frame(x_gyh*x_01)) 
      p_train = ifelse(label_train==-1,0,1)
      pfit1 = cv.glmnet(X1,p_train, family = "binomial",type.measure ="class",nfolds = 5)
      pfit11 <- glmnet(X1,p_train, family = "binomial",lambda = pfit1$lambda.min)
      # pfit11 <- glmnet(X1,p_train, family = "binomial")
      pre11 = predict(pfit11,X1,s = pfit11$lambda,type = "response")
      # result_z <- c(BICAICglm(pfit11),1)
      
      # plot(pre11,xlab="individuals",ylab="Pr(y=1|Emat)",cex.lab=1.5)
      pre11 <- ifelse(pre11 > median(pre11),1,0)
      aa <- 0
      if(min(table(pre11))<n_person/4 | length(table(pre11))==1){
        # pre11 <- pre1
        pre11 <- sample(c(1,0),size=n,replace = TRUE,prob = c(500,500))
        aa <- 1
      }
      pre1 <- pre11
      label_moni[1,] = pre1
    }
  }
}
pre11 = apply(label_moni,2,mean,na.rm=TRUE)

table(pre11)
x = subset(x,pre11!=0.5) 
pre11 = pre11[pre11!=0.5]
pre1 = ifelse(pre11>0.5,1,0)
sum(pre1)
sum(pre1)/length(pre1)

s_matrix = apply(x[,1:(e*m)],2,as.numeric)   

# Divide the screening sample into two parts, CRC (1) or HC (0).
s_matrix = x
s_matrix_HC = subset(s_matrix,pre1 ==0)   
s_matrix_CRC = subset(s_matrix,pre1 ==1) 
cor_x = cor(x)

n_sample_HC = sample(nrow(s_matrix_HC))
n_sample_CRC = sample(nrow(s_matrix_CRC))
s_matrix_HC <- s_matrix_HC[n_sample_HC,]
s_matrix_CRC <- s_matrix_CRC[n_sample_CRC,]

# x_01_HC = subset(x_01,pre1 ==0)  
# x_01_CRC = subset(x_01,pre1 ==1) 
# x_01_HC <- x_01_HC[n_sample_HC,]
# x_01_CRC  <- x_01_CRC[n_sample_CRC,]

######## Model to compare；###################

#  Simulation 1, the real model is 2-order, 3-order;;
a1 <- list()
vari <- seq(50,250,5)  
qs <- matrix(0,2,length(vari))
t3 <- Sys.time() 
for(q in 1:10){ 
  t = 1
  for(g in vari){
    i = 1
    m_len = 1
    n_person = 100       # Individual number
    n_variable = g       # Number of variables 2 5 8 10 15 20
    label_moni = matrix(0,m_len,n_person)  # 1000 individual labels
    # set.seed(5)
    # x = sample2_x(n_person,n_variable)
    # x = sample2_x1(n_person,n_variable,0,0)
    x = sample2_xn(n_person,n_variable,0,0.9)
    # x = as.data.frame(x)
    # x[1,]
    n = dim(x)[1]
    m = dim(x)[2]
    e = 0.1   # 0.1 or 0.5 or 1；
    n_cor = 5
    
    x01 <- x[,1:(g*e)]   # as.data.frame(apply(x,2,gyh1))
    n_items <- (g*e)/2  #  choose(n_variable,2)/2  #n_variable
    X2_i = Generate_interaction_items(x01,2,n_items) 
    X3_i = Generate_interaction_items(x01,3,n_items) 
    # X4_i = Generate_interaction_items(x01,4,n_items)  
    X1_2 <- cbind(x01,X2_i)
    X1_3 <- cbind(x01,X3_i)
    BICACC <- c()
    # data1 <-  list(X1_2,X1_3,X1_4,X1_23,X1_24,X1_34,X1_234)
    data1 <-  list(X1_2,X1_3)
    
    for(ii in 1:length(data1)){
      X1 <-  as.data.frame(data1[[ii]]) # X1_2, X1_3, X1_23, X1_234
      X1 = as.matrix(X1)
      result_z <- c()   # Record BIC and ACC of real model (1)
      t1 <- Sys.time() 
      # pre1 = sample(c(1,0),size=n,replace = TRUE,prob = c(500,500))
      # table(pre1)
      # par(mfrow=c(2,2))
      aa <- 1 # Avoid hashtag problems
      while(aa==1){
        pre1 = sample(c(1,0),size=n,replace = TRUE,prob = c(500,500))
        for(k in 1:4){
          p_train = pre1
          pfit1 = cv.glmnet(X1,p_train, family = "binomial",type.measure ="class",nfolds = 5)
          pfit11 <- glmnet(X1,p_train, family = "binomial",lambda = pfit1$lambda.min)
          # pfit11 <- glmnet(X1,p_train, family = "binomial")
          pre11 = predict(pfit11,X1,s = pfit11$lambda,type = "response") # 
          result_z <- c(BICAICglm(pfit11),1)
          # plot(pre11,xlab="individuals",ylab="Pr(y=1|Emat)",cex.lab=1.5)
          pre11 <- ifelse(pre11 > median(pre11),1,0)
          aa <- 0
          if(min(table(pre11))<n_person/4 | length(table(pre11))==1){
            # pre11 <- pre1
            pre11 <- sample(c(1,0),size=n,replace = TRUE,prob = c(500,500))
            aa <- 1
          }
          pre1 <- pre11
          label_moni[1,] = pre1
        }
      }
      pre11 = apply(label_moni,2,mean,na.rm=TRUE)
      
      data = x  # x original data or x01_ Original data with some disturbance
      label = pre11
      # e = 1
      # n_cor = 3
      output = model_comparison3(data,label,e,n_cor)
      t2 = Sys.time() 
      print(t2-t1)
      rbind_1 <- rbind(output[[1]],result_z)
      rownames(rbind_1) <- c('0','2','3','10','z')
      BICACC <- rbind(BICACC,rbind_1)
      
      rm(list=c('output','data','label','pre1','X1','p_train','rbind_1'))
      gc()
    }
    # BICACC 
    
    if(q==1){
      a1[[t]] = BICACC
    }else{
      a1[[t]] = a1[[t]] + BICACC
    }
    
    t = t+1
    
    rm(list=c('x','X2_i','X3_i','X1_2','X1_3'))
    gc()
  }
}
t4 = Sys.time() 
print(t4-t3)

# Collation of results for Logistic models with first-order and second-order interaction terms.
a1[[1]]
for(i in 1:length(a1)){
  a1[[i]] = a1[[i]]/q   
}
a11<- a1
AIC <- BIC <- ACC <- list()
for(i in 1:length(a11)){
  for(j in 1:length(data1)){
    if(i==1){
      AIC[[j]] <- a11[[i]][((j-1)*5+1):(5*j),1]   # AIC
      BIC[[j]] <- a11[[i]][((j-1)*5+1):(5*j),2]   # BIC
      ACC[[j]] <- a11[[i]][((j-1)*5+1):(5*j),4]   # ACC
    }else{ 
      AIC[[j]] <- rbind(AIC[[j]],a11[[i]][((j-1)*5+1):(5*j),1])   # AIC
      BIC[[j]] <- rbind(BIC[[j]],a11[[i]][((j-1)*5+1):(5*j),2])   # BIC
      ACC[[j]] <- rbind(ACC[[j]],a11[[i]][((j-1)*5+1):(5*j),4])   # ACC
    }
  }
}
# ACC[[1]]
# AIC[[1]]
# BIC[[1]]

# Simulation 2, real model based on Binary matrix;
# a2 <- a1
a2 <- list()
vari <- seq(50,250,5) # seq(10,120,10) 
qs <- rep(0,length(vari))
g_error <- c()
t3=Sys.time()
for(q in 1:10){
  t = 1
  for(g in vari){     #c(5,8,10,15,20,30,50,100,150)
    m_len = 1
    n_person = 100    #Number of observations
    n_variable = g       #Number of variables 2 5 8 10 15 20
    # Match_rate = matrix(0,m_len,4)
    # Match_rate2 = matrix(0,m_len,6)
    label_moni = matrix(0,m_len,n_person)  # 1000 individual labels
    # set.seed(5)
    # x = sample2_x(n_person,n_variable)
    x = sample2_xn(n_person,n_variable,0,0.9)
    n = dim(x)[1]
    m = dim(x)[2]
    e = 0.1     # 0.1 or 0.5 or 1；
    n_cor = 5   #3
    result_z <- c()
    
    t1=Sys.time()
    aa <- 1
    while(aa == 1){
      for(i in 1:m_len){
        ## Artificially constructed tags, two categories (unsupervised classification)
        # (+1,-1) matrix generation;
        
        # x_01 = matrix(0,nrow = n_person,ncol = 10)
        x_gyh <- apply(x[,1:(e*m)],2,gyh1)
        cor_x_01 = cor(x[,1:(e*m)],method = "spearman")  # Correlation matrix
        label_train = sample(c(1,-1),size=n,replace = TRUE,prob = c(500,500))
        pre1 = label_train
        # par(mfrow=c(2,3))
        for(k in 1:4){
          label_train = ifelse(pre1<=0,-1,1)
          x_01 = Calculate_Bij(x[,1:(e*m)],label_train,n_cor)
          x_gen_label = as.data.frame(x_gyh*x_01)  
          p_train = ifelse(label_train==-1,0,1)
          label_train = ifelse(pre1<=0,-1,1)
          x_01 = Calculate_Bij(x[,1:(e*m)],label_train,n_cor)
          X1 <- as.matrix(as.data.frame(x_gyh*x_01))  
          p_train = ifelse(label_train==-1,0,1)
          pfit1 = cv.glmnet(X1,p_train, family = "binomial",type.measure ="class",nfolds = 5)
          pfit11 <- glmnet(X1,p_train, family = "binomial",lambda = pfit1$lambda.min)
          # pfit11 <- glmnet(X1,p_train, family = "binomial")
          pre11 = predict(pfit11,X1,s = pfit11$lambda,type = "response") 
          result_z <- c(BICAICglm(pfit11),1)
          
          # plot(pre11,xlab="individuals",ylab="Pr(y=1|Emat)",cex.lab=1.5)
          pre11 <- ifelse(pre11 > median(pre11),1,0)
          aa <- 0
          if(min(table(pre11))<n_person/4 | length(table(pre11))==1){
            # pre11 <- pre1
            pre11 <- sample(c(1,0),size=n,replace = TRUE,prob = c(500,500))
            aa <- 1
          }
          pre1 <- pre11
          label_moni[1,] = pre1
        } 
      }
    }
    pre11 = apply(label_moni,2,mean,na.rm=TRUE)
    
    data = x
    label = pre11
    # interactions = 3 
    # output <- model_comparison(data,label,x_01,e,n_cor)
    output1 <- model_comparison11(data,label,x_01,e,n_cor)
    t2=Sys.time()
    print(t2-t1)
    
    rbind_1 <- rbind(output1[[1]],result_z)
    rownames(rbind_1) <- c('0','2','3','10','z')
    if(q==1){
      # a[[t]] = output[[1]]
      a2[[t]] = rbind_1# output1[[1]]
    }else{
      # a[[t]] = a[[t]] + output[[1]]
      a2[[t]] = a2[[t]] + rbind_1  # output1[[1]]
    }
    t = t+1
    # rm(list=c('output','output1','x','cor_x_01','x_01','x_gen_label',
    #           'pre11','data','label'))
    rm(list=c('output1','x','cor_x_01','x_01','x_gen_label',
              'pre11','data','label'))
    gc()
  }
}
t4=Sys.time()
print(t4-t3)

# Collation of results for Logistic models with Binary matrix.
a22 <- list()    
for(i in 1:length(a2)){
  # a0[[i]] = a[[i]]/q
  a22[[i]] = a2[[i]]/q
}
a22
# Collect the data first
AIC2 <- BIC2 <- ACC2 <- c()
for(i in 1:length(a22)){
  AIC2 <- rbind(AIC2,a22[[i]][,1])   # AIC
  BIC2 <- rbind(BIC2,a22[[i]][,2])   # BIC
  ACC2 <- rbind(ACC2,a22[[i]][,4])   # ACC
}
AIC2 <- (AIC2-min(AIC2))/(max(AIC2)-min(AIC2))
BIC2 <- (BIC2-min(BIC2))/(max(BIC2)-min(BIC2))
AIC[[3]] <- AIC2
BIC[[3]] <- BIC2
ACC[[3]] <- ACC2

# Drawing of model comparison
pdf(file="result\\model_pic.pdf",width=15, height = 4.5)
for(j in 1:3){
  z_model <- c('X_2-order','X_3-order','X_B')
  # j = 1
  AIC1 <- AIC[[j]]
  BIC1 <- BIC[[j]]
  ACC1 <- ACC[[j]]
  AIC1 <- (AIC1-min(AIC1))/(max(AIC1)-min(AIC1))
  BIC1 <- (BIC1-min(BIC1))/(max(BIC1)-min(BIC1))
  
  rownames(AIC1) <- rownames(BIC1) <- rownames(ACC1) <- vari
  colnames(AIC1) <- colnames(BIC1) <- colnames(ACC1) <- 
    c("Linear","2-order","3-order","LRBmat","Real model")
  
  AIC2 <- melt(AIC1)
  BIC2 <- melt(BIC1)
  ACC2 <- melt(ACC1)
  
  AIC3 <- cbind(AIC2,"AIC")
  BIC3 <- cbind(BIC2,"BIC")
  ACC3 <- cbind(ACC2,"ACC") 
  colnames(AIC3) <- colnames(BIC3) <- colnames(ACC3) <- 
    c("variable","model","value",'criterion')
  
  result <- rbind(AIC3,BIC3,ACC3)
  # result$value[is.na(result$value)]
  # red orange yellow malachite green bluish
  p <- ggplot(data = result,aes(x=variable,y=value,group = model,color=model))+
    geom_point()+
    geom_line()+
    # ylab(z_model[j]) + 
    ylab(NULL)+
    xlab("Number of covariables")+
    ylim(-0.05, 1.05)+
    facet_wrap(~ criterion, scales = "free")+
    scale_color_manual(values = c("#9ec417","#13a983","#44c1f0","#cc340c","#f18800"))+
    theme( axis.text = element_text( size = 12 ),
           # axis.text.y = element_blank(),
           legend.text = element_text( size = 14 ),
           axis.title = element_text( size = 14 ),
           strip.text = element_text(size = 18) )
  print(p)
} 
dev.off()


######## Variable selection and classification；###################
i = 1
leni = 20
r_len = 5    # CV = 5
CV_ncol = 4*3
n_cor = 5   
CV.xgboost <-CV.rf <-CV.logistic <- 
  CV.logistic2 <- CV.logistic3 <- CV.svm  <- matrix(0,nrow=leni,ncol=CV_ncol)
CV.xgboost1<-CV.rf1<-CV.logistic1 <- 
  CV.logistic12 <- CV.logistic13 <- CV.svm1 <- matrix(0,nrow=leni,ncol=CV_ncol)
Correct.xgboost <-Correct.rf <-Correct.logistic <- 
  Correct.logistic2 <- Correct.logistic3 <- Correct.svm  <- matrix(0,nrow=r_len,ncol=CV_ncol)
Correct.xgboost1<-Correct.rf1<-Correct.logistic1 <- 
  Correct.logistic12 <- Correct.logistic13 <- Correct.svm1 <- matrix(0,nrow=r_len,ncol=CV_ncol)
gamma_w = 0.1
cost_w = 10
mtry_w = 7
max_depth_w = 8
eta_w = 0.5

gamma1 = 0.1
cost1 = 10
mtry1 = 8
max_depth1 = 6
eta1 = 0.1

# Add the optimal parameters for training and simulate the data
if(TRUE){  
  t1=Sys.time() 
  for(k in 1:leni){
    # k = 1; i = 1
    # t21=Sys.time() 
    n_sample_HC = sample(nrow(s_matrix_HC))
    n_sample_CRC = sample(nrow(s_matrix_CRC))
    s_matrix_HC1 <- s_matrix_HC[n_sample_HC,]
    s_matrix_CRC1 <- s_matrix_CRC[n_sample_CRC,]
    # output_cvvs = CV_variable_selection(s_matrix_HC,s_matrix_CRC,5)
    # s_matrix_HC1 <- output_cvvs[[1]]
    # s_matrix_CRC1 <- output_cvvs[[2]]
    # cvvs = output_cvvs[[3]]
    # t22=Sys.time() 
    # print(t22-t21)
    for(i in 1:r_len){
      output_divide = divide_train_test(s_matrix_HC1,s_matrix_CRC1,r_len,i)
      train_data = output_divide[[1]]    
      test_data = output_divide[[2]]    
      label_train = output_divide[[3]]  
      label_test = output_divide[[4]]
      
      t3=Sys.time()
      output_bi = binary_matrix_conversion32(train_data,test_data,label_train,label_test,n_cor)
      # t4=Sys.time()
      # print(t4-t3)
      
      train_data1 = as.data.frame(train_data)   
      train_data <- as.data.frame(apply(train_data1,2,gyh1))
      #train_data = as.matrix(train_data1)
      test_data1 = as.data.frame(test_data)
      test_data <- as.data.frame(apply(test_data1,2,gyh1))
      
      train_data_01 = output_bi[[1]]  
      test_data_01 = output_bi[[2]]    
      train_data_w = train_data_01*train_data
      test_data_w = test_data_01*test_data
      
      # Training set result
      output_svm = classification_svm(train_data_w,train_data_w,label_train,label_train,gamma_w,cost_w)
      output_logistic = classification_logistic(train_data_w,train_data_w,label_train,label_train)
      output_logistic2 = classification_logistic2(train_data_w,train_data_w,label_train,label_train,2)
      output_logistic3 = classification_logistic2(train_data_w,train_data_w,label_train,label_train,3)
      output_rf = classification_rf(train_data_w,train_data_w,label_train,label_train,mtry_w)
      output_xgboost = classification_xgboost(train_data_w,train_data_w,label_train,label_train,max_depth_w,eta_w)
      for(j in 1:4){
        Correct.svm1[i,j+8] = output_svm[[j]]
        Correct.logistic1[i,j+8] = output_logistic[[j]]
        Correct.logistic12[i,j+8] = output_logistic2[[j]]
        Correct.logistic13[i,j+8] = output_logistic3[[j]]
        Correct.rf1[i,j+8] = output_rf[[j]]
        Correct.xgboost1[i,j+8] = output_xgboost[[j]]
      }
      output_svm1 = classification_svm(train_data,train_data,label_train,label_train,gamma1,cost1)
      output_logistic1 = classification_logistic(train_data,train_data,label_train,label_train)
      output_logistic12 = classification_logistic2(train_data,train_data,label_train,label_train,2)
      output_logistic13 = classification_logistic2(train_data,train_data,label_train,label_train,3)
      output_rf1 = classification_rf(train_data,train_data,label_train,label_train,mtry1)
      output_xgboost1 = classification_xgboost(train_data,train_data,label_train,label_train,max_depth1,eta1)
      for(j in 1:4){
        Correct.svm1[i,j+4] = output_svm1[[j]]
        Correct.logistic1[i,j+4] = output_logistic1[[j]]
        Correct.logistic12[i,j+4] = output_logistic12[[j]]
        Correct.logistic13[i,j+4] = output_logistic13[[j]]
        Correct.rf1[i,j+4] = output_rf1[[j]]
        Correct.xgboost1[i,j+4] = output_xgboost1[[j]]
      }
      # Test set results,  9-12
      output_svm = classification_svm(train_data_w,test_data_w,label_train,label_test,gamma_w,cost_w)
      output_logistic = classification_logistic(train_data_w,test_data_w,label_train,label_test)
      output_logistic2 = classification_logistic2(train_data_w,test_data_w,label_train,label_test,2)
      output_logistic3 = classification_logistic2(train_data_w,test_data_w,label_train,label_test,3)
      output_rf = classification_rf(train_data_w,test_data_w,label_train,label_test,mtry_w)
      output_xgboost = classification_xgboost(train_data_w,test_data_w,label_train,label_test,max_depth_w,eta_w)
      for(j in 1:4){
        Correct.svm[i,j+8] = output_svm[[j]]
        Correct.logistic[i,j+8] = output_logistic[[j]]
        Correct.logistic2[i,j+8] = output_logistic2[[j]]
        Correct.logistic3[i,j+8] = output_logistic3[[j]]
        Correct.rf[i,j+8] = output_rf[[j]]
        Correct.xgboost[i,j+8] = output_xgboost[[j]]
      }
      output_svm1 = classification_svm(train_data,test_data,label_train,label_test,gamma1,cost1)
      output_logistic1 = classification_logistic(train_data,test_data,label_train,label_test)
      output_logistic12 = classification_logistic2(train_data,test_data,label_train,label_test,2)
      output_logistic13 = classification_logistic2(train_data,test_data,label_train,label_test,3)
      output_rf1 = classification_rf(train_data,test_data,label_train,label_test,mtry1)
      output_xgboost1 = classification_xgboost(train_data,test_data,label_train,label_test,max_depth1,eta1)
      for(j in 1:4){
        Correct.svm[i,j+4] = output_svm1[[j]]
        Correct.logistic[i,j+4] = output_logistic1[[j]]
        Correct.logistic2[i,j+4] = output_logistic12[[j]]
        Correct.logistic3[i,j+4] = output_logistic13[[j]]
        Correct.rf[i,j+4] = output_rf1[[j]]
        Correct.xgboost[i,j+4] = output_xgboost1[[j]]
      }
      t4=Sys.time()
      print(t4-t3)
    }
    CV.svm1[k,] = apply(Correct.svm1,2,mean,na.rm=TRUE)
    CV.logistic1[k,] = apply(Correct.logistic1,2,mean,na.rm=TRUE)
    CV.logistic12[k,] = apply(Correct.logistic12,2,mean,na.rm=TRUE)
    CV.logistic13[k,] = apply(Correct.logistic13,2,mean,na.rm=TRUE)
    CV.rf1[k,] = apply(Correct.rf1,2,mean,na.rm=TRUE)
    CV.xgboost1[k,] = apply(Correct.xgboost1,2,mean,na.rm=TRUE)
    CV.svm[k,] = apply(Correct.svm,2,mean,na.rm=TRUE)
    CV.logistic[k,] = apply(Correct.logistic,2,mean,na.rm=TRUE)
    CV.logistic2[k,] = apply(Correct.logistic2,2,mean,na.rm=TRUE)
    CV.logistic3[k,] = apply(Correct.logistic3,2,mean,na.rm=TRUE)
    CV.rf[k,] = apply(Correct.rf,2,mean,na.rm=TRUE)
    CV.xgboost[k,] = apply(Correct.xgboost,2,mean,na.rm=TRUE)
    print(k)
  }  
  t2=Sys.time() 
  print(t2-t1)
}  

biaoge = matrix(0,nrow=12,ncol=12)
biaoge[1,] = apply(CV.svm1,2,mean,na.rm=TRUE)
biaoge[2,] = apply(CV.logistic1,2,mean,na.rm=TRUE)
biaoge[3,] = apply(CV.logistic12,2,mean,na.rm=TRUE)
biaoge[4,] = apply(CV.logistic13,2,mean,na.rm=TRUE)
biaoge[5,] = apply(CV.rf1,2,mean,na.rm=TRUE)
biaoge[6,] = apply(CV.xgboost1,2,mean,na.rm=TRUE)
biaoge[7,] = apply(CV.svm,2,mean,na.rm=TRUE)
biaoge[8,] = apply(CV.logistic,2,mean,na.rm=TRUE)
biaoge[9,] = apply(CV.logistic2,2,mean,na.rm=TRUE)
biaoge[10,] = apply(CV.logistic3,2,mean,na.rm=TRUE)
biaoge[11,] = apply(CV.rf,2,mean,na.rm=TRUE)
biaoge[12,] = apply(CV.xgboost,2,mean,na.rm=TRUE)
biaoge

rbind(CV.svm1,CV.logistic1,CV.rf1,CV.xgboost1,
      CV.svm,CV.logistic,CV.rf,CV.xgboost)

# setwd("C:/Users/123456/Desktop")
save.image(file = "simulation.RData")
load("simulation.RData")