# Add experiments:
#   1. The characteristics of different stages of colorectal cancer (CRC) patients;
#   2. Distinguish between adenoma (precancerous stage) and normal stages. 
#   3. Two solutions to address the issue of redundancy

############# Functions ########### 
# The following functions are detailed in model_functions.R
divide_train_test <- function(s_matrix_HC,s_matrix_CRC,r_len,k_i) {
  #s_matrix_HC ,s_matrix_CRC
  #r_len<1表示随机选择测试集，反之代表交叉验证，k_i表示交叉验证次数
  m = dim(s_matrix_HC)[2]  # 200
  n0 = dim(s_matrix_HC)[1] # 524
  n1 = dim(s_matrix_CRC)[1] # 503
  s_matrix11 = rbind(s_matrix_HC,s_matrix_CRC)  # 1027*200
  lable_s =c(rep(0, n0),rep(1, n1))   #label
  # s_matrix11 = as.data.frame(t(s_matrix11))  # 200*1027
  # s_matrix1 <- as.data.frame(lapply(s_matrix11,gyh))
  s_matrix1 <- as.data.frame(s_matrix11) 
  if(r_len<1){
    c1 = sample(c(1:n0),size=n0*r_len)
    c2 = sample(c((n0+1):(n0+n1)),size=n1*r_len)
  }else{
    c1 = floor(n0/r_len*(k_i-1)+1):floor(n0/r_len*k_i)
    c2 = floor(n0+n1/r_len*(k_i-1)+1):floor(n0+n1/r_len*k_i)
  }
  test_data = s_matrix1[c(c1,c2),];   p_test = lable_s[c(c1,c2)]
  train_data = s_matrix1[-c(c1,c2),]; p_train = lable_s[-c(c1,c2)]
  
  # 乱序一下
  n_test_data = sample(nrow(test_data))
  n_train_data = sample(nrow(train_data))
  test_data <- test_data[n_test_data,]; p_test =p_test[n_test_data]
  train_data <- train_data[n_train_data,]; p_train =p_train[n_train_data]
  
  output_divide = list(train_data,test_data,p_train,p_test)
  rm(list=c('m','n0','n1','s_matrix11','lable_s','s_matrix1','c1','c2','test_data','train_data'))
  gc()
  return(output_divide)
} # HC在前

CV_variable_selection <- function(s_matrix_HC,s_matrix_CRC,r_len){
  n_sample_HC = sample(nrow(s_matrix_HC))
  n_sample_CRC = sample(nrow(s_matrix_CRC))
  s_matrix_HC <- s_matrix_HC[n_sample_HC,]
  s_matrix_CRC <- s_matrix_CRC[n_sample_CRC,]
  n_select1 = c() 
  n_select2 = c()
  for(i in 1:r_len){
    output_divide = divide_train_test(s_matrix_HC,s_matrix_CRC,r_len,i)
    train_data = output_divide[[1]]    
    test_data = output_divide[[2]]    
    label_train = output_divide[[3]]  
    label_test = output_divide[[4]]
    
    #
    output_chi = contingency_table_test(train_data,test_data,label_train,0.05,10)
    train_data = output_chi[[1]]   
    test_data =output_chi[[2]]      
    n_select = output_chi[[3]]   
    n_select10 = output_chi[[5]]  
    n_select1 = c(n_select1,n_select)
    n_select2 = c(n_select2,n_select10)
  }
  cv.select = names(table(n_select1)[table(n_select1)>=r_len*0.6])
  cv.select = as.numeric(cv.select)
  s_matrix_HC1 = s_matrix_HC[,cv.select]
  s_matrix_CRC1 = s_matrix_CRC[,cv.select]
  output_cvvs = list(s_matrix_HC1,s_matrix_CRC1,cv.select,table(n_select1),n_select2)
  # rm(list=c('train_data','test_data','train_data1','test_data1'))
  rm(list=c('train_data','test_data'))
  gc()
  return(output_cvvs)
}

contingency_table_test <- function(train_data,test_data,label_train,pa,n){
  s_lable_01 = data.frame(train_data)  # ת??Ϊ??��Ϊ??
  # s_lable_01 = ifelse(s_lable_01>0,1,0)
  # s_lable_01 = ifelse(s_lable_01>=0,1,0)
  s_lable_01[s_lable_01>0] = 1
  s_lable_01[s_lable_01!=1] = 0   # ??Ϊ?????????еĸ???
  size = dim(s_lable_01)
  a<-seq(0,0,length=size[2])
  for(i in 1:size[2]){
    if(sum(s_lable_01[[i]])>0){
      if(min(table(s_lable_01[[i]],label_train))<=5){
        p = fisher.test(table(s_lable_01[[i]],label_train))$p.value
        a[i] = p
      }else{
        p = chisq.test(table(s_lable_01[[i]],label_train))$p.value
        a[i] = p
      }
    }else{
      a[i] = 1
    }
  }
  s_p_n = order(a)[1:n]      # ????????
  s_len = length(a[a<=pa])   # Number of selected microorganisms
  #Select samples that have passed the microbiological test;
  s_train=apply(as.data.frame(t(train_data)),2,as.numeric)  #??��Ϊ?? 100*800
  s_test=apply(as.data.frame(t(test_data)),2,as.numeric)
  train_data = as.data.frame(t(subset(s_train,a<=pa)))    # ???н???ɸѡ(???л?????��չ)
  test_data = as.data.frame(t(subset(s_test,a<=pa))) # ??��˳??Ӧ??û?䣿
  train_data_n = as.data.frame(t(s_train[sort(s_p_n),]))  # ??��˳?????ò??䣿
  test_data_n = as.data.frame(t(s_test[sort(s_p_n),]))
  # save(s_matrix,file="s_matrix.Rdata")
  # load("s_matrix.Rdata")  # ????ɸѡ????????200*1027
  s_p_order = c(1:size[2])[a<=pa]   # ????ɸѡ??????΢??????????
  # ??????????��Ϊ??
  output_chi = list(train_data,test_data,s_p_order,s_len,s_p_n,a,train_data_n,test_data_n)
  rm(list=c('size','train_data','test_data','a','train_data_n','test_data_n'))
  gc()
  return(output_chi)
}

binary_matrix_conversion32 <- function(train_data,test_data,label_train,label_test,n_cor){
  # train_data 变量为列
  train_data1 <- train_data[,-c(1,2,3)]
  test_data1 <- test_data[,-c(1,2,3)]
  m = dim(test_data1)[2]   # 200
  p_pre = matrix(0,dim(test_data1)[1],dim(test_data1)[2] )
  t_pre = matrix(0,dim(train_data1)[1],dim(train_data1)[2] )
  cor_s = cor(train_data1,method = "spearman")  # Correlation matrix
  for(j in 1:m){
    sy = order(abs(cor_s[j,]),decreasing=TRUE)[1:n_cor]   #  TRUE
    # sy = order(abs(cor_s[j,]),decreasing=TRUE)[1:5]   #  TRUE
    #sy = order(cor_s[j,],decreasing=TRUE)[1:5]   #  TRUE
    test_data11 = test_data1[,c(sy)]             # 11*n/k
    train_data11 = train_data1[,c(sy)]
    # test_data11 = cbind(test_data[,c(1,2,3)],test_data11,test_data11[,1]*test_data11[,-1])
    # train_data11 = cbind(train_data[,c(1,2,3)],train_data11,train_data11[,1]*train_data11[,-1])
    # test_data11 = cbind(test_data[,c(1,2,3)],test_data11,test_data11[,1]*test_data11[,-1],test_data11[,1]*test_data[,c(1,2,3)])
    # train_data11 = cbind(train_data[,c(1,2,3)],train_data11,train_data11[,1]*train_data11[,-1],train_data11[,1]*train_data[,c(1,2,3)])
    test_data11 = cbind(test_data11,test_data11[,1]*test_data11[,-1])
    train_data11 = cbind(train_data11,train_data11[,1]*train_data11[,-1])
    # pfit = glmnet(train_data11,label_train, family = "binomial")  #0.027s
    # s_lambda = min(pfit$lambda)
    # pred1 = predict(pfit,train_data11,s = s_lambda , type = "class") 
    test_data11 = as.matrix(test_data11)
    train_data11 = as.matrix(train_data11)
    pfit = cv.glmnet(train_data11,label_train, family = "binomial",type.measure ="class",nfolds = 5)  # 0.154
    #pfit$lambda.min
    # coefficients<-coef(pfit, s = "lambda.min")  #如下是变量的系数
    # Active.Index<-which(coefficients!=0) #系数不为0的特征索引
    # Active.coefficients<-coefficients[Active.Index] #系数不为0的特征系数值
    pred1 = predict(pfit,train_data11,s = 'lambda.min',type = "class") 
    # plot(pred) # type = "response"
    pred = predict(pfit,test_data11,s = 'lambda.min',type = "class") 
    t_pre[,j] <- ifelse(pred1==1,1,-1)
    p_pre[,j] <- ifelse(pred==1,1,-1)
  }
  #dim(test_data11)
  test_data = p_pre  # 200*n/k
  train_data = t_pre  #200*(n-n/k)
  output_bi = list(train_data,test_data,label_train,label_test)
  rm(list=c('train_data','test_data','p_pre','t_pre'))
  gc()
  return(output_bi)
}

if(TRUE){
  
  classification_svm <- function(train_data,test_data,label_train,label_test,gamma1,cost1){
    # train_data = train_data_w
    # test_data = train_data_w
    # label_train = label_train
    # label_test = label_train
    # gamma1 = 0.0001
    # cost1 = 1
    #rownames(train_data) = paste('', c(1:dim(train_data)[1]),sep = "", collapse = NULL)
    #rownames(test_data) = paste('', c(1:dim(test_data)[1]),sep = "", collapse = NULL)
    d_train <- data.frame(train_data,y=as.factor(label_train))
    d_test <- data.frame(test_data,y=as.factor(label_test))
    #tuned <- tune.svm(y~.,data=d_train,gamma=10^(-4:1),cost=10^(-2:2))
    # $best.parameters
    # gamma cost
    # 22   0.1   10
    # gamma越小第一个值可能越大？
    train_svm = svm(y~.,data=d_train,kernel='radial',gamma=gamma1,cost=cost1,scale=F)
    #pred_train=predict(train_svm,d_train)
    #table(pred=pred_train,truth=label_train)
    pred=predict(train_svm,d_test,decision.values = TRUE)    # testing set
    # pred=predict(train_svm,d_train,decision.values = TRUE)  # training set
    #lebal_test = lebal_train
    if(length(pred)!=length(label_test)){
      print('length(pred)!=length(label_test)')
      label_test = label_test[c(1:length(label_test)) %in% names(pred)]
      label_test = label_test[!is.na(pred)]
      pred = pred[!is.na(pred)]
    }
    test_table = table(pred=pred,truth=label_test)
    pred=predict(train_svm,d_test,decision.values = TRUE)    # testing set
    auc = classification_auc(-attr(pred, "decision.values"),label_test)
    auc = ifelse(auc<0.5,1-auc,auc)
    #pred CRC HC
    TP = test_table[1]
    FP = test_table[2]
    FN = test_table[3]
    TN = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_svm = list(TPR,FPR,ACC,auc)
    rm(list=c('d_train','d_test','train_svm','pred','test_table'))
    gc()
    return(output_svm)
  }
  
  classification_svm1 <- function(train_data,test_data,label_train,label_test,sigma1,C1){
    #colnames(train_data) = paste('', c(1:dim(train_data)[2]),sep = "", collapse = NULL)
    #colnames(test_data) = paste('', c(1:dim(test_data)[2]),sep = "", collapse = NULL)
    d_train <- data.frame(train_data,y=as.factor(label_train))
    d_test <- data.frame(test_data,y=as.factor(label_test))
    train_svm <- ksvm(y~., data = d_train,type = "C-svc", kernel = "rbfdot",
                      kpar=  list(sigma=sigma1), # "automatic",#  # #
                      C=C1,cross=3,prob.model = TRUE) #训练  ,cross=5
    # train_svm <- ksvm(y~., data = d_train) #训练  ,cross=5
    # train_svm <- ksvm(y~., data = d_train, kernel = "rbfdot",
    #                   kpar=  list(sigma=sigma1), # "automatic",#  # #
    #                   C=C1,cross=3,prob.model = TRUE) #训练  ,cross=5
    #pred_train=predict(train_svm,d_train)
    #table(pred=pred_train,truth=label_train)
    pred = predict(train_svm,d_train,type = "response")#    # testing set
    #pred=predict(train_svm,d_train,decision.values = TRUE)  # training set
    #lebal_test = lebal_train
    if(length(pred)!=length(label_test)){
      print('length(pred)!=length(label_test)')
      #lebal_test = lebal_test[c(1:length(label_test)) %in% names(pred)]
      label_test = label_test[!is.na(pred)]
      pred = pred[!is.na(pred)]
    }
    test_table = table(pred=pred,truth=label_test)
    pred=predict(train_svm,d_test,type = "probabilities")    # testing set
    auc = classification_auc(pred[,1],label_test)
    #pred=predict(train_svm,d_test,type = "decision")
    #auc = classification_auc(pred,lebal_test)
    #pred CRC HC
    TP = test_table[1]
    FP = test_table[2]
    FN = test_table[3]
    TN = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_svm = list(TPR,FPR,ACC,auc)
    rm(list=c('d_train','d_test','train_svm','pred','test_table'))
    gc()
    return(output_svm)
  } 
  
  classification_logistic1 <- function(train_data,test_data,label_train,label_test){
    #label_train = ifelse(label_train=='HC',0,1)
    #label_test = ifelse(label_test=='HC',0,1)
    d_train <- data.frame(train_data,y=as.factor(label_train))
    d_test <- data.frame(test_data,y=as.factor(label_test))
    glm.fit = glm(y~.,family = binomial,d_train)   # LR
    #summary(glm.fit)
    glm.pre1 = predict(glm.fit,d_train,type="response") #“link”, “response”, “terms”
    glm.pre1=1/(1+exp(glm.pre1)) 
    #plot(sort(glm.pre1))
    #plot(glm.pre1)
    #table(pred=pred_train,truth=label_train)
    ii = 0.01
    errorq = seq(0,0,length=1/ii)
    for(q in seq(ii,1,ii)){
      qq = seq(0,0,length=length(glm.pre1))
      qq[glm.pre1 < q] = 0
      qq[glm.pre1!=0] = 1
      errorq[q/ii] = length(label_train[qq==label_train])/length(glm.pre1)
    }
    p0 = median(seq(ii,1,ii)[errorq==max(errorq)])
    
    pred = predict(glm.fit,d_test,type = 'response')
    pred=1/(1+exp(pred)) 
    auc = classification_auc(pred,label_test)
    
    #plot(pred)
    p0 = 1/(1+exp(p0)) 
    pred = ifelse(pred < p0,0,1)
    
    test_table = table(pred=pred,truth=label_test)
    # auc = classification_auc(pred,label_test)
    TP = test_table[1]
    FP = test_table[2]
    FN = test_table[3]
    TN = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_logistic = list(TPR,FPR,ACC,auc)
    rm(list=c('d_train','d_test','glm.fit','pred','test_table'))
    gc()
    return(output_logistic)
  }
  
  classification_logistic <- function(train_data,test_data,label_train,label_test){
    train_data = as.matrix(train_data)
    test_data = as.matrix(test_data)
    cvfit = cv.glmnet(train_data,label_train,family = "binomial",type.measure ="class")     # "class" "auc"
    glm.pre1 = predict(cvfit, train_data, s = "lambda.min", type = "class")
    pred = predict(cvfit,test_data, s = "lambda.min", type = "response")
    auc = classification_auc(1-pred,label_test)
    auc = ifelse(auc<0.5,1-auc,auc)
    pred = predict(cvfit,test_data, s = "lambda.min", type = "class")
    test_table = table(pred=pred,truth=label_test)
    TP = test_table[1]
    FP = test_table[2]
    FN = test_table[3]
    TN = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_logistic = list(TPR,FPR,ACC,auc)
    rm(list=c('cvfit','glm.pre1','pred','test_table'))
    gc()
    return(output_logistic)
  }
  
  classification_logistic2 <- function(train_data,test_data,label_train,label_test,interactions){
    train_data <- model_interactions(train_data,0.1,interactions)
    test_data <- model_interactions(test_data,0.1,interactions)
    train_data = as.matrix(train_data)
    test_data = as.matrix(test_data)
    
    cvfit = cv.glmnet(train_data,label_train,family = "binomial",type.measure ="class")     # "class" "auc"
    glm.pre1 = predict(cvfit, train_data, s = "lambda.min", type = "class")
    pred = predict(cvfit,test_data, s = "lambda.min", type = "response")
    auc = classification_auc(1-pred,label_test)
    auc = ifelse(auc<0.5,1-auc,auc)
    pred = predict(cvfit,test_data, s = "lambda.min", type = "class")
    test_table = table(pred=pred,truth=label_test)
    TP = test_table[1]
    FP = test_table[2]
    FN = test_table[3]
    TN = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_logistic = list(TPR,FPR,ACC,auc)
    rm(list=c('cvfit','glm.pre1','pred','test_table'))
    gc()
    return(output_logistic)
  }
  
  classification_tree <- function(train_data,test_data,label_train,label_test){
    d_train <- data.frame(train_data,y=as.factor(label_train))
    d_test <- data.frame(test_data,y=as.factor(label_test))
    train_tree = tree(y~.,data=d_train)
    pred_train=predict(train_tree,d_train,type = 'class')
    table(pred=pred_train,truth=label_train)
    pred = predict(train_tree,d_test,type = 'class')
    test_table = table(pred=pred,truth=label_test)
    pred_auc =predict(train_tree,d_test,type = 'vector')
    auc = classification_auc(pred_auc[,1],label_test)
    TP = test_table[1]
    FP = test_table[2]
    FN = test_table[3]
    TN = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_tree = list(TPR,FPR,ACC,auc)
    rm(list=c('d_train','d_test','train_tree','pred','test_table'))
    gc()
    return(output_tree)
  }
  
  classification_rf <- function(train_data,test_data,label_train,label_test,mtry1){
    rownames(train_data) = paste('', c(1:dim(train_data)[1]),sep = "", collapse = NULL)
    rownames(test_data) = paste('', c(1:dim(test_data)[1]),sep = "", collapse = NULL)
    d_train <- data.frame(train_data,y=as.factor(label_train))
    d_test <- data.frame(test_data,y=as.factor(label_test))
    # mtry=(length(d_test))^(1/2),
    train_rf = randomForest(y~.,data=d_train,mtry=mtry1,na.action = na.omit)
    #pred_train = predict(train_rf,d_train,type = 'class')
    #table(pred=pred_train,truth=lebal_train)
    #d_test = d_train               # training set
    #lebal_test = lebal_train        # training set
    pred=predict(train_rf,d_test,type = 'class')    # testing set
    pred=na.omit(pred)
    if(length(pred)!=length(label_test)){
      print('length(pred)!=length(lebal_test)')
      label_test=label_test[c(1:length(label_test)) %in% names(pred)]  # quchu NA
    }
    test_table = table(pred=pred,truth=label_test)
    pred_auc = predict(train_rf,d_test,type="prob",)
    pred_auc = na.omit(pred_auc)
    #lebal_test2=lebal_test[c(1:length(label_test)) %in% row.names(pred_auc)]  # 其实去一次就行
    auc = classification_auc(pred_auc[,1],label_test)
    auc = ifelse(auc<0.5,1-auc,auc)
    TP = test_table[1]
    FP = test_table[2]
    FN = test_table[3]
    TN = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_rf = list(TPR,FPR,ACC,auc)
    rm(list=c('d_train','d_test','train_rf','pred','test_table'))
    gc()
    return(output_rf)
  }
  
  classification_xgboost <- function(train_data,test_data,label_train,label_test,max_depth1,eta1){
    # label_test[label_test=='HC']=0
    # label_test[label_test=='CRC']=1
    # label_train[label_train=='HC']=0
    # label_train[label_train=='CRC']=1
    train_data = as.matrix(train_data)
    test_data = as.matrix(test_data)
    dtrain <- xgb.DMatrix(data=train_data,label=label_train) 
    dtest <- xgb.DMatrix(data=test_data,label=label_test)
    # max_depth=2,eta=0.1
    train_xgboost <- xgboost(data=dtrain,max_depth=max_depth1,eta=eta1,
                             objective='binary:logistic',nround=300,verbose=0)
    #在测试集上预测
    pred_xgb = predict(train_xgboost,dtest)
    auc = classification_auc(pred_xgb,label_test)
    pred_xgb1 = round(pred_xgb)
    test_table = table(pred_xgb1,label_test)
    # pred_xgb1  0  1
    TN = test_table[1]
    FN = test_table[2]
    FP = test_table[3]
    TP = test_table[4]
    TPR =  TP/(TP+FN)
    FPR =  FP/(FP+TN)
    ACC =  (TP+TN)/(TP+FP+FN+TN)
    output_xgboost = list(TPR,FPR,ACC,auc)
    rm(list=c('dtrain','dtest','train_xgboost','pred_xgb','test_table'))
    gc()
    return(output_xgboost)
  }
  
} 

classification_auc <- function(pred_test,label_test){
  # library(ROCR)  # 计算AUC值
  pred_test = as.matrix(pred_test) # 预测结果
  # pred_test[pred_test=='CRC']=1
  # pred_test[pred_test=='HC']=0
  pred_test = as.numeric(pred_test)
  label_test = as.matrix(label_test) # 真实标签
  # label_test[label_test=='CRC']=1
  # label_test[label_test=='HC']=0
  label_test = as.numeric(label_test)
  
  pred = prediction(pred_test,label_test)
  perf = performance(pred,"tpr","fpr")
  #plot(perf,col="blue",lty=3,lwd=3,cex.lab=1.5,cex.axis=2,cex.main=1.5,main="ROC plot")
  auc = performance(pred,"auc")
  auc = unlist(slot(auc, "y.values"))
  return(auc)
}

vartest_CRC2 <- function(data,label,n_cv){
  # Description: Test criteria for variable heterogeneity (heterogeneity test).
  # Input:
  #   data: Initial covariate data.
  #   label: Corresponding labels.
  #   n_cv: Number of observation groups.
  # Output:
  #   matrix_coef: Model fitting results.
  #   n_coef: The count that satisfies the test.
  #   n_coef/m: The percentage that satisfies the test.
  #   n_coef_mean: The number of covariates tested by different groups converged
  #   n_byz: The ordinal number of the covariate that satisfies the test.
  
  m = dim(data)[2]   
  n1 = dim(data)[1]  
  n_for = 10
  n_coef_mean = seq(0,length.out = n_for)
  n_byz = c()
  for(k in 1:10){
    order  =  sample(n1)
    data <- data[order,]
    label <- label[order]
    matrix_coef = matrix(0,nrow = n_cv,ncol = m+1 ) 
    for(i in 1:n_cv){
      train_data = as.matrix(data)
      c1 = floor(n1/n_cv*(i-1)+1):floor(n1/n_cv*i)
      train_data = train_data[c1,]
      label_train =label[c1]
      cvfit = cv.glmnet(train_data,label_train,family = "binomial",type.measure ="class") 
      matrix_coef[i,] = coef(cvfit, s = "lambda.min")[,1]
    }
    matrix_coef = matrix_coef[,-1]
    n_coef = 0
    for(i in 1:dim(matrix_coef)[2]){
      z = length(matrix_coef[,i][matrix_coef[,i]>0])  
      f = length(matrix_coef[,i][matrix_coef[,i]<=0]) 
      if(z>0 & f>0){ 
        n_coef = n_coef + 1
        n_byz = c(n_byz,i)
      }
    }
    n_coef_mean[k] = n_coef
  }
  n_coef = sum(n_coef_mean)/n_for
  output_vartest = list(matrix_coef,n_coef,n_coef/m,n_coef_mean,n_byz)
  return(output_vartest)
}

gyh<-function(x){
  if(length(table(x))!=1){  # sum(x)!=0，因为只有这种情况
    # x = (x-min(x))/(max(x)-min(x))  # 0 1
    mean_x = (max(x)+min(x))/2
    x = (x-mean_x)*3/(max(x)-min(x))   # -1 1 *2
  }
  return(x)
}

gyh1<-function(x){
  if(length(table(x))!=1){  # sum(x)!=0，因为只有这种情况
    x = (x-min(x))/(max(x)-min(x))  # 0 1
    # mean_x = (max(x)+min(x))/2
    # x = (x-mean_x)*3/(max(x)-min(x))   # -1 1 *2
  }
  return(x)
}

arules_data01_label <- function(test_data_01,label_test,vs_select,support,confidence){
  
  # Description: Generation and filtering of association rules.
  # Input:
  #   test_data_01: The binary matrix corresponding to the original data.
  #   label_test: Corresponding labels.
  #   vs_select: Results of covariate selection.
  #   support, confidence: Initial screening parameters of association rules.
  # Output:
  #   rules1_js: Association rules related to CRC or 1.
  #   rules2_js: Association rules related to HC or 0.
  
  # library(ROCR)  # Calculating AUC value
  # test_data_01_po = test_data_01   # positive
  # test_data_01_ne = test_data_01  # negative
  if(length(label_test[label_test=='CRC'])>0){
    label_test1 = ifelse(label_test=='CRC',1,0)  # CRC
    label_test2 = ifelse(label_test=='HC',1,0)   # HC
  }else{
    label_test1 = ifelse(label_test==1,1,0)  # CRC
    label_test2 = ifelse(label_test==0,1,0)   # HC
  }
  # label_test1 = ifelse(label_test=='CRC',1,0)   # CRC
  # label_test2 = ifelse(label_test=='HC',1,0)    # HC
  n_select = vs_select                            # output_cvvs[[3]]
  test_data_01_po = ifelse(test_data_01==1,1,0)   # positive
  colnames(test_data_01_po) = paste('po',n_select,sep = "", collapse = NULL)
  test_data_01_ne = ifelse(test_data_01==-1,1,0)  # negative
  colnames(test_data_01_ne) = paste('ne',n_select,sep = "", collapse = NULL)
  
  ######## Screening CRC-related association rules.
  data_CRC_rules = cbind(test_data_01_po,test_data_01_ne,label_test1)
  rules1 <- apriori(data_CRC_rules,parameter=list(supp=support,conf=confidence,target = "rules" ),
                    control = list ( verbose = FALSE ),
                    appearance = list(rhs ="label_test1"))  #list(items ="label_test1") 
  rules1.sorted <- sort(rules1, by = "confidence")  # by = "lift"
  quality(rules1.sorted)$improvement <- interestMeasure(rules1.sorted, measure = "improvement")
  # inspect(rules1.sorted)
  # is.redundant(rules1.sorted)
  ## non-redundant rules
  rules1.pruned <- rules1.sorted[!is.redundant(rules1.sorted)]
  # inspect(subset(rules1.pruned, items %in% 'label_test1' ))   
  
  rules1.pruned_label = inspect(sort(rules1.pruned,by="lift"))
  rules1.pruned_label$len_lhs  = sapply(rules1.pruned_label$lhs,n_var)
  rules1.pruned_label$len_ne  = sapply(rules1.pruned_label$lhs,n_ne)   
  #rules1_results = rules1.pruned_label[order(rules1.pruned_label$len_lhs,decreasing=T),]
  rules1_results = rules1.pruned_label[order(rules1.pruned_label$len_ne,decreasing=T),]
  # rules1_results[rules1_results$rhs=="{label_test1}",]
  # rules1_js = rules1_results[rules1_results$rhs=="{label_test1}",]   
  rules1_js = rules1_results
  
  ######## Screening HC-related association rules.
  data_HC_rules = cbind(test_data_01_po,test_data_01_ne,label_test2)
  rules2 <- apriori(data_HC_rules,parameter=list(supp=support,conf=confidence,target = "rules" ),
                    control = list(verbose = FALSE),
                    appearance = list(rhs ="label_test2"))
  rules2.sorted <- sort(rules2,by = "confidence")  # by = "lift" ,by = "confidence"
  quality(rules2.sorted)$improvement <- interestMeasure(rules2.sorted, measure = "improvement")
  rules2.pruned <- rules2.sorted[!is.redundant(rules2.sorted)]
  # rules2.pruned@rhs@itemInfo
  # items
  # rules2.pruned_label = inspect(sort(subset(rules2.pruned, items %in% 'label_test2' ),by="lift"))
  rules2.pruned_label = inspect(sort(rules2.pruned,by="lift"))
  rules2.pruned_label$len_lhs  = sapply(rules2.pruned_label$lhs,n_var)
  rules2.pruned_label$len_ne  = sapply(rules2.pruned_label$lhs,n_ne)
  #rules2_results = rules2.pruned_label[order(rules2.pruned_label$len_lhs,decreasing=T),]
  rules2_results = rules2.pruned_label[order(rules2.pruned_label$len_ne,decreasing=T),]
  # rules2_results[rules2_results$rhs=="{label_test2}",]
  # rules2_js = rules2_results[rules2_results$rhs=="{label_test2}",]  
  rules2_js = rules2_results
  output_rules = list(rules1_js,rules2_js)
  rm(list=c('rules1','rules2','rules1.pruned_label','rules2.pruned_label'))
  gc()
  return(output_rules)
}

n_var <- function(a){ 
  # Description: Filter the number of covariates to the left of the association rule.
  # Input:
  #   a: The expression to the left of the association rule.
  # Output:
  #   n_var: the number of covariates to the left of the association rule.
  a = strsplit(a,',')
  n_var = length(a[[1]])
  return(n_var)
}

n_ne <-function(a){
  # Description: To screen for covariates with negative effects.
  # Input:
  #   a: The expression to the left of the association rule.
  # Output:
  #   n_var: the number of covariates with negative effects.
  a = strsplit(a,'ne')
  n_ne = length(a[[1]])-1
  return(n_ne)
}

library(gbm)
library('e1071')  #SVM
library("kernlab")
library(tree)
library(randomForest)
library(xgboost)
# library(Matrix)   
library(ROCR)     
library(EnvStats) 
library(mvtnorm)   
# library(GenSA)    # 模拟退火
library(glmnet)
library(mniw)   
library(arules)
library(arulesViz)


# The adenoma (precancerous stage) and normal stage of CRC were studied
if(FALSE){
  # setwd("D:/OneDrive - tshaha/研究生工作/肠道微生物分类/Microbial_dataa")
  s = read.csv("C:/Users/duorou/Desktop/Microbial_data/s_label.csv")
  g = read.table("C:/Users/duorou/Desktop/Microbial_data/group_1224_112.txt",header=T,na.strings = c("NA"))
  bmi <- read.csv("C:/Users/duorou/Desktop/Microbial_data/BMI1.csv")
  bmi <- unique(bmi)
  g <- merge(bmi,g,by="ID")  
  s_lable_ID <- merge(g,s,by="ID") 
  s_lable_ID$Gender <- ifelse(s_lable_ID$Gender == "M",0,1)
  s_lable_ID <- subset(s_lable_ID,s_lable_ID$BMI > 0 )
  # table(s_lable_ID$CRCorHC,s_lable_ID$COUNTRY) 数据不平衡
  # table(s_lable_ID$CRCorHC,s_lable_ID$sex)
  # table(s_lable_ID$Country,s_lable_ID$CRCorHC)
  # table(s_lable_ID$CRCorHC)
  
  # Study of CRC adenomas (precancerous stage), using French data for analysis; 59+27 CRC patients;
  if(FALSE){
    table(s_lable_ID$AJCC_stage)
    s_lable_ID_CRC <- subset(s_lable_ID,
                             s_lable_ID$AJCC_stage %in% c("Normal","Small adenoma"))
    g_label = s_lable_ID_CRC$Disease.Stage  
    s_lable01_CRC = s_lable_ID_CRC[,-c(1,5:13)]   #-c(1,5,6,7,8) 1224  845
    s_lable01_CRC[,c(1,2,3)] = apply(s_lable01_CRC[,c(1,2,3)],2,gyh)
    dim(s_lable01_CRC)   #[1] 1224  846; 86 848
    s_matrix_HC = subset(s_lable01_CRC,g_label =='Normal')  # 619*845
    s_matrix_CRC = subset(s_lable01_CRC,g_label =='Small adenoma') # 605*845
    
  }
  
  if(TRUE){
    g_label = s_lable_ID$CRCorHC  
    s_lable01 = s_lable_ID[,-c(1,5:13)]   #-c(1,5,6,7,8) 1224  845
    s_lable01[,c(1,2,3)] = apply(s_lable01[,c(1,2,3)],2,gyh)
    # write.csv(s_lable01,file = "s_lable01.csv",row.names = F)
    dim(s_lable01)   #[1] 1224  846
    s_matrix_HC = subset(s_lable01,g_label =='HC')  # 619*845
    s_matrix_CRC = subset(s_lable01,g_label =='CRC') # 605*845
  }

}

# Article theme experiment. And dealing with the classification of different stages of CRC.
if(FALSE){
  # setwd("D:/OneDrive - tshaha/研究生工作/肠道微生物分类/Microbial_dataa")
  s = read.csv("C:/Users/duorou/Desktop/Microbial_data/s_label.csv")
  g = read.table("C:/Users/duorou/Desktop/Microbial_data/group_1224_112.txt",header=T,na.strings = c("NA"))
  bmi <- read.csv("C:/Users/duorou/Desktop/Microbial_data/BMI2.csv")
  bmi_id <- read.csv("C:/Users/duorou/Desktop/Microbial_data/BMI_id.csv")
  bmi_id <- unique(bmi_id)
  colnames(bmi_id)
  bmi <- merge(bmi_id,bmi,by="SampleID")  
  
  g <- merge(bmi,g,by="ID")  
  s_lable_ID <- merge(g,s,by="ID") 
  s_lable_ID <- subset(s_lable_ID, s_lable_ID$Gender %in% c("M","F"))
  
  s_lable_ID$Gender <- ifelse(s_lable_ID$Gender == "M",0,1)
  s_lable_ID <- subset(s_lable_ID,s_lable_ID$BMI > 0 )
  # table(s_lable_ID$COUNTRY.x,s_lable_ID$CRCorHC)
  # table(s_lable_ID$CRCorHC)
  # s_lable_ID[is.na(s_lable_ID)==TRUE]
  
  # The early and Advanced outcomes of CRC were studied. 459 CRC patients;
  if(TRUE){
    s_lable_ID_CRC <- subset(s_lable_ID,
                             s_lable_ID$Disease.Stage %in% c("Advanced","Early"))
    g_label = s_lable_ID_CRC$Disease.Stage  
    s_lable01_CRC = s_lable_ID_CRC[,-c(1:11,15:18)]   #-c(1,5,6,7,8) 1224  845
    s_lable01_CRC[,c(1,2,3)] = apply(s_lable01_CRC[,c(1,2,3)],2,gyh)
    dim(s_lable01)   #[1] 1224  846
    s_matrix_HC = subset(s_lable01_CRC,g_label =='Early')  # 619*845
    s_matrix_CRC = subset(s_lable01_CRC,g_label =='Advanced') # 605*845
    g_label <- ifelse(g_label=='Early',"HC","CRC")
  }
  
  
  if(FALSE){
    g_label = s_lable_ID$CRCorHC  
    s_lable01 = s_lable_ID[,-c(1:11,15:18)]   #-c(1,5,6,7,8) 1224  845
    s_lable01[,c(1,2,3)] = apply(s_lable01[,c(1,2,3)],2,gyh)
    
    # write.csv(s_lable01,file = "s_lable01.csv",row.names = F)
    dim(s_lable01)   #[1] 1224  846
    s_matrix_HC = subset(s_lable01,g_label =='HC')  # 619*845
    s_matrix_CRC = subset(s_lable01,g_label =='CRC') # 605*845
  }
}

# The solution one to the redundancy problem of association rules;
# The BMI2 data were manipulated, but at the genus level, the columns were grouped and summed
if(TRUE){
  # setwd("D:/OneDrive - tshaha/研究生工作/肠道微生物分类/Microbial_dataa")
  s = read.csv("C:/Users/duorou/Desktop/Microbial_data/s_label.csv")
  g = read.table("C:/Users/duorou/Desktop/Microbial_data/group_1224_112.txt",header=T,na.strings = c("NA"))
  bmi <- read.csv("C:/Users/duorou/Desktop/Microbial_data/BMI2.csv")
  bmi_id <- read.csv("C:/Users/duorou/Desktop/Microbial_data/BMI_id.csv")
  bmi_id <- unique(bmi_id)
  colnames(bmi_id)
  bmi <- merge(bmi_id,bmi,by="SampleID")  
  
  g <- merge(bmi,g,by="ID")  
  s_lable_ID <- merge(g,s,by="ID") 
  s_lable_ID <- subset(s_lable_ID, s_lable_ID$Gender %in% c("M","F"))
  
  s_lable_ID$Gender <- ifelse(s_lable_ID$Gender == "M",0,1)
  s_lable_ID <- subset(s_lable_ID,s_lable_ID$BMI > 0 )
  # table(s_lable_ID$COUNTRY.x,s_lable_ID$CRCorHC)
  # table(s_lable_ID$CRCorHC)
  # s_lable_ID[is.na(s_lable_ID)==TRUE]
  g_label = s_lable_ID$CRCorHC  
  s_lable01_Genus = s_lable_ID[,-c(1:18)]   #-c(1,5,6,7,8) 1224  845
  group_gm1 = read.csv("C:/Users/duorou/Desktop/Microbial_data/group_gm1.csv") 
  b <- c()
  Genus_names <- names(table(group_gm1$Genus))
  for(i in 1:length(Genus_names)){ 
    number <- group_gm1$number.y[group_gm1$Genus == Genus_names[i]]
    if(length(number)>1){
      a <- apply(s_lable01_Genus[,number],1,sum)
      b <- cbind(b,a)
    }else{
      a <- s_lable01_Genus[,number] 
      b <- cbind(b,a)
    }
  } 
  s_lable01_Genus <- as.data.frame(cbind(s_lable_ID$Gender,s_lable_ID$Age,s_lable_ID$BMI,b))
  colnames(s_lable01_Genus) <- c("Gender","Age","BMI",Genus_names)
  s_lable01 <- s_lable01_Genus
  s_lable01[,c(1,2,3)] = apply(s_lable01[,c(1,2,3)],2,gyh)
  s_matrix_HC = subset(s_lable01,g_label =='HC')  # 619*845
  s_matrix_CRC = subset(s_lable01,g_label =='CRC') # 605*845
}

#####  Variable selection  #################### 

# Perform 100 variable selections in parallel
if(TRUE){
  library(foreach)
  library(parallel)
  library(iterators)
  library(doParallel)
  cores <- detectCores(logical=F)-1
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  cvvs1 = c()   
  n_cv = 100
  # s_matrix_HC = subset(s_lable01,g_label =='HC')  # 619*845
  # s_matrix_CRC = subset(s_lable01,g_label =='CRC') # 605*845
  s_matrix_HC2 = s_matrix_HC[,-c(1,2,3)]
  s_matrix_CRC2 = s_matrix_CRC[,-c(1,2,3)]
  
  t1=Sys.time()
  system.time(
    res_cvvs <- foreach(k = 1:n_cv, #.combine='rbind',
                        .packages = c('ROCR','glmnet','kernlab','randomForest','xgboost')) %dopar%
      {
        # k = 1; i = 1
        # n_CRC <- dim(s_matrix_CRC2)[1]
        s_matrix_HC3 <- s_matrix_HC2[sample(nrow(s_matrix_HC2)),] # [1:n_CRC]
        s_matrix_CRC3 <- s_matrix_CRC2[sample(nrow(s_matrix_CRC2)),]
        output_cvvs = CV_variable_selection(s_matrix_HC3,s_matrix_CRC3,5)
        # s_matrix_HC1 <- output_cvvs[[1]]
        # s_matrix_CRC1 <- output_cvvs[[2]]
        cvvs = output_cvvs[[3]]
      }
  )
  t2=Sys.time() 
  print(t2-t1)
  stopImplicitCluster()
  stopCluster(cl)
  
  for(i in 1:n_cv){
    cvvs1 = c(cvvs1,res_cvvs[[i]])
  }
  # names(table(cvvs1)[table(cvvs1)>=100*0.3 & table(cvvs1)<100*0.4])
  cvvs1_1 = names(table(cvvs1)[table(cvvs1)>=n_cv*0.05])   # 0.5 162??0.95 138??0.6 159
  cvvs1_2 = as.numeric(cvvs1_1)
  length(cvvs1_2)
  cvvs1_2
  Genus_names[cvvs1_2]
  # group_gm3 = group_gm1[cvvs1_2,]
  # group_gm = group1[cvvs,]
  # write.csv(group_gm3,file = "group_gm3.csv",row.names = F)
  
  # 看看变量选择后的属种由多少时单独由种类选择出来的？
  group_gm_specials <- read.csv(file = "group_gm3.csv")
  Genus_names[cvvs1_2]
  
  number_genes_specials <- group_gm_specials[group_gm_specials$Genus %in% Genus_names[cvvs1_2],]
  group_gm_specials$Species.1
  
  s_matrix_CRC1 <- s_matrix_CRC[,c(1,2,3,cvvs1_2+3)]
  s_matrix_HC1 <-s_matrix_HC[,c(1,2,3,cvvs1_2+3)]
  # write.csv(s_matrix_CRC1,file = "s_matrix_CRC1.csv",row.names = F)
  # write.csv(s_matrix_HC1,file = "s_matrix_HC1.csv",row.names = F)
  
  
  
  # table(cvvs1)[table(cvvs1)<60]  
}
# group_gm3 = read.csv("C:/Users/duorou/Desktop/Microbial_data/group_gm3.csv")
# cvvs1_2 = group_gm3$number.y
# s_matrix_CRC1 = read.csv("C:/Users/duorou/Desktop/Microbial_data/s_matrix_CRC1.csv")
# s_matrix_HC1 = read.csv('C:/Users/duorou/Desktop/Microbial_data/s_matrix_HC1.csv')

######### heterogeneity test ####

if(FALSE){
  
  group_gm1 = read.csv("C:/Users/duorou/Desktop/Microbial_data/group_gm1.csv")  
  group_gm3 = group_gm1[cvvs1_2,]
  
  # Group regression was used to test for the presence of heterogeneity.
  s_matrix_jy = s_lable01[,group_gm3$number.y + 3]
  dim(s_matrix_jy)  #  1224 159 [1] 1064  139
  t1=Sys.time()
  output_vartest2 = vartest_CRC2(s_matrix_jy,s_lable_ID$CRCorHC,5)
  t2=Sys.time()
  print(t2-t1)
  matrix_coef2 = output_vartest2[[1]]
  n_coef2 = output_vartest2[[2]]
  n_coef_mean2 = output_vartest2[[4]]
  n_byz = output_vartest2[[5]]
  table(n_byz)
  n_byz_order = names(table(n_byz)[table(n_byz)>=6])
  n_byz_order = as.numeric(n_byz_order)
  length(n_byz_order)   # >=5 63 >=6 53; 1064 >=6 54
  group_gm3$number.y[n_byz_order]   #  917 
  
  s_matrix_z = s_matrix_jy
  
  s_matrix_jy1 = as.data.frame(s_matrix_jy)  
  s_matrix_jy1 <- as.data.frame(apply(s_matrix_jy1,2,gyh))
  #train_data = as.matrix(train_data1)
  s_matrix_z = as.data.frame(s_matrix_jy1)
  
  s_matrix_HC_z = subset(s_matrix_z,s_lable_ID$CRCorHC =='HC')  # 524*200
  s_matrix_CRC_z = subset(s_matrix_z,s_lable_ID$CRCorHC =='CRC') # 503*200
  dim(s_matrix_HC_z) 
  
  # The t-test was used to test for the existence of heterogeneity.
  m = length(cvvs1_2)
  p = rep(1,m)
  for(i in 1:m){
    s_CRC = s_matrix_CRC_z[,i]
    s_HC = s_matrix_HC_z[,i]
    x <- c(s_CRC,s_HC)
    group <- c(rep("CRC",524),rep("HC",503))
    test_value = t.test(s_CRC,s_HC,paired = FALSE,var.equal = F)
    p[i] = test_value$p.value
  } 
  length(p[p>0.05])  # 83/159; 1064 72/139
  group_gm3$number.y[p>0.05]
  
  # m = 159
  # p1 = rep(1,m)
  # for(i in 1:m){
  #   s_CRC = s_matrix_CRC_z[,i]
  #   s_HC = s_matrix_HC_z[,i]
  #   s_HC = mean(s_matrix_HC_z[,i])
  #   test_value1 = t.test(s_CRC, mu=s_HC,paired = FALSE,var.equal = F)
  #   p1[i] = test_value1$p.value
  # } 
  # length(p1[p1>0.05])  # 83/159
  # group_gm3$number.y[p1>0.05]
  
  order_test1 = group_gm3$number.y[n_byz_order] 
  order_test2 = group_gm3$number.y[p>0.05]
  length(order_test1[order_test1 %in% order_test2]) # ?ص??ı?��
  order_test_ch = order_test1[order_test1 %in% order_test2]
  group_gm4 = subset(group_gm3,group_gm3$number.y %in% order_test_ch) 
  write.csv(group_gm4,'group_gm4.csv',row.names = F)
  group_gm5 = group_gm3[p>0.05,]   
  group_gm4$Class
  table(group_gm5$Order)
  
  
} 

#####  Classification model  ####################  

if(TRUE){
  
  library(foreach)
  library(parallel)
  library(iterators)
  library(doParallel)
  cores <- detectCores(logical=F)-1
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)
  i = 1
  leni = 20
  r_len = 5     # CV = 5
  CV_ncol = 4*3  
  n_cor = 5   
  CV.xgboost <-CV.rf <-CV.logistic <-CV.svm  <- matrix(0,nrow=leni,ncol=CV_ncol)
  CV.xgboost1<-CV.rf1<-CV.logistic1<-CV.svm1 <- matrix(0,nrow=leni,ncol=CV_ncol)
  Correct.xgboost <-Correct.rf <-Correct.logistic <-Correct.svm  <- matrix(0,nrow=r_len,ncol=CV_ncol)
  Correct.xgboost1<-Correct.rf1<-Correct.logistic1<-Correct.svm1 <- matrix(0,nrow=r_len,ncol=CV_ncol)
  cvvs1 = c()
  # s_matrix_CRC1 <- s_matrix_CRC[,cvvs1_2]
  # s_matrix_HC1 <-s_matrix_HC[,cvvs1_2]
  t1=Sys.time()
  system.time(
    #.combine 
    res <- foreach(k = 1:leni, #.combine='rbind',
                   .packages = c('ROCR','glmnet','kernlab','randomForest','xgboost','e1071')) %dopar%
      { 
        s_matrix_HC1 <- s_matrix_HC1[sample(nrow(s_matrix_HC1)),] # [1:n_CRC]
        s_matrix_CRC1 <- s_matrix_CRC1[sample(nrow(s_matrix_CRC1)),]
        for(i in 1:r_len){
          output_divide = divide_train_test(s_matrix_HC1,s_matrix_CRC1,r_len,i)
          train_data = output_divide[[1]]   # 801*100 ??��Ϊ??,ԭʼ????
          test_data = output_divide[[2]]    # 
          label_train = output_divide[[3]]  #    801
          label_test = output_divide[[4]]
          
          output_bi = binary_matrix_conversion32(train_data,test_data,label_train,label_test,n_cor)
          
          train_data12 = as.data.frame(train_data)  # ?Ա?��??һ??
          train_data <- as.data.frame(apply(train_data12,2,gyh))
          #train_data = as.matrix(train_data1)
          test_data12 = as.data.frame(test_data)
          test_data <- as.data.frame(apply(test_data12,2,gyh))
          
          train_data_01 = output_bi[[1]]  # 801 10  ?任?????? 823*197
          test_data_01 = output_bi[[2]]   # 199 10             204*167
          train_data_w = train_data_01*train_data[,-c(1,2,3)]
          test_data_w = test_data_01*test_data[,-c(1,2,3)]
          train_data_w = cbind(train_data[,c(1,2,3)],train_data_w)
          test_data_w = cbind(test_data[,c(1,2,3)],test_data_w)
          
          output_svm = classification_svm(train_data_w,test_data_w,label_train,label_test,0.001,2)
          output_logistic = classification_logistic(train_data_w,test_data_w,label_train,label_test)
          output_rf = classification_rf(train_data_w,test_data_w,label_train,label_test,13) #17
          output_xgboost = classification_xgboost(train_data_w,test_data_w,label_train,label_test,5,0.02) # 2,0.05
          for(j in 1:4){
            Correct.svm[i,j+4] = output_svm[[j]]
            Correct.logistic[i,j+4] = output_logistic[[j]]
            Correct.rf[i,j+4] = output_rf[[j]]
            Correct.xgboost[i,j+4] = output_xgboost[[j]]
          }
          output_svm1 = classification_svm(train_data,test_data,label_train,label_test,0.01,1)
          output_logistic1 = classification_logistic(train_data,test_data,label_train,label_test)
          output_rf1 = classification_rf(train_data,test_data,label_train,label_test,12)
          output_xgboost1 = classification_xgboost(train_data,test_data,label_train,label_test,5,0.02)  # 6,0.05) # 5,0.1)
          for(j in 1:4){
            Correct.svm[i,j] = output_svm1[[j]]
            Correct.logistic[i,j] = output_logistic1[[j]]
            Correct.rf[i,j] = output_rf1[[j]]
            Correct.xgboost[i,j] = output_xgboost1[[j]]
          }
        }
        CV.svm[k,] = apply(Correct.svm,2,mean,na.rm=TRUE)      # testing set
        CV.logistic[k,] = apply(Correct.logistic,2,mean,na.rm=TRUE)
        CV.rf[k,] = apply(Correct.rf,2,mean,na.rm=TRUE)
        CV.xgboost[k,] = apply(Correct.xgboost,2,mean,na.rm=TRUE)
        # print(k)
        CV1 = list(CV.svm1[k,],CV.logistic1[k,],CV.rf1[k,],CV.xgboost1[k,])
        CV  = list(CV.svm[k,],CV.logistic[k,],CV.rf[k,],CV.xgboost[k,])
        # out_foreach = list(cvvs,CV1,CV)
        out_foreach = list(CV1,CV)
        result <- out_foreach 
      }
  )
  t2=Sys.time() 
  print(t2-t1)
  stopImplicitCluster()
  stopCluster(cl)
  
  # Data sorting
  for(k in 1:leni){
    # k = 1,j = 1
    CV.svm1[k,]      = res[[k]][[1]][[1]]   # training set 
    CV.logistic1[k,] = res[[k]][[1]][[2]]
    CV.rf1[k,]       = res[[k]][[1]][[3]]
    CV.xgboost1[k,]  = res[[k]][[1]][[4]]
    CV.svm[k,]       = res[[k]][[2]][[1]]      # testing set
    CV.logistic[k,]  = res[[k]][[2]][[2]]
    CV.rf[k,]        = res[[k]][[2]][[3]]
    CV.xgboost[k,]   = res[[k]][[2]][[4]]
  }
  apply(CV.svm1,2,mean,na.rm=TRUE)
  apply(CV.logistic1,2,mean,na.rm=TRUE)
  apply(CV.rf1,2,mean,na.rm=TRUE)
  apply(CV.xgboost1,2,mean,na.rm=TRUE)
  apply(CV.svm,2,mean,na.rm=TRUE)
  apply(CV.logistic,2,mean,na.rm=TRUE)
  apply(CV.rf,2,mean,na.rm=TRUE)
  apply(CV.xgboost,2,mean,na.rm=TRUE)
  
  biaoge = matrix(0,nrow=8,ncol=12)
  biaoge[1,] = apply(CV.svm1,2,mean,na.rm=TRUE)
  biaoge[2,] = apply(CV.logistic1,2,mean,na.rm=TRUE)
  biaoge[3,] = apply(CV.rf1,2,mean,na.rm=TRUE)
  biaoge[4,] = apply(CV.xgboost1,2,mean,na.rm=TRUE)
  biaoge[5,] = apply(CV.svm,2,mean,na.rm=TRUE)
  biaoge[6,] = apply(CV.logistic,2,mean,na.rm=TRUE)
  biaoge[7,] = apply(CV.rf,2,mean,na.rm=TRUE)
  biaoge[8,] = apply(CV.xgboost,2,mean,na.rm=TRUE)
  # setwd("C:/Users/Administrator/Desktop")
  biaoge
  # write.csv(biaoge,file = "biaoge_moni.csv",row.names = F)
  
}

#####  Association Rules Mining  ####################

if(TRUE){ 
  n_cor = 5
  r_len = 5
  label_test1<-label_train1<-test_data1<-train_data1 <- list()
  test_data_w1<-train_data_w1<-test_data_011<-train_data_011 <- list()
  for(k in 1:r_len){
    output_divide = divide_train_test(s_matrix_HC1,s_matrix_CRC1,r_len,k)
    train_data = output_divide[[1]]   
    test_data = output_divide[[2]]    
    label_train = output_divide[[3]]  
    label_test = output_divide[[4]]
    # t9=Sys.time()
    output_bi = binary_matrix_conversion32(train_data,test_data,label_train,label_test,n_cor)
    # t10=Sys.time()
    # print(t10-t9)
    
    train_data_11 = as.data.frame(train_data)  # ?Ա?��??һ??
    train_data <- as.data.frame(apply(train_data_11,2,gyh)) 
    test_data_11 = as.data.frame(test_data)
    test_data <- as.data.frame(apply(test_data_11,2,gyh))
    
    train_data_01 = output_bi[[1]]  # 801 10  ?任?????? 823*197
    test_data_01 = output_bi[[2]]   # 199 10             204*167
    train_data_w = train_data_01*train_data[,-c(1,2,3)]
    test_data_w = test_data_01*test_data[,-c(1,2,3)]
    
    train_data1[[k]] = train_data
    test_data1[[k]] = test_data
    label_train1[[k]] = label_train
    label_test1[[k]] = label_test
    train_data_011[[k]] = train_data_01
    test_data_011[[k]] = test_data_01
    train_data_w1[[k]] = train_data_w
    test_data_w1[[k]] = test_data_w
  }
  
  #  test_data_01 
  test_data_01_rbind = rbind(test_data_011[[1]],
                             test_data_011[[2]],
                             test_data_011[[3]],
                             test_data_011[[4]],
                             test_data_011[[5]])
  label_test1_bind = c(label_test1[[1]],label_test1[[2]],
                       label_test1[[3]],label_test1[[4]],
                       label_test1[[5]])
  
  dim(test_data_01_rbind)
  vs_select = cvvs1_2
  rules1_js_5 = c()
  rules2_js_5 = c()
  rules1_js_list = list()
  rules2_js_list = list()
  t3=Sys.time()
  for(i in 1:5){
    n_diso = sample(nrow(test_data_01_rbind))
    test_data_01_rbind1 = test_data_01_rbind[n_diso,]
    label_test1_bind1 = label_test1_bind[n_diso]
    t31=Sys.time()
    output_rules1 = arules_data01_label(test_data_01_rbind,label_test1_bind,
                                        vs_select,0.01,0.8) # output_cvvs[[3]][1:10]
    t41=Sys.time()
    print(t41-t31)
    rules1_js = output_rules1[[1]]  
    rules2_js = output_rules1[[2]]
    rules1_js_list[[i]] = rules1_js
    rules2_js_list[[i]] = rules2_js
    rules1_js_5 = c(rules1_js_5,rules1_js$lhs)
    rules2_js_5 = c(rules2_js_5,rules2_js$lhs)
  }
  t4=Sys.time()
  print(t4-t3)
  
  length(table(rules1_js_5)[table(rules1_js_5)>=3])
  length(table(rules2_js_5)[table(rules2_js_5)>=3])
 
  ch_lhs_CRC= names(table(rules1_js_5)[table(rules1_js_5)>=5])
  a_rules = rules1_js_list[[1]]
  CRC_rules58 = a_rules[a_rules$lhs %in% ch_lhs_CRC,] 
  
  ch_lhs_HC= names(table(rules2_js_5)[table(rules2_js_5)>=5])
  b_rules = rules2_js_list[[1]]
  HC_rules58 = b_rules[b_rules$lhs %in% ch_lhs_HC,] 
  
  # The redundancy problem is reduced by modifying different screening rules
  table_CRC <- c()
  table_HC <- c()
  for (i in c(0.01,.05,.1)) {
    for(j in c(.8,.85,.90,.95)){
      table_CRC <- c(table_CRC,dim(CRC_rules58[CRC_rules58$support>=i & CRC_rules58$confidence >=j,])[1])
      table_HC <- c(table_HC,dim(HC_rules58[HC_rules58$support>=i & HC_rules58$confidence >=j,])[1])
    }
  }
  cbind(table_CRC,table_HC)
}
