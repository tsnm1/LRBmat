# model functions

sample2_xn <- function(n,p,u,rho){
  # Description: simulated sampling function.
  # Input:  n: the number of samples.
  #         p: the number of variables.
  #         u: the mean.
  #         rho: the correlation.
  # Output: x1: the generated data.
  
  # n=100;p=50;rho = 0.9;u = 0;
  sigma1 <- matrix(0,nrow=p,ncol=p)  
  for(i in 1:p){
    for(j in 1:p){
      sigma1[i,j] = rho^(abs(i-j))
    }
  }
  # sigma = cov2cor(sigma1)
  x1 <- rmvnorm(n,mean = rep(u,p),sigma = sigma1)  
  x1 = as.data.frame(x1)
  #cor_x = cor(x,method = "spearman")
  rm(list=c('sigma1'))
  gc()
  return(x1)
}

Calculate_Bij <- function(data,label,n_cor){
  # Description: binary matrix generation algorithm.
  # Input:  data: the initial data matrix, Individual * covariates;.
  #         label: labels of observations.
  #         n_cor: number of most relevant covariates.
  # Output: x_01: the binary matrix corresponding to the initial data.
  
  x = data
  label = ifelse(label<=0,-1,1)
  x_01 = matrix(0,nrow = dim(x)[1],ncol = dim(x)[2])
  cor_x_01 = cor(x,method = "spearman")  # Correlation matrix
  # label_train = sample(c(1,-1),size=n,replace = TRUE,prob = c(500,500))
  # pre1 = label_train
  for(j in 1:dim(x)[2]){
    sy = order(abs(cor_x_01[j,]),decreasing=TRUE)[1:n_cor]   #  TRUE
    train_data11 = x[,sy]
    train_data11 = cbind(train_data11,train_data11[,1]*train_data11[,-1])
    train_data11 = as.matrix(train_data11)
    # colnames(train_data11) <- paste('a', c(1:dim(train_data11)[2]),sep = "", collapse = NULL)
    pfit = cv.glmnet(train_data11,label,family = "binomial",type.measure ="class",nfolds = 5)  # 0.154
    pred1 = predict(pfit,train_data11,s = 'lambda.min',type = "class") 
    x_01[,j] <- pred1
  }
  x_01 = apply(x_01,2,as.numeric)
  rm(list=c('data','x','cor_x_01','train_data11'))
  gc()
  return(x_01)
}

divide_train_test <- function(s_matrix_HC,s_matrix_CRC,r_len,k_i) {
  # Description: data set partition function.
  # Input:  s_matrix_HC: observations classified as HC or 0;
  #         s_matrix_CRC: observations classified as CRC or 1;
  #         r_len: r_len<1 denotes random selection of test set, otherwise denotes cross-validation.
  #         k_i: the ith times of the cross-validation.
  # Output: train_data, test_data: training and test data in part i.
  #         p_train, p_test: corresponding labels.
  
  m = dim(s_matrix_HC)[2]   
  n0 = dim(s_matrix_HC)[1]  
  n1 = dim(s_matrix_CRC)[1]  
  s_matrix11 = rbind(s_matrix_HC,s_matrix_CRC)  
  lable_s =c(rep(0, n0),rep(1, n1))   #label
  # s_matrix11 = as.data.frame(t(s_matrix11)) 
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
  
  n_test_data = sample(nrow(test_data))
  n_train_data = sample(nrow(train_data))
  test_data <- test_data[n_test_data,]; p_test =p_test[n_test_data]
  train_data <- train_data[n_train_data,]; p_train =p_train[n_train_data]
  
  output_divide = list(train_data,test_data,p_train,p_test)
  rm(list=c('m','n0','n1','s_matrix11','lable_s','s_matrix1','c1','c2','test_data','train_data'))
  gc()
  return(output_divide)
}  

CV_variable_selection <- function(s_matrix_HC,s_matrix_CRC,r_len){
  # Description: binary matrix generation algorithm.
  # Input: 
  #   s_matrix_HC, s_matrix_CRC: Input data before variable selection.
  #   r_len: Cross validation parameters.
  # Output: 
  #   s_matrix_HC1, s_matrix_CRC1: Output data after variable selection.
  #   cv.select: Summary of covariates selected by p-value.
  #   table(n_select1): Covariates selected by p-value. 
  #   n_select2: The first n covariates selected.
  
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
    n_select = output_chi[[3]]      # s_p_order
    n_select10 = output_chi[[5]]    # s_p_n  
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
  # Description: binary matrix generation algorithm.
  # Input:  
  #   train_data, test_data: Input data before variable selection.
  #   label_train: Correspanding labels.
  #   pa: p-value.
  #   n: Select the most likely n covariates.
  # Output: 
  #   train_data, test_data: Output data after variable selection.
  #   s_p_order: The ordinal number of the covariate satisfying the p-value condition.
  #   s_len: Number of covariates that satisfy the p-value condition.
  #   s_p_n: The first n covariates that satisfy the condition.
  #   a: The p value of the covariate.
  #   train_data_n, test_data_n: The output data corresponding to the selected N covariates.
  
  s_lable_01 = data.frame(train_data)  #  
  # s_lable_01 = ifelse(s_lable_01>0,1,0)
  # s_lable_01 = ifelse(s_lable_01>=0,1,0)
  s_lable_01[s_lable_01>0] = 1
  s_lable_01[s_lable_01!=1] = 0       
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
  s_p_n = order(a)[1:n]       
  s_len = length(a[a<=pa])   # Number of selected microorganisms
  #Select samples that have passed the microbiological test;
  s_train=apply(as.data.frame(t(train_data)),2,as.numeric)  #  
  s_test=apply(as.data.frame(t(test_data)),2,as.numeric)
  train_data = as.data.frame(t(subset(s_train,a<=pa)))      
  test_data = as.data.frame(t(subset(s_test,a<=pa)))        
  train_data_n = as.data.frame(t(s_train[sort(s_p_n),]))   
  test_data_n = as.data.frame(t(s_test[sort(s_p_n),]))
  # save(s_matrix,file="s_matrix.Rdata")
  # load("s_matrix.Rdata")          
  s_p_order = c(1:size[2])[a<=pa]    
  # 
  output_chi = list(train_data,test_data,s_p_order,s_len,s_p_n,a,train_data_n,test_data_n)
  rm(list=c('size','train_data','test_data','a','train_data_n','test_data_n'))
  gc()
  return(output_chi)
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

binary_matrix_conversion32 <- function(train_data,test_data,label_train,label_test,n_cor){
  # Description: the binary matrix generation algorithm corresponds to the training and test data.
  # Input:  
  #   train_data, test_data: training and test data.
  #   label_train, label_test: corresponding labels.
  #   n_cor: number of most relevant covariates.
  # Output:
  #   train_data,test_data: training and test data.
  #   label_train,label_test: corresponding labels.
  
  # train_data 
  m = dim(test_data)[2]  # 200
  p_pre = matrix(0,dim(test_data)[1],dim(test_data)[2])
  t_pre = matrix(0,dim(train_data)[1],dim(train_data)[2])
  cor_s = cor(train_data,method = "spearman")  # Correlation matrix
  train_data1 = as.data.frame(train_data)
  test_data1 = as.data.frame(test_data) 
  for(j in 1:m){
    sy = order(abs(cor_s[j,]),decreasing=TRUE)[1:n_cor]   #  TRUE
    # sy = order(abs(cor_s[j,]),decreasing=TRUE)[1:5]   #  TRUE
    #sy = order(cor_s[j,],decreasing=TRUE)[1:5]   #  TRUE
    test_data11 = test_data[,sy]             # 11*n/k
    train_data11 = train_data[,sy]
    test_data11 = cbind(test_data11,test_data11[,1]*test_data11[,-1])
    train_data11 = cbind(train_data11,train_data11[,1]*train_data11[,-1])
    # pfit = glmnet(train_data11,label_train, family = "binomial")  #0.027s
    # s_lambda = min(pfit$lambda)
    # pred1 = predict(pfit,train_data11,s = s_lambda , type = "class") 
    test_data11 = as.matrix(test_data11)
    train_data11 = as.matrix(train_data11)
    pfit = cv.glmnet(train_data11,label_train, family = "binomial",type.measure ="class",nfolds = 5)  # 0.154
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

BICAICglm <- function(fit){
  # Description: Calculation of AIC and BIC.
  # Input:  
  #   fit: Model training results.
  # Output:
  #   AIC_, BIC, AICc: the result of model.
  
  tLL <- fit$nulldev -  deviance(fit)
  # tLL <- -deviance(fit) # 2*log-likelihood
  k <- fit$df        #  length(which(fit$beta!=0))
  n <- nobs(fit)
  AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  AIC_ <- -tLL+2*k
  BIC<-log(n)*k - tLL
  
  res=c(AIC_, BIC, AICc)
  names(res)=c("AIC","BIC","AICc")
  return(res)
}  

Generate_interaction_items <- function(x,order,n_items){
  # Description: Generate the interaction term matrix corresponding to the order of the data matrix.
  # Input:
  #   x: the data matrix.
  #   order: the order of generated interaction items.
  #   n_items: the number of interaction items required.
  # Output:
  #   x1: Interaction terms matrix of corresponding order.
  
  # order = 2
  x1 <- x
  m = dim(x)[2]
  a <- combn(1:m, order)
  a <- a[,sample(1:dim(a)[2],n_items)]
  for(i in 1:dim(a)[2]){ #dim(a)[2]
    x_z <- apply(x[,a[,i]],1,cumprod)
    x1 <- cbind(x1,x_z[dim(x_z)[1],])
  }
  x1 <- x1[,-c(1:m)]
  colnames(x1) <- paste('a',order, c(1:dim(x1)[2]),sep = "", collapse = NULL)
  return(x1)
}

model_comparison11 <- function(data,label,x_01,e,n_cor){
  # Description: Model comparison when the real model is the logistic model with interaction terms.
  # Input:
  #   data: Initial observation data.
  #   label: Data labels.
  #   x_01: The corresponding binary matrix after data generation.
  #   e: The proportion of covariates associated with the interaction term.
  #   n_cor: Number of associated covariates.
  # Output:
  #   r: The results summary.
  #   result0: Results of logistic model without interaction terms.
  #   result1: Results of logistic model with 2-order interaction terms.
  #   result1_1: Results of logistic model with 3-order interaction terms.
  #   result2: Results of logistic model based on binary matrix.
  
  data <- as.data.frame(apply(data,2,gyh1))
  # Data to construct
  X1 = data
  n_v1 = dim(data)[2] * e*2   
  for(i in 1:(n_v1-1)){
    for(j in (i+1):n_v1){
      X1 = cbind(X1,data[,i]*data[,j])
    }
  }
  X1_1 = X1
  if(dim(data)[2]>=3 & dim(data)[2]<=30){
    for(i in 1:(n_v1-2)){
      for(j in (i+1):(n_v1-1)){
        for(k in (j+1):n_v1){
          # print(c(i,j,k))
          X1_1 = cbind(X1_1,data[,i]*data[,j]*data[,k])
          # X1_0 = cbind(X1_0,data[,i]*data[,j]*data[,k])
        }
      }
    }
  }else{
    n_v1 = 30
    for(i in 1:(n_v1-2)){
      for(j in (i+1):(n_v1-1)){
        for(k in (j+1):n_v1){
          # print(c(i,j,k))
          X1_1 = cbind(X1_1,data[,i]*data[,j]*data[,k])
          # X1_0 = cbind(X1_0,data[,i]*data[,j]*data[,k])
        }
      }
    }
  }
  
  B2 = Calculate_Bij(data,label,n_cor)
  # data <- as.data.frame(apply(data,2,gyh1))
  X2 = B2 * data 
  
  X2_1 = x_01 * data
  
  
  colnames(data) <- paste('a', c(1:dim(data)[2]),sep = "", collapse = NULL)
  colnames(X1) <- paste('a', c(1:dim(X1)[2]),sep = "", collapse = NULL)
  colnames(X1_1) <- paste('a', c(1:dim(X1_1)[2]),sep = "", collapse = NULL)
  colnames(X2) <- paste('a', c(1:dim(X2)[2]),sep = "", collapse = NULL)
  # colnames(X2_1) <- paste('a', c(1:dim(X2_1)[2]),sep = "", collapse = NULL)
  result0 <- result1 <- result1_1 <- result2 <- result2_1 <- c()
  
  for(i in 1:2){
    X0 <- as.matrix(data)
    pfit1 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)
    pfit11 <- glmnet(X0,pre11, family = "binomial",lambda = pfit1$lambda.min)
    pred1 = predict(pfit11,X0,s = pfit11$lambda , type = "class")
    table1 = table(pred=pred1,truth=pre11)
    ACC =  (table1[1]+table1[4])/sum(table1)
    BIC1 <- BICAICglm(pfit11)
    result0 <- rbind(result0,c(BIC1,ACC))
    
    X0 <- as.matrix(X1)
    pfit1 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)
    pfit11 <- glmnet(X0,pre11, family = "binomial",lambda = pfit1$lambda.min)
    pred1 = predict(pfit11,X0,s = pfit11$lambda , type = "class")
    table1 = table(pred=pred1,truth=pre11)
    ACC =  (table1[1]+table1[4])/sum(table1)
    BIC1 <- BICAICglm(pfit11)
    result1 <- rbind(result1,c(BIC1,ACC))
    
    X0 <- as.matrix(X1_1)
    pfit1 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)
    pfit11 <- glmnet(X0,pre11, family = "binomial",lambda = pfit1$lambda.min)
    pred1 = predict(pfit11,X0,s = pfit11$lambda , type = "class")
    table1 = table(pred=pred1,truth=pre11)
    ACC =  (table1[1]+table1[4])/sum(table1)
    BIC1 <- BICAICglm(pfit11)
    result1_1 <- rbind(result1_1,c(BIC1,ACC))
    
    X0 <- as.matrix(X2)
    pfit2 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)  
    pfit22 <- glmnet(X0,pre11, family = "binomial",lambda = pfit2$lambda.min)
    pred2 = predict(pfit22,X0,s = pfit22$lambda , type = "class")
    table2 = table(pred=pred2,truth=pre11)
    ACC =  (table2[1]+table2[4])/sum(table2)
    BIC2 <- BICAICglm(pfit22)
    result2 <- rbind(result2,c(BIC2,ACC))
    
    X0 <- as.matrix(X2_1)
    pfit2 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)  
    pfit22 <- glmnet(X0,pre11, family = "binomial",lambda = pfit2$lambda.min)
    pred2 = predict(pfit22,X0,s = pfit22$lambda , type = "class")
    table2 = table(pred=pred2,truth=pre11)
    ACC =  (table2[1]+table2[4])/sum(table2)
    BIC2 <- BICAICglm(pfit22)
    result2_1 <- rbind(result2_1,c(BIC2,ACC))
  }
  # r = rbind(apply(result0, 2, mean,na.rm = TRUE),
  #           apply(result1, 2, mean,na.rm = TRUE),
  #           apply(result1_1, 2, mean,na.rm = TRUE),
  #           apply(result2, 2, mean,na.rm = TRUE),
  #           apply(result2_1, 2, mean,na.rm = TRUE))
  r = rbind(apply(result0, 2, mean,na.rm = TRUE),
            apply(result1, 2, mean,na.rm = TRUE),
            apply(result1_1, 2, mean,na.rm = TRUE),
            apply(result2, 2, mean,na.rm = TRUE))
  output <- list(r,result0,result1,result1_1,result2,result2_1)
  rm(list=c('data','X1','X1_1','X2','X2_1','B2','x_01','X0',
            'result0','result1','result1_1','result2','result2_1'))
  gc()
  return(output)
}

model_comparison3 <- function(data,label,e,n_cor){
  # Description: Model comparison when the real model is the logistic model based on Binary matrix.
  # Input:
  #   data: Initial observation data.
  #   label: Data labels.
  #   e: The proportion of covariates associated with the interaction term.
  #   n_cor: Number of associated covariates.
  # Output: 
  #   r: The results summary.
  #   result0: Results of logistic model without interaction terms.
  #   result1: Results of logistic model with 2-order interaction terms.
  #   result1_1: Results of logistic model with 3-order interaction terms.
  #   result2: Results of logistic model based on binary matrix.
  
  # Data to construct
  X1 = data
  n_v1 = dim(data)[2] * e*2    
  for(i in 1:(n_v1-1)){
    for(j in (i+1):n_v1){
      X1 = cbind(X1,data[,i]*data[,j])
    }
  }
  X1_1 = X1
  if(dim(data)[2]>=3 & dim(data)[2]<=30){
    for(i in 1:(n_v1-2)){
      for(j in (i+1):(n_v1-1)){
        for(k in (j+1):n_v1){
          # print(c(i,j,k))
          X1_1 = cbind(X1_1,data[,i]*data[,j]*data[,k])
          # X1_0 = cbind(X1_0,data[,i]*data[,j]*data[,k])
        }
      }
    }
  }else{
    n_v1 = 30
    for(i in 1:(n_v1-2)){
      for(j in (i+1):(n_v1-1)){
        for(k in (j+1):n_v1){
          # print(c(i,j,k))
          X1_1 = cbind(X1_1,data[,i]*data[,j]*data[,k])
          # X1_0 = cbind(X1_0,data[,i]*data[,j]*data[,k])
        }
      }
    }
  }
  
  # X2 <- data
  # for(jj in 1:2){
  #   B2 = Calculate_Bij(X2,label,n_cor)
  #   X2 = B2 * X2
  # }
  B2 = Calculate_Bij(data,label,n_cor)
  data <- as.data.frame(apply(data,2,gyh1))
  X2 = B2 * data
  
  # X2_1 = x_01 * data
  
  colnames(data) <- paste('a', c(1:dim(data)[2]),sep = "", collapse = NULL)
  colnames(X1) <- paste('a', c(1:dim(X1)[2]),sep = "", collapse = NULL)
  colnames(X1_1) <- paste('a', c(1:dim(X1_1)[2]),sep = "", collapse = NULL)
  colnames(X2) <- paste('a', c(1:dim(X2)[2]),sep = "", collapse = NULL)
  # colnames(X2_1) <- paste('a', c(1:dim(X2_1)[2]),sep = "", collapse = NULL)
  result0 <- result1 <- result1_1 <- result2 <- result2_1 <- c()
  
  for(i in 1:2){
    X0 <- as.matrix(data)
    pfit1 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)
    pfit11 <- glmnet(X0,pre11, family = "binomial",lambda = pfit1$lambda.min)
    pred1 = predict(pfit11,X0,s = pfit11$lambda , type = "class")
    table1 = table(pred=pred1,truth=pre11)
    ACC =  (table1[1]+table1[4])/sum(table1)
    BIC1 <- BICAICglm(pfit11)
    result0 <- rbind(result0,c(BIC1,ACC))
    
    X0 <- as.matrix(X1)
    pfit1 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)
    pfit11 <- glmnet(X0,pre11, family = "binomial",lambda = pfit1$lambda.min)
    pred1 = predict(pfit11,X0,s = pfit11$lambda , type = "class")
    table1 = table(pred=pred1,truth=pre11)
    ACC =  (table1[1]+table1[4])/sum(table1)
    BIC1 <- BICAICglm(pfit11)
    result1 <- rbind(result1,c(BIC1,ACC))
    
    X0 <- as.matrix(X1_1)
    pfit1 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)
    pfit11 <- glmnet(X0,pre11, family = "binomial",lambda = pfit1$lambda.min)
    pred1 = predict(pfit11,X0,s = pfit11$lambda , type = "class")
    table1 = table(pred=pred1,truth=pre11)
    ACC =  (table1[1]+table1[4])/sum(table1)
    BIC1 <- BICAICglm(pfit11)
    result1_1 <- rbind(result1_1,c(BIC1,ACC))
    
    X0 <- as.matrix(X2)
    pfit2 = cv.glmnet(X0,pre11, family = "binomial",type.measure ="class",nfolds = 5)  
    pfit22 <- glmnet(X0,pre11, family = "binomial",lambda = pfit2$lambda.min)
    pred2 = predict(pfit22,X0,s = pfit22$lambda , type = "class")
    table2 = table(pred=pred2,truth=pre11)
    ACC =  (table2[1]+table2[4])/sum(table2)
    BIC2 <- BICAICglm(pfit22)
    result2 <- rbind(result2,c(BIC2,ACC))
    
  }
  r = rbind(apply(result0, 2, mean,na.rm = TRUE),
            apply(result1, 2, mean,na.rm = TRUE),
            apply(result1_1, 2, mean,na.rm = TRUE),
            apply(result2, 2, mean,na.rm = TRUE))
  output <- list(r,result0,result1,result1_1,result2)
  rm(list=c('data','X1','X1_1','X2','B2','X0',
            'result0','result1','result1_1','result2'))
  gc()
  return(output)
}

# Classification function
if(TRUE){
  
  classification_svm <- function(train_data,test_data,label_train,label_test,gamma1,cost1){
    # Description: classification algorithm for SVM.
    # Input:
    #   train_data, test_data: training and test data.
    #   label_train, label_test: corresponding labels.
    #   gamma1, cost1: the SVM parameters.
    # Output: 
    #   TPR: True Positive Rate.
    #   FPR: False Positive Rate.
    #   ACC: Accuracy. 
    #   auc: Area Under Curve.
    
    # train_data = train_data_w
    # test_data = train_data_w
    # label_train = label_train
    # label_test = label_train
    # gamma1 = 0.01
    # cost1 = 3
    #rownames(train_data) = paste('', c(1:dim(train_data)[1]),sep = "", collapse = NULL)
    #rownames(test_data) = paste('', c(1:dim(test_data)[1]),sep = "", collapse = NULL)
    d_train <- data.frame(train_data,y=as.factor(label_train))
    d_test <- data.frame(test_data,y=as.factor(label_test))
    #tuned <- tune.svm(y~.,data=d_train,gamma=10^(-4:1),cost=10^(-2:2))
    # $best.parameters
    # gamma cost
    # 22   0.1   10 
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
    # Description: classification algorithm for SVM.
    # Input:
    #   train_data, test_data: training and test data.
    #   label_train, label_test: corresponding labels.
    #   sigma1,C1: the SVM parameters.
    # Output: 
    #   TPR: True Positive Rate.
    #   FPR: False Positive Rate.
    #   ACC: Accuracy. 
    #   auc: Area Under Curve.
    
    
    #colnames(train_data) = paste('', c(1:dim(train_data)[2]),sep = "", collapse = NULL)
    #colnames(test_data) = paste('', c(1:dim(test_data)[2]),sep = "", collapse = NULL)
    d_train <- data.frame(train_data,y=as.factor(label_train))
    d_test <- data.frame(test_data,y=as.factor(label_test))
    train_svm <- ksvm(y~., data = d_train,type = "C-svc", kernel = "rbfdot",
                      kpar=  list(sigma=sigma1), # "automatic",#  # #
                      C=C1,cross=3,prob.model = TRUE) #  ,cross=5
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
  
  classification_logistic <- function(train_data,test_data,label_train,label_test){
    # Description: classification algorithm for logistic model without interaction terms.
    # Input:
    #   train_data, test_data: training and test data.
    #   label_train, label_test: corresponding labels.
    # Output: 
    #   TPR: True Positive Rate.
    #   FPR: False Positive Rate.
    #   ACC: Accuracy. 
    #   auc: Area Under Curve.
    
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
    # Description: classification algorithm for logistic model with interaction terms.
    # Input:
    #   train_data, test_data: training and test data.
    #   label_train, label_test: corresponding labels.
    #   interactions: order of interaction terms.
    # Output: 
    #   TPR: True Positive Rate.
    #   FPR: False Positive Rate.
    #   ACC: Accuracy. 
    #   auc: Area Under Curve.
    
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
  
  classification_rf <- function(train_data,test_data,label_train,label_test,mtry1){
    # Description: classification algorithm for random forest.
    # Input:
    #   train_data, test_data: training and test data.
    #   label_train, label_test: corresponding labels.
    #   mtry1: the parameters.
    # Output: 
    #   TPR: True Positive Rate.
    #   FPR: False Positive Rate.
    #   ACC: Accuracy. 
    #   auc: Area Under Curve.
    
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
    # Description: classification algorithm SVM.
    # Input:
    #   train_data, test_data: training and test data.
    #   label_train, label_test: corresponding labels.
    #   max_depth1,eta1: the parameters.
    # Output: 
    #   TPR: True Positive Rate.
    #   FPR: False Positive Rate.
    #   ACC: Accuracy. 
    #   auc: Area Under Curve.
    
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

model_interactions <- function(data,e,interactions){
  # Description: To construct covariate interaction term data for the input of the logistic model.
  # Input:
  #   data: the initial data matrix.
  #   e: proportion of covariates selected for interaction terms.
  #   interactions: order of interaction terms.
  # Output:
  #   X1_1: generated interaction term data.
  
  X1_1 = data
  n_v1 = dim(data)[2] * e   # It's either 10 or 50
  for(i in 1:(n_v1-1)){
    for(j in (i+1):n_v1){
      X1_1 = cbind(X1_1,data[,i]*data[,j])
    }
  }
  
  if(interactions==3){
    for(i in 1:(n_v1-2)){
      for(j in (i+1):(n_v1-1)){
        for(k in (j+1):n_v1){
          X1_1 = cbind(X1_1,data[,i]*data[,j]*data[,k])
        }
      }
    }
  }
  # output <- list(r,result0,result1,result10,result2)
  return(X1_1)
}

classification_auc <- function(pred_test,label_test){
  # Description: Obtain the AUC result.
  # Input:
  #   pred_test: The prediction results of the classifier.
  #   label_test: Data truth label.
  # Output:
  #   auc: the auc result.
  # library(ROCR)  # Calculating AUC value
  
  pred_test = as.matrix(pred_test) # Predicted results
  # pred_test[pred_test=='CRC']=1
  # pred_test[pred_test=='HC']=0
  pred_test = as.numeric(pred_test)
  label_test = as.matrix(label_test) # Actual label
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

# Other functions
gyh<-function(x){
  # Description: (-1.5, 1.5) Normalize the function
  # Input:
  #   x: The initial data of observations.
  # Output:
  #   x: The normalized data.
  
  if(length(table(x))!=1){  
    mean_x = (max(x)+min(x))/2
    x = (x-mean_x)*3/(max(x)-min(x))   # -1 1 *2
  }
  return(x)
}

gyh1<-function(x){
  # Description: (0, 1) Normalize the function
  # Input:
  #   x: The initial data of observations.
  # Output:
  #   x: The normalized data.
  
  if(length(table(x))!=1){   
    x = (x-min(x))/(max(x)-min(x))  
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
