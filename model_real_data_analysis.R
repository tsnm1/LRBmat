####
# real data analysis.


# Some package.
library(gbm)
library('e1071')   
library("kernlab")
library(tree)
library(randomForest)
library(xgboost)
library(Matrix)    
library(ROCR)     
library(EnvStats) 
library(mvtnorm)  
library(GenSA)     
library(glmnet)
library(mniw)   
library(arules)  
library(arulesViz)

############# Data import and preprocessing  ############

if(TRUE){
  setwd("C:/Users/123456/Desktop")
  s = read.csv("s_label.csv")
  g = read.table("group_1224_112.txt",header=T,na.strings = c("NA"))
  s_lable_ID <- merge(s,g,by="ID")  
  g_label = s_lable_ID$CRCorHC  
  s_lable01 = s_lable_ID[,c(-1,-847,-848)]   # 1224  845
  # write.csv(s_lable01,file = "s_lable01.csv",row.names = F)
  dim(s_lable01)   #[1] 1224  846
  s_matrix_HC = subset(s_lable01,g_label =='HC')  # 619*845
  s_matrix_CRC = subset(s_lable01,g_label =='CRC') # 605*845
  # write.csv(s_matrix_CRC,file = "s_matrix_CRC.csv",row.names = F)
  # write.csv(s_matrix_HC,file = "s_matrix_HC.csv",row.names = F)
  
  # s_matrix_HC = read.csv("s_matrix_HC.csv")
  # s_matrix_CRC = read.csv("s_matrix_CRC.csv")
  
  # group2 = read.csv("group2.csv")
  # group_1 = read.csv("group_1.csv")
  # names(group_1)[1] = 'Species'
  # group_gm1 <- merge(group2,group_1,by="Species",all.y=TRUE)   
  # group_gm1 <- group_gm1[order(group_gm1[,10],decreasing=F),]  
  # write.csv(group_gm1,file = "group_gm1.csv",row.names = F)
  group_gm1 = read.csv("group_gm1.csv")   
  
  
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
  cvvs1 = c()   
  n_cv = 100
  t1=Sys.time()
  system.time(
    res_cvvs <- foreach(k = 1:n_cv, #.combine='rbind',
                        .packages = c('ROCR','glmnet','kernlab','randomForest','xgboost')) %dopar%
      {
        # k = 1; i = 1
        s_matrix_HC <- s_matrix_HC[sample(nrow(s_matrix_HC)),]
        s_matrix_CRC <- s_matrix_CRC[sample(nrow(s_matrix_CRC)),]
        output_cvvs = CV_variable_selection(s_matrix_HC,s_matrix_CRC,5)
        s_matrix_HC1 <- output_cvvs[[1]]
        s_matrix_CRC1 <- output_cvvs[[2]]
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
  cvvs1_1 = names(table(cvvs1)[table(cvvs1)>=100*0.4])   # 0.5 162??0.95 138??0.6 159
  cvvs1_2 = as.numeric(cvvs1_1)
  length(cvvs1_2)
  cvvs1_2
  group_gm3 = group_gm1[cvvs1_2,]
  # group_gm = group1[cvvs,]
  write.csv(group_gm3,file = "group_gm3.csv",row.names = F)
  
  # s_matrix_CRC1 <- s_matrix_CRC[,cvvs1_2]
  # s_matrix_HC1 <-s_matrix_HC[,cvvs1_2]
  # write.csv(s_matrix_CRC1,file = "s_matrix_CRC1.csv",row.names = F)
  # write.csv(s_matrix_HC1,file = "s_matrix_HC1.csv",row.names = F)
  
  # table(cvvs1)[table(cvvs1)<60]  
}
# group_gm3 = read.csv("group_gm3.csv") 
# cvvs1_2 = group_gm1$number
# s_matrix_CRC1 = read.csv("s_matrix_CRC1.csv")
# s_matrix_HC1 = read.csv('s_matrix_HC1.csv')

######### heterogeneity test ####

if(TRUE){
  # group_gm3 = read.csv('group_gm3.csv')  
  # group1 = read.csv('group2.csv')   
  # s_lable01 = read.csv("s_lable01.csv")   
  # s_matrix_CRC1 = read.csv("s_matrix_CRC1.csv")
  # s_matrix_HC1 = read.csv('s_matrix_HC1.csv')
  
  # Group regression was used to test for the presence of heterogeneity.
  s_matrix_jy = s_lable01[,group_gm3$number.y]
  dim(s_matrix_jy)  #  1224 159
  t1=Sys.time()
  output_vartest2 = vartest_CRC2(s_matrix_jy,group1$CRCorHC,5)
  t2=Sys.time()
  print(t2-t1)
  matrix_coef2 = output_vartest2[[1]]
  n_coef2 = output_vartest2[[2]]
  n_coef_mean2 = output_vartest2[[4]]
  n_byz = output_vartest2[[5]]
  table(n_byz)
  n_byz_order = names(table(n_byz)[table(n_byz)>=6])
  n_byz_order = as.numeric(n_byz_order)
  length(n_byz_order)   # >=5 63 >=6 53
  group_gm3$number.y[n_byz_order]   #  917 
  
  s_matrix_z = s_matrix_jy
  
  s_matrix_jy1 = as.data.frame(s_matrix_jy)  
  s_matrix_jy1 <- as.data.frame(apply(s_matrix_jy1,2,gyh))
  #train_data = as.matrix(train_data1)
  s_matrix_z = as.data.frame(s_matrix_jy1)
  
  s_matrix_HC_z = subset(s_matrix_z,group1$CRCorHC =='HC')  # 524*200
  s_matrix_CRC_z = subset(s_matrix_z,group1$CRCorHC =='CRC') # 503*200
  dim(s_matrix_HC_z) 
  
  # The t-test was used to test for the existence of heterogeneity.
  m = 159
  p = rep(1,m)
  for(i in 1:m){
    s_CRC = s_matrix_CRC_z[,i]
    s_HC = s_matrix_HC_z[,i]
    x <- c(s_CRC,s_HC)
    group <- c(rep("CRC",524),rep("HC",503))
    test_value = t.test(s_CRC,s_HC,paired = FALSE,var.equal = F)
    p[i] = test_value$p.value
  } 
  length(p[p>0.05])  # 83/159
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
                   .packages = c('ROCR','glmnet','kernlab','randomForest','xgboost')) %dopar%
      {
        # k = 1; i = 1
        # s_matrix_HC <- s_matrix_HC[sample(nrow(s_matrix_HC)),]
        # s_matrix_CRC <- s_matrix_CRC[sample(nrow(s_matrix_CRC)),]
        # output_cvvs = CV_variable_selection(s_matrix_HC,s_matrix_CRC,5)
        # s_matrix_HC1 <- output_cvvs[[1]]
        # s_matrix_CRC1 <- output_cvvs[[2]]
        # cvvs = output_cvvs[[3]]
        # cvvs1 = c(cvvs1,cvvs)
        s_matrix_HC1 <- s_matrix_HC1[sample(nrow(s_matrix_HC1)),]
        s_matrix_CRC1 <- s_matrix_CRC1[sample(nrow(s_matrix_CRC1)),]
        for(i in 1:r_len){
          output_divide = divide_train_test(s_matrix_HC1,s_matrix_CRC1,r_len,i)
          train_data = output_divide[[1]]   # 801*100 ??��Ϊ??,ԭʼ????
          test_data = output_divide[[2]]    # 
          label_train = output_divide[[3]]  #    801
          label_test = output_divide[[4]]
          
          output_bi = binary_matrix_conversion32(train_data,test_data,label_train,label_test,n_cor)
          
          train_data1 = as.data.frame(train_data)  # ?Ա?��??һ??
          train_data <- as.data.frame(apply(train_data1,2,gyh))
          #train_data = as.matrix(train_data1)
          test_data1 = as.data.frame(test_data)
          test_data <- as.data.frame(apply(test_data1,2,gyh))
          
          train_data_01 = output_bi[[1]]  # 801 10  ?任?????? 823*197
          test_data_01 = output_bi[[2]]   # 199 10             204*167
          train_data_w = train_data_01*train_data
          test_data_w = test_data_01*test_data
          # t5=Sys.time()
          output_svm = classification_svm1(train_data_w,train_data_w,label_train,label_train,0.0001,2)
          output_logistic = classification_logistic(train_data_w,train_data_w,label_train,label_train)
          output_rf = classification_rf(train_data_w,train_data_w,label_train,label_train,13)
          output_xgboost = classification_xgboost(train_data_w,train_data_w,label_train,label_train,2,0.05)
          for(j in 1:4){
            Correct.svm1[i,j+8] = output_svm[[j]]
            Correct.logistic1[i,j+8] = output_logistic[[j]]
            Correct.rf1[i,j+8] = output_rf[[j]]
            Correct.xgboost1[i,j+8] = output_xgboost[[j]]
          }
          output_svm = classification_svm1(train_data_01,train_data_01,label_train,label_train,0.0001,2)
          output_logistic = classification_logistic(train_data_01,train_data_01,label_train,label_train)
          output_rf = classification_rf(train_data_01,train_data_01,label_train,label_train,10)
          output_xgboost = classification_xgboost(train_data_01,train_data_01,label_train,label_train,3,0.05)
          for(j in 1:4){
            Correct.svm1[i,j] = output_svm[[j]]
            Correct.logistic1[i,j] = output_logistic[[j]]
            Correct.rf1[i,j] = output_rf[[j]]
            Correct.xgboost1[i,j] = output_xgboost[[j]]
          }
          output_svm1 = classification_svm1(train_data,train_data,label_train,label_train,0.0001,1)
          output_logistic1 = classification_logistic(train_data,train_data,label_train,label_train)
          output_rf1 = classification_rf(train_data,train_data,label_train,label_train,11)
          output_xgboost1 = classification_xgboost(train_data,train_data,label_train,label_train,5,0.02)
          for(j in 1:4){
            Correct.svm1[i,j+4] = output_svm1[[j]]
            Correct.logistic1[i,j+4] = output_logistic1[[j]]
            Correct.rf1[i,j+4] = output_rf1[[j]]
            Correct.xgboost1[i,j+4] = output_xgboost1[[j]]
          }
          
          output_svm = classification_svm1(train_data_w,test_data_w,label_train,label_test,0.0001,2)
          output_logistic = classification_logistic(train_data_w,test_data_w,label_train,label_test)
          output_rf = classification_rf(train_data_w,test_data_w,label_train,label_test,13) #17
          output_xgboost = classification_xgboost(train_data_w,test_data_w,label_train,label_test,2,0.05) 
          for(j in 1:4){
            Correct.svm[i,j+8] = output_svm[[j]]
            Correct.logistic[i,j+8] = output_logistic[[j]]
            Correct.rf[i,j+8] = output_rf[[j]]
            Correct.xgboost[i,j+8] = output_xgboost[[j]]
          }
          # ??Ԫ???? 1-4
          output_svm = classification_svm1(train_data_01,test_data_01,label_train,label_test,0.0001,2)
          output_logistic = classification_logistic(train_data_01,test_data_01,label_train,label_test)
          output_rf = classification_rf(train_data_01,test_data_01,label_train,label_test,10)
          output_xgboost = classification_xgboost(train_data_01,test_data_01,label_train,label_test,3,0.05)
          for(j in 1:4){
            Correct.svm[i,j] = output_svm[[j]]
            Correct.logistic[i,j] = output_logistic[[j]]
            Correct.rf[i,j] = output_rf[[j]]
            Correct.xgboost[i,j] = output_xgboost[[j]]
          }
          # ԭ???? 5-8
          output_svm1 = classification_svm1(train_data,test_data,label_train,label_test,0.0001,1)
          output_logistic1 = classification_logistic(train_data,test_data,label_train,label_test)
          output_rf1 = classification_rf(train_data,test_data,label_train,label_test,11)
          output_xgboost1 = classification_xgboost(train_data,test_data,label_train,label_test,5,0.02)# 6,0.05) # 5,0.1) 
          for(j in 1:4){
            Correct.svm[i,j+4] = output_svm1[[j]]
            Correct.logistic[i,j+4] = output_logistic1[[j]]
            Correct.rf[i,j+4] = output_rf1[[j]]
            Correct.xgboost[i,j+4] = output_xgboost1[[j]]
          }
          
          # t6=Sys.time()
          # print(t6-t5)
          #cat('t6-t5 = ',t6-t5)
        }
        CV.svm1[k,] = apply(Correct.svm1,2,mean,na.rm=TRUE)   # training set 
        CV.logistic1[k,] = apply(Correct.logistic1,2,mean,na.rm=TRUE)
        CV.rf1[k,] = apply(Correct.rf1,2,mean,na.rm=TRUE)
        CV.xgboost1[k,] = apply(Correct.xgboost1,2,mean,na.rm=TRUE)
        CV.svm[k,] = apply(Correct.svm,2,mean,na.rm=TRUE)      # testing set
        CV.logistic[k,] = apply(Correct.logistic,2,mean,na.rm=TRUE)
        CV.rf[k,] = apply(Correct.rf,2,mean,na.rm=TRUE)
        CV.xgboost[k,] = apply(Correct.xgboost,2,mean,na.rm=TRUE)
        # print(k)
        CV1 = list(CV.svm1[k,],CV.logistic1[k,],CV.rf1[k,],CV.xgboost1[k,])
        CV  = list(CV.svm[k,],CV.logistic[k,],CV.rf[k,],CV.xgboost[k,])
        # out_foreach = list(cvvs,CV1,CV)
        out_foreach = list(CV1,CV)
        result <- out_foreach   # R???᷵??????һ???ı?��
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
  write.csv(biaoge,file = "biaoge_moni.csv",row.names = F)
  
}

#####  Association Rules Mining  ####################

if(TRUE){
  # time difference of 12.9171 mins
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
    train_data_w = train_data_01*train_data
    test_data_w = test_data_01*test_data
    
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
  
  # test_data_01 = test_data_011[[5]]
  # label_test  = label_test1[[5]]
  # vs_select = cvvs  # output_cvvs[[3]]
  # group_gm3 = read.csv('group_gm3.csv')
  # vs_select = group_gm3$number
  vs_select = cvvs1_2
  t3=Sys.time()
  output_rules1 = arules_data01_label(test_data_01_rbind,label_test1_bind,
                                      vs_select,0.01,0.8)# output_cvvs[[3]][1:10]
  t4=Sys.time()
  print(t4-t3)
  
  rules1_js = output_rules1[[1]]  
  rules2_js = output_rules1[[2]]
  write.csv(rules1_js,file = "zhenshi1_rules1.csv",row.names = F)
  write.csv(rules2_js,file = "zhenshi1_rules2.csv",row.names = F)
  
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
  # table(rules1_js_5)[table(rules1_js_5)==5]
  
  #  .05 0.8 
  # > length(table(rules1_js_5)[table(rules1_js_5)>=3])
  # [1] 64786
  # > length(table(rules2_js_5)[table(rules2_js_5)>=3])
  # [1] 2
  ch_lhs_CRC= names(table(rules1_js_5)[table(rules1_js_5)>=5])
  a_rules = rules1_js_list[[1]]
  CRC_rules58 = a_rules[a_rules$lhs %in% ch_lhs_CRC,]
  write.csv(CRC_rules58,file = "CRC_rules18.csv",row.names = F)
  
  ch_lhs_HC= names(table(rules2_js_5)[table(rules2_js_5)>=5])
  b_rules = rules2_js_list[[1]]
  HC_rules58 = b_rules[b_rules$lhs %in% ch_lhs_HC,]
  write.csv(HC_rules58,file = "HC_rules18.csv",row.names = F)
  
  t(group_gm4[,c(9,10)])
  group_gm4$number.x
  group_gm4$number.y
}
