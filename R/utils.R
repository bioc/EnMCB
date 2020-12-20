# env parameter
# platform       x86_64-redhat-linux-gnu
# arch           x86_64
# os             linux-gnu
# system         x86_64, linux-gnu
# status
# major          3
# minor          5.0
# year           2018
# month          04
# day            23
# svn rev        74626
# language       R
# version.string R version 3.5.0 (2018-04-23)
# nickname       Joy in Playing

print_as_data <- function(variables,file) {
  sink(file)
  print(variables)
  sink()
}

auc_cal_cv<- function(single_one,train_set,y_surv_train,nfold=10,seed=NULL) {
  requireNamespace('survival')
  
  #remove all the NA and zero survival data
  na_or_zero_data<-(is.na(y_surv_train[,1])|y_surv_train[,1]==0)
  train_set = train_set[,!na_or_zero_data]
  y_surv_train = y_surv_train[!na_or_zero_data]
  
  cox<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'cox',silent = T,seed = seed)
  svr<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'svm',silent = T,seed = seed)
  enet<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'enet',silent = T,seed = seed)
  em<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'ensemble',silent = T,seed = seed)
  
  auc_COX_cv<-auc_roc(y_surv = y_surv_train,cox$MCB_matrix)
  auc_SVR_cv<-auc_roc(y_surv = y_surv_train,svr$MCB_matrix)
  auc_eNet_cv<-auc_roc(y_surv = y_surv_train,enet$MCB_matrix)
  auc_em_cv<-auc_roc(y_surv = y_surv_train,em$MCB_matrix)
  
  return(c(single_one,
           auc_COX_cv,auc_SVR_cv,auc_eNet_cv,auc_em_cv))
}

auc_cal_cv_tvt <- function(single_one,
                           train_set,validation_set,test_set,
                           y_surv_train,y_surv_validation,y_surv_test,
                           seed=9){
  requireNamespace(survival)
  #remove all the NA and zero survival data
  na_or_zero_data<-(is.na(y_surv_train[,1])|y_surv_train[,1]==0)
  train_set = train_set[,!na_or_zero_data]
  y_surv_train = y_surv_train[!na_or_zero_data]
  
  na_or_zero_data<-(is.na(y_surv_validation[,1])|y_surv_validation[,1]==0)
  validation_set = validation_set[,!na_or_zero_data]
  y_surv_validation = y_surv_validation[!na_or_zero_data]
  
  na_or_zero_data<-(is.na(y_surv_test[,1])|y_surv_test[,1]==0)
  test_set = test_set[,!na_or_zero_data]
  y_surv_test = y_surv_test[!na_or_zero_data]
  
  
  cox<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'cox',silent = T,seed = seed)
  svr<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'svm',silent = T,seed = seed)
  enet<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'enet',silent = T,seed = seed)
  em<-metricMCB.cv(t(single_one),train_set,y_surv_train,Method = 'ensemble',silent = T,seed = seed)
  
  auc_COX_train_cv<-auc_roc(y_surv = y_surv_train,cox$MCB_matrix)
  auc_SVR_train_cv<-auc_roc(y_surv = y_surv_train,svr$MCB_matrix)
  auc_eNet_train_cv<-auc_roc(y_surv = y_surv_train,enet$MCB_matrix)
  auc_em_train_cv<-auc_roc(y_surv = y_surv_train,em$MCB_matrix)
  
  
  vt_data<-cbind(validation_set,test_set)
  vt_survival<-c(y_surv_validation,y_surv_test)
  vt_flag<-c(rep('validation',ncol(validation_set)),
             rep('test',ncol(test_set)))
  
  cox = metricMCB(t(single_one),
                  training_set = train_set,
                  Surv = y_surv_train,
                  testing_set = vt_data,
                  Surv.new = vt_survival,
                  Method = 'cox',silent = T)
  
  svr = metricMCB(t(single_one),
                  training_set = train_set,
                  Surv = y_surv_train,
                  testing_set = vt_data,
                  Surv.new = vt_survival,
                  Method = 'svm',silent = T) 
  
  enet = metricMCB(t(single_one),
                   training_set = train_set,
                   Surv = y_surv_train,
                   testing_set = vt_data,
                   Surv.new = vt_survival,
                   Method = 'enet',silent = T)  
  
  em = EnMCB::ensemble_model(t(single_one),
                             training_set = train_set,
                             Surv_training = y_surv_train,
                             testing_set = vt_data,
                             Surv_testing = vt_survival)
  prid_em<- EnMCB::ensemble_prediction(em,train_set)
  test_and_validation_em<- EnMCB::ensemble_prediction(em,vt_data)
  
  auc_COX_train<-auc_roc(y_surv = y_surv_train,cox$MCB_cox_matrix_training)
  auc_SVR_train<-auc_roc(y_surv = y_surv_train,svr$MCB_svm_matrix_training)
  auc_eNet_train<-auc_roc(y_surv = y_surv_train,enet$MCB_enet_matrix_training)
  auc_em_train<-auc_roc(y_surv = y_surv_train,prid_em)
  
  
  auc_COX_validation<-auc_roc(y_surv = y_surv_validation,cox$MCB_cox_matrix_test_set[vt_flag == 'validation'])
  auc_SVR_validation<-auc_roc(y_surv = y_surv_validation,svr$MCB_svm_matrix_test_set[vt_flag == 'validation'])
  auc_eNet_validation<-auc_roc(y_surv = y_surv_validation,enet$MCB_enet_matrix_test_set[vt_flag == 'validation'])
  auc_em_validation<-auc_roc(y_surv = y_surv_validation,test_and_validation_em[vt_flag == 'validation'])
  
  auc_COX_test<-auc_roc(y_surv = y_surv_test,cox$MCB_cox_matrix_test_set[vt_flag == 'test'])
  auc_SVR_test<-auc_roc(y_surv = y_surv_test,svr$MCB_svm_matrix_test_set[vt_flag == 'test'])
  auc_eNet_test<-auc_roc(y_surv = y_surv_test,enet$MCB_enet_matrix_test_set[vt_flag == 'test'])
  auc_em_test<-auc_roc(y_surv = y_surv_test,test_and_validation_em[vt_flag == 'test'])
  
  return(c(single_one,
           auc_COX_train_cv,auc_SVR_train_cv,auc_eNet_train_cv,auc_em_train_cv,
           auc_COX_train,auc_SVR_train,auc_eNet_train,auc_em_train,
           auc_COX_validation,auc_SVR_validation,auc_eNet_validation,auc_em_validation,
           auc_COX_test,auc_SVR_test,auc_eNet_test,auc_em_test))
}

auc_roc <- function(y_surv,predictions) {
  survivalROC::survivalROC.C(Stime = y_surv[,1],
                             status = y_surv[,2],
                             marker = predictions,
                             predict.time = 5,
                             span =0.25*length(y_surv[,1])^(-0.20))$AUC
}


test_logrank_p_in_KM<-function(mcb_matrix,y_surv){
  p.val_all<-NULL
  for (mcb in rownames(mcb_matrix)){
    Test<-mcb_matrix[mcb,]
    group_sur<-factor(1*(Test>median(Test)),
                      levels=0:1,
                      labels=c("High risk","Low risk"))
    diff_tcga<-survival::survdiff(y_surv ~ group_sur)
    p.val <- 1 - stats::pchisq(diff_tcga$chisq, length(diff_tcga$n) - 1)
    p.val_all<-rbind(p.val_all,p.val)
  }
  names(p.val_all)<-rownames(mcb_matrix)
  p.val_all
}


give_me_the_sen_spc<-function(pred,print=TRUE) {
  tpr<-(pred@tp[[1]]/(pred@tp[[1]]+pred@fn[[1]])) #sensitivity
  fpr<-(pred@fp[[1]]/(pred@tn[[1]]+pred@fp[[1]])) #specificity
  best_one<-which.min(sqrt( (1-tpr)^2+ fpr^2 ))
  if (print == TRUE)  cat(paste(" sensitivity:",tpr[best_one],"\n","specificity:",1-fpr[best_one],"\n"))
  pred@cutoffs[[1]][best_one]#specificity
}

random_matrix <- function(x) {
  x_random=matrix(x[sample(seq_along(x),length(x))],ncol = ncol(x),nrow = nrow(x))
  rownames(x_random)<-rownames(x)
  colnames(x_random)<-colnames(x)
  x_random
}


#'@title create demo matrix
#'@author Xin Yu
#'@description Demo matrix for methylation matrix.
#'@param model Two options, 'all' or 'short' for creating full dataset or very brief demo.
#'@return This function will generate a demo data.
#'@export
#'@examples 
#'demo_set<-create_demo()
create_demo<-function(model=c('all','short')[1]){
  utils::data(demo_data)
  rmax<-matrix(data = rnorm(length(demo_data$rownames)*length(demo_data$colnames), mean=0.5, sd=0.08),
           nrow = length(demo_data$rownames),
           ncol=length(demo_data$colnames),dimnames = list(demo_data$rownames,demo_data$colnames)
           )
  rmax[rownames(demo_data$realdata),]<-demo_data$realdata
  if (model=='short') {
    rmax<-rmax[ceiling(nrow(rmax)*0.8):nrow(rmax),]
  }
  rmax
}



as.numeric_matrix<-function(met_matrixf){
  met_matrixff<-matrix(as.numeric(met_matrixf),
                       nrow = NROW(met_matrixf),
                       ncol = NCOL(met_matrixf))
  rownames(met_matrixff)<-rownames(met_matrixf)
  colnames(met_matrixff)<-colnames(met_matrixf)
  met_matrixff
}

makeEset<-function(met_matrix,cli_data){
  all_col<-NULL
  for (i in seq_len(ncol(met_matrix))) {
    if (length(grep(pattern=tolower(colnames(met_matrix)[i]),cli_data))>0) {
      loca<-grep(pattern=tolower(colnames(met_matrix)[i]),cli_data)[1]
      col<-loca%/%nrow(cli_data)+1
      row<-loca%%nrow(cli_data)
      all_col<-c(all_col,col)
    }
  }
  cli_data_f<-cli_data[,all_col]
  colnames(cli_data_f)<-colnames(met_matrix)
  eSet<-make_eset(met_matrix,data.frame(t(cli_data_f)))
  eSet
}


multiple_time_ROC <- function(Test,y_surv,genesel) {
  requireNamespace("prognosticROC") #version 0.7
  sroclong_all<-NULL
  for (n_time in c(3,5,7)) {
    ROC_res= survivalROC::survivalROC(Stime=y_surv[,1],
                         status=y_surv[,2],
                         marker =as.numeric(Test),
                         lambda = NULL,
                         predict.time = n_time,method = "NNE",span =0.25*length(y_surv[,1])^(-0.20))#
    # AUC at 7 years: 0·670 (0·578–0·762)
    sroclong_all<-ROCdata_save(sroclong_all,ROC_res,mark = paste("AUC at",n_time,"years:",round(ROC_res$AUC,2),collapse = " "))
  }
  gg<-ggplot(sroclong_all, aes(x = FPF, y = TPF, label = c, color = groups))+
    coord_cartesian(ylim=c(0, 1.05),xlim=c(0, 1.05),expand=FALSE)+
    geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_light())+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    theme(legend.position=c(0.8,0.35), legend.justification=c(1,1),legend.title = element_blank())+
    geom_abline(slope=1, colour="black")
  ggsave(filename = paste("Time ROC of ",genesel,".jpeg",sep=""),plot = gg,device ="jpeg" ,
         path = getwd(),dpi = 300,units = "in",width = 5, height = 4.5,
         limitsize=FALSE)
}




draw_pheatmap_matrix<-function(eSet_matrix,group,savename="draw_pheatmap_matrix",scale = FALSE){
  requireNamespace(pheatmap)
  requireNamespace(ggplot2)
  group<-as.character(group)
  df<-as.numeric_matrix(eSet_matrix)
  df<-df[,order(group,decreasing = TRUE)]
  group<-group[order(group,decreasing = TRUE)]
  df[is.na(df)]<-0
  aka2 = data.frame(ID = group)
  rownames(aka2)<-colnames(df)
  aka3 = list(ID = c('TRUE' = "blue", 'FALSE'="red"))

  df<-df[apply(df, 1,stats::sd)>0,]
  if (scale == TRUE) {
    df=scale(df)
  }
   gg<-pheatmap(df,
           annotation_col = aka2,
           annotation_colors = aka3[1],
           annotation_legend = FALSE,
           show_colnames = FALSE, show_rownames = FALSE, cluster_rows = TRUE,
           cluster_cols = FALSE, legend = TRUE,
           clustering_distance_rows = "euclidean", border_color = FALSE)
   ggplot2::ggsave(filename = paste(savename,".jpeg",sep=""),plot = gg,device ="jpeg" ,
          path = getwd(),dpi = 300,units = "in",width = 10, height = 5,
          limitsize=FALSE)
}


give_me_the_lasso_gene<-function(cvfit,s=NULL,plot_figure=TRUE){
  if (plot_figure ==TRUE) {
    plot(cvfit)
    # cat("min: ",cvfit$lambda.min,"\n")
    # cat("1se: ",cvfit$lambda.1se,"\n")
  }
  if (is.null(s)) {
    s=cvfit$lambda.min
  }
  Coefficients <- coef(cvfit, s = s) #0.20
  Active.Index <- which(as.logical(Coefficients != 0))
  Active.Coefficients <- Coefficients[Active.Index]
  #print(Active.Coefficients)
  gene_lasso<-Coefficients@Dimnames[[1]][Active.Index]
  if (plot_figure ==TRUE) {
    plot(cvfit$glmnet.fit,xvar="lambda")
  }
  list(gene_lasso=gene_lasso,Active.Coefficients=Active.Coefficients)
}


plot_ROC<-function(perf){
  sroclong<-data.frame(TPF = perf@y.values[[1]], FPF =perf@x.values[[1]],
                       c = perf@alpha.values[[1]] #,
                       #group = rep("M stage", length(Mayo5[["FP"]]))
  )
  ggplot(sroclong, aes(x = FPF, y = TPF, label = c
                       #, color = group
  )) +
  geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_grey)
}



mutiple_time_ROC <- function(Test,y_surv,genesel) {
  requireNamespace(survivalROC) #version 1.0.3
  requireNamespace(prognosticROC) #version 0.7
  sroclong_all<-NULL
  for (n_time in c(3,5,7)) {
    ROC_res= survivalROC::survivalROC(Stime=y_surv[,1],
                         status=y_surv[,2],
                         marker =as.numeric(Test),
                         lambda = NULL,
                         predict.time = n_time,method = "NNE",span =0.25*length(y_surv[,1])^(-0.20))#
    # AUC at 7 years: 0·670 (0·578–0·762)

    sroclong<-data.frame(TPF = ROC_res[["TP"]], FPF = ROC_res[["FP"]],
                         c = ROC_res[["cut.values"]],
                         groups = rep(paste("AUC at",n_time,"years:",round(ROC_res$AUC,2),collapse = " "),
                                     length(ROC_res[["FP"]])))
    sroclong_all<-rbind(sroclong_all,sroclong)
  }
  gg<-ggplot(sroclong_all, aes(x = FPF, y = TPF, label = c, color = groups))+
    coord_cartesian(ylim=c(0, 1.05),xlim=c(0, 1.05),expand=FALSE)+
    geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_light())+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    theme(legend.position=c(0.8,0.35), legend.justification=c(1,1),legend.title = element_blank())+
    geom_abline(slope=1, colour="black")
  ggsave(filename = paste("Time ROC of ",genesel,".jpeg",sep=""),plot = gg,device ="jpeg" ,
         path = getwd(),dpi = 300,units = "in",width = 5, height = 4.5,
         limitsize=FALSE)
}




#'@title draw survival curve
#'@author Xin Yu
#'@description Draw a survival curve based on survminer package. This is a wrapper function of ggsurvplot.
#'@param exp expression level for gene.
#'@param living_days The survival time (days) for each individual.
#'@param living_events The survival event for each individual, 0 indicates alive and 1 indicates death.
#'Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
#'@param write_name The name for pdf file which contains the result figure.
#'@param title_name The title for the result figure.
#'@param threshold Threshold used to indicate the high risk or low risk.
#'
#'
#'@return This function will generate a pdf file with 300dpi which compare survival curves using the Kaplan-Meier (KM) test.
#'@export
#'@examples 
#' data(demo_survival_data)
#' library(survival)
#' demo_set<-create_demo()
#' draw_survival_curve(demo_set[1,],
#'     living_days = demo_survival_data[,1],
#'     living_events =demo_survival_data[,2],
#'     write_name = "demo_data" )
#'
draw_survival_curve<-function(exp,living_days,living_events,write_name,title_name="",threshold=NA){
  if (is.na(threshold)) threshold=median(exp)
  group_sur<-factor(1*(exp >threshold) ,
                    levels=0:1,
                    labels=c("low risk", "high risk"))
  pdata_gene <- list(time=living_days,
                     event=living_events,
                     gene=exp,
                     group_sur=group_sur)
  fit <- survival::survfit(survival::Surv(time, event) ~ group_sur, data = pdata_gene)
  data.survdiff <- survival::survdiff(survival::Surv(time, event) ~ group_sur, data = pdata_gene)
  p.val = 1 - stats::pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  if (p.val<0.00001){
    p_value="p value < 0.00001"
  }else{
    p_value=paste("p =",round(p.val,5),collapse = " ")
  }
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  #require(ggfortify)
  gg<-survminer::ggsurvplot(fit,data = pdata_gene,  risk.table = TRUE,risk.table.y.text = TRUE,break.time.by=1,
                            fillLineSize=3,survSize=5,surv.linesize=5,pval = FALSE,
                            fontsize=2.5,
                            legend.labs = c("low risk","high risk"),xlab="Time (years)",
                        censor.alpha=0.5,palette=c("#00BFC4","#F8766D"),
                        title=paste(title_name ,"\n",p_value," HR ",round(HR,2),
                                         "(",round(up95,2) ,"-", round(low95,2),  ")",collapse = ""),
  ggtheme=ggplot2::theme(title =ggplot2::element_text(size=8, face='bold'),
          panel.grid.major =ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"),
    legend.position=c(1,1), legend.justification=c(1,1),
    legend.background = ggplot2::element_rect(fill = 'white', colour = 'white'))
  )
  ggplot2::ggsave(filename = paste("survival of ",write_name,".jpeg",sep=""),plot = print(gg),device ="jpeg" ,
         path = getwd(),dpi = 300,units = "in",width = 5, height = 5,
         limitsize=FALSE)
}


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

metricMCB.mean<-function(MCBset,MCB_matrix,Surv,data_set,show_bar=T){
  FunctionResults<-list()
  MCB_model_res<-NULL
  if (show_bar) {
    bar<-utils::txtProgressBar(min = 1,max = nrow(MCBset),char = "#",style = 3)
  }
  for (mcb in seq(nrow(MCBset))) {
    utils::setTxtProgressBar(bar, mcb)
    write_MCB<-rep(NA,5)
    #save the mcb number
    write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
    write_MCB[2]<-MCBset[mcb,'CpGs']
    CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
    MCB_matrix[mcb,]<-colMeans(data_set[CpGs,])
    AUC_value<-survivalROC::survivalROC(Stime = Surv[,1],
                                        status = Surv[,2],
                                        marker = MCB_matrix[mcb,],
                                        predict.time = 5,
                                        method = "NNE",
                                        span =0.25*length(Surv)^(-0.20))$AUC
    write_MCB[3]<-AUC_value
    cindex<-survival::survConcordance(Surv ~ MCB_matrix[mcb,])
    write_MCB[4]<-cindex$concordance
    write_MCB[5]<-cindex$std.err
    MCB_model_res<-rbind(MCB_model_res,write_MCB)
  }
  colnames(MCB_model_res)<-c("MCB_no","CpGs","auc","C-index","C-index_SE")
  rownames(MCB_matrix)<-MCB_model_res[,'MCB_no']
  FunctionResults$MCB_matrix<-MCB_matrix
  FunctionResults$auc_results<-MCB_model_res
  return(FunctionResults)
}

ROCdata_save<-function(origin=NULL,perf,mark="none"){
  sroclong<-data.frame(TPF = perf$TP, FPF =perf$FP, 
                       c = perf$cut.values, 
                       group = rep(mark, length(perf$TP)) )
  sroclong<-rbind(origin,sroclong)
  sroclong
}

ensemble_prediction.m<- function(ensemble_model,predition_data) {
  predition_data<-predition_data[rownames(predition_data) %in% names(coef(ensemble_model$cox$cox_model)),]
  if (nrow(predition_data)!=length(rownames(predition_data))) {
    stop("ERROR: The predition data and the model have wrong dimensions!")
  }
  svm<-stats::predict(ensemble_model$svm$svm_model, data.frame(t(predition_data)))$predicted
  cox<-stats::predict(ensemble_model$cox$cox_model, data.frame(t(predition_data)))
  enet<-stats::predict(ensemble_model$enet$`enet model`,t(predition_data),s=ensemble_model$enet$`corrected lambda(min)`)
  data<-rbind(cox,
              svm,
              t(enet)
  )
  rownames(data)<-c('cox','svm','enet')
  data<-t(data)
  data_f<-as.data.frame(data)
  ensemble=stats::predict(ensemble_model$stacking, data_f)
  return(t(cbind(data,ensemble)))
}