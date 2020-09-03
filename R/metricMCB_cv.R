#' @title Calculation of model AUC for Methylation Correlation Blocks using cross validation
#'
#' @description To enable quantitative analysis of the methylation patterns
#' within individual Methylation Correlation Blocks across many samples, a single metric to
#' define the methylated pattern of multiple CpG sites within each block.
#' Compound scores which calculated all CpGs within individual Methylation Correlation Blocks by SVM model
#' were used as the compound methylation values of Methylation Correlation Blocks.
#' @usage metricMCB(MCBset,training_set,Surv,testing_set,Surv.new,Method,silent)
#' @export
#' @param MCBset Methylation Correlation Block information returned by the IndentifyMCB function.
#' @param data_set methylation matrix used for training the model in the analysis.
#' @param Surv Survival function contain the survival information for training.
#' @param nfold fold used in the cross validation precedure.
#' @param Method model used to calculate the compound values for multiple Methylation correlation blocks. Options include "svm" "cox" and "lasso". The default option is "svm" for SVR method.
#' @param seed an integer vector, containing the random number generator (RNG) state for random number generation.
#' @param silent Ture indicates that processing information and progress bar will be shown.
#' @author Xin Yu
#' @keywords Methylation Correlation
#' @examples
#' #import datasets
#' data(demo_survival_data)
#' datamatrix<-create_demo()
#' data(demo_MCBinformation)
#' #select MCB with at least 3 CpGs.
#' demo_MCBinformation<-demo_MCBinformation[demo_MCBinformation[,"CpGs_num"]>2,]
#'
#' trainingset<-colnames(datamatrix) %in% sample(colnames(datamatrix),0.6*length(colnames(datamatrix)))
#' testingset<-!trainingset
#' #create the results using Cox regression. 
#' mcb_cox_res<-metricMCB.cv(MCBset = demo_MCBinformation,
#'                data_set = datamatrix,
#'                Surv = demo_survival_data,
#'                Method = "cox")
#'
#' @return Object of class \code{list} with elements (XXX will be replaced with the model name you choose):
#'  \tabular{ll}{
#'    \code{MCB_matrix} \tab Prediction results of model
#'    \code{auc_results} \tab AUC results for each model. \cr
#'  }
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#'
metricMCB.cv<-function(
  MCBset,
  data_set,
  Surv,
  nfold=10,
  Method=c("mean","svm","cox","lasso","ensemble")[1],
  seed=NA,
  silent=FALSE
){
  if (!is.na(seed)&is.integer(seed)){
    set.seed(seed)
  }
  requireNamespace("stats")
  requireNamespace("survival")
  if (!silent) {
    cat("Start anaylsis, this may take a while...\n")
    show_bar=nrow(MCBset)>1
  }else{
    show_bar=FALSE
  }
  if (show_bar) {
    bar<-utils::txtProgressBar(min = 1,max = nrow(MCBset),char = "#",style = 3)
  }
  if (is.null(Surv)) {
    stop(paste("You must have a survival function to train the data."))
  }
  if (is.integer0(grep("MCB_no|CpGs",colnames(MCBset)))){
    stop(paste("Methylation Correlation Block information in your result must have columns of MCB_no and CpGs. Please check your results."))
  }
  na_or_zero_data<-(is.na(Surv[,1])|Surv[,1]==0)
  if (sum(na_or_zero_data)>0){
    data_set<-data_set[,!na_or_zero_data]
    Surv<-Surv[!na_or_zero_data]
    warning("survival data contains NAs or zero survival times, NAs or data with zero survival times are remove automaticly.")
  }
  # constuction of MCB Method matrix for SVM
  MCB_matrix<-matrix(0,nrow = nrow(MCBset),ncol = ncol(data_set))
  colnames(MCB_matrix)<-colnames(data_set)
  rownames(MCB_matrix)<-as.numeric(MCBset[,'MCB_no'])
  if (!Method %in% c("mean","svm","cox","lasso","ensemble")){
    stop(paste("Method:",Method,"is not supported, see hlep files for the details.",collapse = " "))
  }
  if (Method == "mean"){
    return(metricMCB.mean(MCBset,MCB_matrix,Surv,data_set,show_bar=!silent))
  }
  sp<-sample(1:ncol(data_set),replace = F)
  order_sp<-order(sp)
  data_set<-data_set[,sp]
  folds <- cut(seq(1,ncol(data_set)),breaks=nfold,labels=FALSE)
  #if it has a independent test set create the test_set res set
  FunctionResults<-NULL
  best_auc<-0
  best_model<-NULL
  MCB_model_res<-NULL
  for (mcb in seq_len(nrow(MCBset))) {
    #if (nrow(MCBset)>1){}
    if (show_bar&!silent) {
      utils::setTxtProgressBar(bar, mcb)
    }
    write_MCB<-rep(NA,5)
    #save the mcb number
    write_MCB[1]<-as.numeric(MCBset[mcb,'MCB_no'])
    write_MCB[2]<-MCBset[mcb,'CpGs']
    # build temp variable for saving the results.
    # MCB number
    # aquire information for CpG sites in MCB
    CpGs<-strsplit(MCBset[mcb,'CpGs']," ")[[1]]
    #cat(CpGs)
    model<-NULL
    for (i in seq(unique(folds))) {
      rz<- which(folds==i,arr.ind=TRUE)
      data_used_for_training<-data.frame(t(data_set[CpGs,-rz]))
      data_used_for_testing <-data.frame(t(data_set[CpGs,rz]))
      # train a svm model
      times = Surv[-rz]
      if (Method=="svm") {
        model<-tryCatch(survivalsvm::survivalsvm(times ~ .,
                                                 data_used_for_training,
                                                 gamma.mu = 0.1,
                                                 type = "regression"),
                        error = function(e){warning(paste('SVR can not be built, error occurs:', e));return(NULL)})
      }else if(Method=="cox"){
        model<-tryCatch(survival::coxph(times ~ ., 
                                        data_used_for_training),
                        error = function(e){warning(paste('Cox can not be built, error occurs:', e));return(NULL)})
      }else if(Method=="lasso"){
        model<-tryCatch( glmnet::glmnet(as.matrix(data_used_for_training),
                                           times,
                                           #cox model in lasso was used, note that here cox and lasso penalty were used.
                                           family="cox",
                                           alpha=0.5,
                                           # The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
                                           # (1-alpha)/2||beta||_2^2+alpha||beta||_1
                                           # alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
                                           # type.measure = "AUC"
                                           type.measure= "deviance"
                                           # It uses AUC as the criterion for 10-fold cross-validation.
                                           #foldid = 10
        ),error = function(e){warning(paste('lasso model can not be built, error occurs:', e));return(NULL)})
        # correctional_value=1
        # while ( sum(stats::coef(model, s = model$lambda.min-0.001*(correctional_value-1))>0)<1 &
        #         (model$lambda.min-0.001*(correctional_value-1))>0 ) {
        #   correctional_value=correctional_value*1.25
        # }
        # lambda_min_corrected<-model$lambda.min-0.001*(correctional_value-1)
        lambda_min_corrected<-min(model$lambda)
      }else if(Method=="ensemble"){
        models<-list()
        models$cox<-tryCatch(survival::coxph(times ~ ., 
                                        data_used_for_training),
                        error = function(e){warning(paste('Cox can not be built, error occurs:', e));return(NULL)})
        models$svr<-tryCatch(survivalsvm::survivalsvm(times ~ .,
                                                      data_used_for_training,
                                                      gamma.mu = 0.1,
                                                      type = "regression"),
                             error = function(e){warning(paste('SVR can not be built, error occurs:', e));return(NULL)})
        models$lasso<-tryCatch(glmnet::glmnet(as.matrix(data_used_for_training),
                                        times,
                                        family="cox",
                                        alpha=0.5,
                                        type.measure= "deviance"),
                                error = function(e){warning(paste('lasso model can not be built, error occurs:', e));return(NULL)})
        if (!is.null(models$lasso)){
          models$lasso$lambda_min_corrected<-min(models$lasso$lambda)          
        }
        if ((!is.null(models$cox))&
            (!is.null(models$svr))&
            (!is.null(models$lasso))){
        tr_da<-data.frame(cox=stats::predict(models$cox,data_used_for_training),
                          svr=as.numeric(stats::predict(models$svr,data_used_for_training)$predicted),
                          lasso=stats::predict(models$lasso,as.matrix(data_used_for_training),s=models$lasso$lambda_min_corrected))
        models$ensemble<-survival::coxph(times ~.,tr_da)
        MCB_matrix[mcb,rz]<-stats::predict(models$ensemble,
                                            data.frame(cox=stats::predict(models$cox,data_used_for_testing),
                                                        svr=as.numeric(stats::predict(models$svr,data_used_for_testing)$predicted),
                                                        lasso=stats::predict(models$lasso,as.matrix(data_used_for_testing),s=models$lasso$lambda_min_corrected)))
        }else{
          MCB_matrix[mcb,rz]<-NA
        }
        }
      
      if (!is.null(model)) {
        #predictions
        if (Method=="svm") MCB_matrix[mcb,rz]<-stats::predict(model,data_used_for_testing)$predicted
        if (Method=="cox") MCB_matrix[mcb,rz]<-stats::predict(model,data_used_for_testing)
        if (Method=="lasso"){
          #if you use lambda.1se instead, the penalty of lasso would be larger, leading that most of covariates were removed form the final model.
          MCB_matrix[mcb,rz]<-stats::predict(model,as.matrix(data_used_for_testing),s=lambda_min_corrected)
        }
      }else if (Method %in% c("svm","cox","lasso")&is.null(model)){
        MCB_matrix[mcb,rz]<-NA
      }
    }
    MCB_matrix[mcb,]<-MCB_matrix[mcb,order_sp]
    if (sum(is.na(MCB_matrix[mcb,])) == 0){
      AUC_value<-survivalROC::survivalROC(Stime = Surv[,1],
                                          status = Surv[,2],
                                          marker = MCB_matrix[mcb,],
                                          predict.time = 5,
                                          method = "NNE",
                                          span =0.25*length(Surv)^(-0.20))$AUC
      # if (AUC_value<0.5) AUC_value = 1 - AUC_value
      write_MCB[3]<-AUC_value
      cindex<-survival::survConcordance(Surv ~ MCB_matrix[mcb,])
      write_MCB[4]<-cindex$concordance
      write_MCB[5]<-cindex$std.err
    }else{
      write_MCB[3:5]<-NA
    }
    MCB_model_res<-rbind(MCB_model_res,write_MCB)
  }
  colnames(MCB_model_res)<-c("MCB_no","CpGs","auc","C-index","C-index_SE")
  rownames(MCB_matrix)<-MCB_model_res[,'MCB_no']
  FunctionResults$MCB_matrix<-MCB_matrix
  FunctionResults$auc_results<-MCB_model_res
  return(FunctionResults)
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
