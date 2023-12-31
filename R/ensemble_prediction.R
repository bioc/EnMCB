#' @title fitting function using stacking ensemble model for Methylation Correlation Block
#'
#' @description predict is a generic function for predictions from the results of stacking ensemble model fitting functions.
#' The function invokes particular methods which is the ensemble model described in the reference.
#' @param ensemble_model ensemble model which built by ensemble_model() function
#' @param prediction_data A vector, matrix, list, or data frame containing the predictions (input).
#' @param multiple_results Boolean vector, True for including the single model results.
#' @references
#' Xin Yu et al. 2019 Predicting disease progression in lung adenocarcinoma patients based on methylation correlated blocks using ensemble machine learning classifiers (under review)
#' @export
#' @return Object of numeric class \code{double}
#' @examples 
#' library(survival)
#' #import datasets
#' data(demo_survival_data)
#' datamatrix<-create_demo()
#' data(demo_MCBinformation)
#' #select MCB with at least 3 CpGs.
#' demo_MCBinformation<-demo_MCBinformation[demo_MCBinformation[,"CpGs_num"]>2,]
#' trainingset<-colnames(datamatrix) %in% sample(colnames(datamatrix),0.6*length(colnames(datamatrix)))
#' testingset<-!trainingset
#' #select one MCB
#' select_single_one=1
#' em<-ensemble_model(t(demo_MCBinformation[select_single_one,]),
#'     training_set=datamatrix[,trainingset],
#'     Surv_training=demo_survival_data[trainingset])
#'
#' em_prediction_results<-ensemble_prediction(ensemble_model = em,
#' prediction_data = datamatrix[,testingset])
#'
ensemble_prediction <- function(ensemble_model,prediction_data, multiple_results = FALSE) {
  if (multiple_results) {
    return(ensemble_prediction.m(ensemble_model,prediction_data))
  }
  prediction_data<-prediction_data[ensemble_model$cox$cox_model$CpGs,]
  if (nrow(prediction_data)!=length(rownames(prediction_data))) {
    stop("ERROR: The predition data and the model have wrong dimensions!")
  }
  svm<-stats::predict(ensemble_model$svm$svm_model, data.frame(t(prediction_data)))$predicted
  cox<-stats::predict(ensemble_model$cox$cox_model, data.frame(t(prediction_data)))
  enet<-stats::predict(ensemble_model$enet$enet_model,t(prediction_data),s=ensemble_model$enet$`corrected_lambda(min)`)
  mboost<-stats::predict(ensemble_model$mboost$mboost_model, t(prediction_data))[,1]
  data<-rbind(cox,
              svm,
              t(enet),
              mboost
  )
  rownames(data)<-c('cox','svm','enet','mboost')
  data<-t(data)
  data_f<-as.data.frame(data)
  if (class(ensemble_model$stacking)[1] == "cv.glmnet"){
    if (is.null(ensemble_model$stacking$ensemble_type)){
      return(stats::predict(ensemble_model$stacking, as.matrix(data_f)))
    }else{
      f1<-apply(data_f,1,function(x)e1071::kurtosis(x))
      f2<-apply(data_f,1,function(x)(-log(sd(x),2)))
      data_f<-data.frame(f1cox=data_f$cox*f1,f2cox=data_f$cox*f2,
                         f1svm=data_f$svm*f1,f2svm=data_f$svm*f2,
                         f1enet=data_f$enet*f1,f2enet=data_f$enet*f2,
                         f1mboost=data_f$mboost*f1,f2mboost=data_f$mboost*f2)
      return(stats::predict(ensemble_model$stacking, as.matrix(data_f)))
    }
  }
  else if (class(ensemble_model$stacking)[1] == "cph")
    return(stats::predict(ensemble_model$stacking, as.matrix(data_f)))
  else
    stop("ERROR: The ensemble predition model is invaild !")
    
}
