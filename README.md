# EnMCB
Package: EnMCB

Type: Package

Title: Predicting Disease Progression Based on Methylation Correlated Blocks using Ensemble Models
        
Version: 1.11.1

Author: Xin Yu

Maintainer: Xin Yu <whirlsyu@gmail.com> <yuxin@qlnu.edu.cn>

Description: This package is designed to help you to create the methylation correlated blocks using methylation profiles. A stacked ensemble of machine learning models, which combined the Cox regression, support vector regression, Coxboost and elastic-net regression model, can be constructed using this package. You also can choose one of them to build DNA methylation signatures associated with disease progression (survival).

License: GPL-3

Citation:
Xin Yu, De-Xin Kong, EnMCB: an R/bioconductor package for predicting disease progression based on methylation correlated blocks using ensemble models, Bioinformatics, 2021, btab415


Followings are brief instructions for using this package:

You can install and our package via BiocManager as following or downloading source from github.

<pre>
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnMCB")
</pre>


First, you need a methylation data set, currently only most common platform 'Illumina Infinium Human Methylation 450K' is supported. For other platform, for instance, BS sequencing data, users need to 
introduce the annonation file for other own.

You can use your own datasets:

<pre>
eset_met<-some_methylation_datamatrix
</pre>

or use our demo data.

<pre>
data('demo_data',package = "EnMCB")
methylation_dataset<-demo_data$realdata
</pre>


Then, you can automatically run following:

<pre>
library(SummarizedExperiment)

res<-IdentifyMCB(methylation_dataset)

</pre>

You can extract the MCB information,

<pre>

MCB<-res$MCBinformation

</pre>

and select some of MCBs for further modeling.

<pre>

MCB<-MCB[MCB[,"CpGs_num"]>5,]

</pre>

In order to get differentially methylated blocks, one may run following:

<pre>
#simulation for the group data
groups = c(rep("control",200),rep("dis",255))

diffMCB_results <- DiffMCB(methylation_dataset, groups, MCB)

</pre>

In order to build survival models, one may run following:

<pre>
data(demo_survival_data)
survival_data <- demo_survival_data
# sample the dataset into training set and testing set
trainingset<-colnames(methylation_dataset) %in% sample(colnames(methylation_dataset),0.6*length(colnames(methylation_dataset)))

testingset<-!trainingset

#build the models

models<-metricMCB(MCB,
                    training_set=methylation_dataset[,training_set],
                    Surv=survival_data)

#select the best
onemodel<-models$best_cox_model

</pre>                    

Then, you can predict the risk by the model you build:

<pre>
predict(onemodel,methylation_dataset[,testing_set])
</pre>

In order to build ensemble model, one may run following:

<pre>
# You can choose one of MCBs:
select_single_one=1

em<-ensemble_model(MCB[select_single_one,],
                    training_set=methylation_dataset[,training_set],
                    Surv=survival_data)
                    
</pre>
Note that this function only can be used for single MCB only, otherwise the precessing time could be very long.

Then, you can predict the risk by the model you build:

<pre>
ensemble_prediction(ensemble_model = em,
                    predition_data = methylation_dataset[,testing_set])
</pre>

For detailed information, you can find at our references.

