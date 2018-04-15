setwd("~/Dropbox (UNC Charlotte)/poplar_gwas_data/GS_data-03272018_Wellington/geno_pheno/");
#install.packages("rrBLUP");
#install.packages("BGLR");
library(rrBLUP);
rm(list = ls())
geno <- read.table("Jmax25.geno",header = T);
pheno <- read.table("pheno_Jmax25",header = T);

geno <- as.matrix(geno);
geno_code <- matrix(".",nrow(geno),ncol(geno));
colnames(geno_code) = colnames(geno);
geno_code[,1:6] = geno[,1:6];

###Convert 0/0 -> -1,0/1 -> 0, 1/1 -> 1;
for (i in 1:nrow(geno)){
  geno_code[i,which(geno[i,] == "0/0")] = -1;
  geno_code[i,which(geno[i,] == "0/1")] = 0;
  geno_code[i,which(geno[i,] == "1/1")] = 1;
}

geno_code_genotype_filter <- NULL;
for(i in 1:nrow(geno_code)){
  tmp = unique(geno_code[i,7:ncol(geno_code)]);
  if(length(tmp) != 1){
    geno_code_genotype_filter <- rbind(geno_code_genotype_filter,geno_code[i,]);
  }
}

geno_used <- t(geno_code_genotype_filter[,7:ncol(geno_code_genotype_filter)]);
p = ncol(geno_used);
print(p);
pheno_used <- pheno[,2,drop=FALSE];
rownames(pheno_used) <- pheno[,1];

##Set the Train and Validation population, we use training (70%) and validation (30%) population sets:
###rrBlup:
predict_accuracy_iter <- NULL;
predict_accuracy_iter_bias <- NULL;
predict_accuracy_pval_iter <- NULL;
for(iter in 1:1000){
  print(iter);
  train_percent = 0.7;
  valid_percent = 0.3;
  train = sample(1:nrow(geno_used),(trunc(nrow(geno_used) * train_percent) +1),replace = FALSE);
  validation = setdiff(1:nrow(geno_used),train);
  pheno_train = as.matrix(pheno_used[train,,drop=FALSE]);
  geno_train = as.matrix(geno_used[train,]);
  class(geno_train) <- "numeric";
  
  pheno_valid = as.matrix(pheno_used[validation,,drop=FALSE]);
  geno_valid = as.matrix(geno_used[validation,]);
  class(geno_valid) <- "numeric";
  
  pheno_trait_res = mixed.solve(pheno_train,Z=geno_train,K=NULL,SE=FALSE,return.Hinv = FALSE);
  marker_effect = matrix(pheno_trait_res$u,ncol=1);
  
  predict_pheno_valid = geno_valid %*% marker_effect;
  predict_pheno_valid = predict_pheno_valid + as.numeric(pheno_trait_res$beta);
  bias = predict_pheno_valid - pheno_valid;
  
  predict_accuracy = cor.test(predict_pheno_valid,pheno_valid,use="complete")$estimate;
  predict_accuracy_pval = cor.test(predict_pheno_valid,pheno_valid,use="complete")$p.value;
  
  predict_accuracy_iter <- c(predict_accuracy_iter,predict_accuracy);
  predict_accuracy_iter_bias <- cbind(predict_accuracy_iter_bias,bias);
  predict_accuracy_pval_iter <- c(predict_accuracy_pval_iter,predict_accuracy_pval);
}

bias_final = colMeans(predict_accuracy_iter_bias);
write.table(predict_accuracy_iter,"predict_accuracy_iter_pheno_Jmax25_rrBlup",quote = F,sep="\t",col.names = F,row.names = F);
write.table(predict_accuracy_iter_bias,"predict_accuracy_iter_bias_Jmax25_rrBlup",quote = F,sep="\t",col.names = F,row.names = F);
write.table(bias_final,"bias_final_Jmax25_rrBlup",quote = F,sep="\t",col.names = F,row.names = F);
write.table(predict_accuracy_pval_iter,"predict_accuracy_pval_iter_Jmax25_rrBlup",quote = F,sep="\t",col.names = F,row.names = F);


####
library(BGLR);
rm(list=ls());
setwd("~/Dropbox (UNC Charlotte)/poplar_gwas_data/");
geno <- read.table("Jmax25.geno",header = T);
pheno <- read.table("pheno_Jmax25",header = T);

geno <- as.matrix(geno);
geno_code <- matrix(".",nrow(geno),ncol(geno));
colnames(geno_code) = colnames(geno);
geno_code[,1:6] = geno[,1:6];

###Convert 0/0 -> -1,0/1 -> 0, 1/1 -> 1;
for (i in 1:nrow(geno)){
  geno_code[i,which(geno[i,] == "0/0")] = -1;
  geno_code[i,which(geno[i,] == "0/1")] = 0;
  geno_code[i,which(geno[i,] == "1/1")] = 1;
}

geno_code_genotype_filter <- NULL;
for(i in 1:nrow(geno_code)){
  tmp = unique(geno_code[i,7:ncol(geno_code)]);
  if(length(tmp) != 1){
    geno_code_genotype_filter <- rbind(geno_code_genotype_filter,geno_code[i,]);
  }
}

geno_used <- t(geno_code_genotype_filter[,7:ncol(geno_code_genotype_filter)]);
p = ncol(geno_used);
print(p);
geno_used <- t(geno_code[,7:ncol(geno_code)]);
class(geno_used) <- "numeric";
pheno_used <- pheno[,2,drop=FALSE];
rownames(pheno_used) <- pheno[,1];
n<-nrow(geno_used); 
p<-ncol(geno_used);

##Set the Train and Validation population, we use training (70%) and validation (30%) population sets:
###Bayesian Lasso, Reference Pérez, Paulino, and Gustavo de Los Campos. "Genome-wide regression & prediction with the BGLR statistical package." Genetics (2014): genetics-114.:
predict_accuracy_iter <- NULL;
predict_accuracy_iter_bias <- NULL;
predict_accuracy_pval_iter <- NULL;

for(iter in 1:1000){
  print(iter);
  train_percent = 0.7;
  valid_percent = 0.3;
  train = sample(1:nrow(geno_used),(trunc(nrow(geno_used) * train_percent) +1),replace = FALSE);
  validation = setdiff(1:nrow(geno_used),train);
  pheno_train = as.matrix(pheno_used[train,,drop=FALSE]);
  geno_train = as.matrix(geno_used[train,]);
  class(geno_train) <- "numeric";
  pheno_valid = as.matrix(pheno_used[validation,,drop=FALSE]);
  geno_valid = as.matrix(geno_used[validation,]);
  class(geno_valid) <- "numeric";
  
  yNA <- pheno_used;
  yNA[validation,] <- NA;
  yNA= as.matrix(yNA);
  
  ETA = list(list(X=geno_used,model='BL'));
  fm<-BGLR(y=yNA,ETA=ETA,nIter=1500, burnIn=500,saveAt="BL");
  predict_pheno_valid = fm$yHat[validation];
  bias = predict_pheno_valid - pheno_valid;
  predict_accuracy = cor.test(predict_pheno_valid,pheno_valid,use="complete")$estimate;
  predict_accuracy_pval = cor.test(predict_pheno_valid,pheno_valid,use="complete")$p.value;
  
  predict_accuracy_iter <- c(predict_accuracy_iter,predict_accuracy);
  predict_accuracy_iter_bias <- cbind(predict_accuracy_iter_bias,bias);
  predict_accuracy_pval_iter <- c(predict_accuracy_pval_iter,predict_accuracy_pval);
  rm(fm);
}  
bias_final = colMeans(predict_accuracy_iter_bias);
write.table(predict_accuracy_iter,"predict_accuracy_iter_Jmax25_BL",quote = F,sep="\t",col.names = F,row.names = F);
write.table(predict_accuracy_iter_bias,"predict_accuracy_iter_bias_Jmax25_BL",quote = F,sep="\t",col.names = F,row.names = F);
write.table(bias_final,"bias_final_Jmax25_BL",quote = F,sep="\t",col.names = F,row.names = F);
write.table(predict_accuracy_pval_iter,"predict_accuracy_pval_iter_Jmax25_BL",quote = F,sep="\t",col.names = F,row.names = F);


#### Random forest:
library(randomForest);
setwd("~/Dropbox (UNC Charlotte)/poplar_gwas_data/");
#install.packages("rrBLUP");
#install.packages("BGLR");
rm(list = ls());
geno <- read.table("Jmax25.geno",header = T);
pheno <- read.table("pheno_Jmax25",header = T);

geno <- as.matrix(geno);
geno_code <- matrix(".",nrow(geno),ncol(geno));
colnames(geno_code) = colnames(geno);
geno_code[,1:6] = geno[,1:6];

###Convert 0/0 -> -1,0/1 -> 0, 1/1 -> 1;
for (i in 1:nrow(geno)){
  geno_code[i,which(geno[i,] == "0/0")] = -1;
  geno_code[i,which(geno[i,] == "0/1")] = 0;
  geno_code[i,which(geno[i,] == "1/1")] = 1;
}
geno_code_genotype_filter <- NULL;
for(i in 1:nrow(geno_code)){
  tmp = unique(geno_code[i,7:ncol(geno_code)]);
  if(length(tmp) != 1){
    geno_code_genotype_filter <- rbind(geno_code_genotype_filter,geno_code[i,]);
  }
}

geno_used <- t(geno_code_genotype_filter[,7:ncol(geno_code_genotype_filter)]);
p = ncol(geno_used);
print(p);
geno_used <- t(geno_code[,7:ncol(geno_code)]);
pheno_used <- pheno[,2,drop=FALSE];
rownames(pheno_used) <- pheno[,1];

##Set the Train and Validation population, we use training (70%) and validation (30%) population sets:
###Random forest:
predict_accuracy_iter <- NULL;
predict_accuracy_iter_bias <- NULL;
predict_accuracy_pval_iter <- NULL;

for(iter in 1:1000){
  print(iter);
  train_percent = 0.7;
  valid_percent = 0.3;
  train = sample(1:nrow(geno_used),(trunc(nrow(geno_used) * train_percent) +1),replace = FALSE);
  validation = setdiff(1:nrow(geno_used),train);
  pheno_train = as.matrix(pheno_used[train,,drop=FALSE]);
  geno_train = as.matrix(geno_used[train,]);
  class(geno_train) <- "numeric";
  
  pheno_valid = as.matrix(pheno_used[validation,,drop=FALSE]);
  geno_valid = as.matrix(geno_used[validation,]);
  class(geno_valid) <- "numeric";
  poplar_trait = cbind(geno_used,pheno_used);
  colnames(poplar_trait) = paste("X", colnames(poplar_trait), sep="");
  
  set.seed(131);
  pheno_trait_res = randomForest(XJmax25 ~ .,data = poplar_trait,subset = train,ntree=500);
  # mtry <- tuneRF(poplar_trait[-1],poplar_trait$XJmax25, ntreeTry=500,
  #                stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE);
  # best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
  # print(mtry)
  # print(best.m)
  
  #plot(pheno_trait_res);
  pred<-predict(pheno_trait_res,poplar_trait[-train,]);
  bias = pred - pheno_valid;
  
  predict_accuracy = cor.test(pred,pheno_valid,use="complete")$estimate;
  predict_accuracy_pval = cor.test(pred,pheno_valid,use="complete")$p.value;
  
  predict_accuracy_iter <- c(predict_accuracy_iter,predict_accuracy);
  predict_accuracy_iter_bias <- cbind(predict_accuracy_iter_bias,bias);
  predict_accuracy_pval_iter <- c(predict_accuracy_pval_iter,predict_accuracy_pval);
}

bias_final = colMeans(predict_accuracy_iter_bias);
write.table(predict_accuracy_iter,"predict_accuracy_iter_pheno_Jmax25_rf",quote = F,sep="\t",col.names = F,row.names = F);
write.table(predict_accuracy_iter_bias,"predict_accuracy_iter_bias_Jmax25_rf",quote = F,sep="\t",col.names = F,row.names = F);
write.table(bias_final,"bias_final_Jmax25_rf",quote = F,sep="\t",col.names = F,row.names = F);
write.table(predict_accuracy_pval_iter,"predict_accuracy_pval_iter_Jmax25_rf",quote = F,sep="\t",col.names = F,row.names = F);


###Box plot:
data <- read.table("plotdata_rrBlup",header=T);
pdf(file = "predict_accuracy_traits_rrBlup.pdf");
boxplot(data, ylab = "Prediction accuracy", xlab ="Different traits",col = c("red","sienna","palevioletred1","royalblue2","purple"),
        names = c("RdLight25","Jmax25","WUEref","Resistwp25","DBH"));
dev.off();

####ggplot:Compare different methods:
library(ggplot2);
#data(ToothGrowth);
data <- read.table("plotdata_comparemethods",header=T);

pdf("GS_comparemethods.pdf",width=12,height = 6);
p1 <- ggplot(data, aes(x=trait, y=value, fill=subp)) +
      geom_boxplot(position=position_dodge(1)) + 
      facet_wrap(~ trait, scale="free") + xlab("Traits") + ylab("prediction accuracy") 
p1
dev.off()    

