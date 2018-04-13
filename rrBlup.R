setwd("~/Dropbox (UNC Charlotte)/poplar_gwas_data/GS_data-03272018_Wellington/geno_pheno/");
library(rrBLUP);
rm(list = ls());
thresholds = c(1e-05,1e-04,0.001);
traits = c("Jmax25","Rdlight25","Resistwp25","WUEref","Av_Diameter_mm");

for (threshold in 1:length(thresholds)){
  for(trait in 1:length(traits)){
    cat("threshold is",thresholds[threshold],";trait is",traits[trait],"\n");
    geno <- read.table(paste(traits[trait],"_",thresholds[threshold],".geno",sep=""),header = T);
    pheno <- read.table(paste("pheno_",traits[trait],"_",thresholds[threshold],sep=""),header = T);
    
    geno <- as.matrix(geno);
    geno_code <- matrix(".",nrow(geno),ncol(geno));
    colnames(geno_code) = colnames(geno);
    geno_code[,1:6] = geno[,1:6];
    
    ###Convert 0/0 -> -1,0/1 -> 0, 1/1 -> 1;
    for (i in 1:nrow(geno)){
      geno_code[i,which(geno[i,] == "0/0")] = -1;
      geno_code[i,which(geno[i,] == "0/1" | geno[i,] == "1/0")] = 0;
      geno_code[i,which(geno[i,] == "1/1")] = 1;
    }
    
    geno_used <- t(geno_code[,7:ncol(geno_code)]);
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
      cat("threshold is",thresholds[threshold],";trait is",traits[trait],"iteration is",iter,"\n");
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
    write.table(predict_accuracy_iter,paste("predict_accuracy_iter_",traits[trait],"_",thresholds[threshold],"_rrBlup",sep=""),quote = F,sep="\t",col.names = F,row.names = F);
    write.table(predict_accuracy_iter_bias,paste("predict_accuracy_iter_bias_",traits[trait],"_",thresholds[threshold],"_rrBlup",sep=""),quote = F,sep="\t",col.names = F,row.names = F);
    write.table(bias_final,paste("bias_final_",traits[trait],"_",thresholds[threshold],"_rrBlup",sep=""),quote = F,sep="\t",col.names = F,row.names = F);
    write.table(predict_accuracy_pval_iter,paste("predict_accuracy_pval_iter_",traits[trait],"_",thresholds[threshold],"_rrBlup",sep=""),quote = F,sep="\t",col.names = F,row.names = F);
  }
}

