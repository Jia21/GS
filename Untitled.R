setwd("~/Dropbox (UNC Charlotte)/poplar_gwas_data/");
id_to_name <- read.table("id_to_name.txt",header = F);
id_to_name <- as.matrix(id_to_name);
sampleid <- read.table("sample_plink_geno",header=F);
id_to_name [,2] = toupper(id_to_name[,2]);
id = matrix(id_to_name[(id_to_name[,2] %in% sampleid[,1]),2],ncol=1);
write.table(id_up,"id_pheno",quote=F,col.names = F,row.names = F,sep="\t");

########
rm(list=ls());
snp <- read.table("Jmax_sig_GT_plink_geno_header",header = F);
snp_sampleid <- unlist(strsplit(t(snp[1,]),"_"));
snp_sampleid <- matrix(unlist(unique(snp_sampleid)),ncol = 1);
colnames(snp) = c("CHROM:POS","CHROM","POS","ID","REF","ALT",snp_sampleid[7:nrow(snp_sampleid),]);
# id_to_name <- read.table("id_to_name.txt",header = F);
# id_to_name <- as.matrix(id_to_name);
# id_to_name_upcase = data.frame(toupper(id_to_name));
# id = id_to_name_upcase[id_to_name_upcase[,3] %in% t(snp[1,4:ncol(snp)]),];
# id = data.frame(id);
# index_snp = as.vector(t(snp[1,4:ncol(snp)]));
# id_sort =  id[match(index_snp, id$V3),];
# write.table(id_sort,"id_sort",quote=F,col.names = F,row.names = F,sep="\t");
# snp_convertid = rbind(cbind(matrix(unlist(snp[1,1:3]),nrow=1),matrix(t(id_sort[,2]),nrow = 1)),snp[2:nrow(snp),]);
# snp_convertid = as.matrix(snp_convertid);

pheno <- read.csv("Photosynthesis_Water_Traits_01.19.18.csv");
pheno = as.matrix(pheno);
sampleid_pheno = pheno[which(!is.na(pheno[,"Jmax25"])),1];
inter_sample = intersect(sampleid_pheno,snp_sampleid[7:nrow(snp_sampleid),1]);
print(length(inter_sample));
# colnames(snp_convertid) = snp_convertid[1,];
rownames(pheno) = pheno[,1];
pheno_specific = pheno[inter_sample,"Jmax25",drop=FALSE];
pheno_GT = snp[2:nrow(snp),c("CHROM:POS","CHROM","POS","ID","REF","ALT",inter_sample)];
write.table(pheno_specific,"pheno_Av_Diameter_mm",quote = F,sep="\t",row.names = T,col.names = T);

miss = matrix("NA",nrow(pheno_GT),1);
for(sl in 1:nrow(pheno_GT)){
  slice = pheno_GT[sl,7:ncol(pheno_GT)];
  miss[sl] = sum(slice == "./.");
}
fre_miss = matrix(as.numeric(miss)/ncol(pheno_GT),ncol=1);
pheno_GT_filter_missGT = pheno_GT[which(fre_miss < 0.20),];
print(nrow(pheno_GT_filter_missGT));

pheno_GT_filter_missGT <- as.matrix(pheno_GT_filter_missGT);
maf = matrix("NA",nrow = nrow(pheno_GT_filter_missGT),1);
#########only one genotype????
for(sl in 1:nrow(pheno_GT_filter_missGT)){
  slice = pheno_GT_filter_missGT[sl,7:ncol(pheno_GT_filter_missGT),drop=FALSE];
  tab = data.frame(table(slice));
 # if(dim(tab)[[1]] != 1){
    tmp1 = (tab[which(tab[,1] == "0/0"),2] *2 + tab[which(tab[,1] == "0/1"),2]) /(sum(pheno_GT_filter_missGT[sl,] != "./.") *2);
    tmp2 = (tab[which(tab[,1] == "1/1"),2] *2 + tab[which(tab[,1] == "0/1"),2]) /(sum(pheno_GT_filter_missGT[sl,] != "./.") *2);
  #}else{
  #   tmp1 = 0;
  #   tmp2 = 1;
  # }
  maf[sl] = min(tmp1,tmp2);
}
##########??????
pheno_GT_filter_missGT_MAF = pheno_GT_filter_missGT[maf > 0.05,];
print(nrow(pheno_GT_filter_missGT_MAF));
write.table(pheno_GT_filter_missGT_MAF,"Av_Diameter_mm.geno",quote = F,sep = "\t",col.names = T,row.names = F)


########New plink genotype data:
tped <- read.table("temp020_12.tped",header = F);
tped <- as.matrix(tped);
tfam <- read.table("temp020_12.tfam",header = F);
tfam <- as.matrix(tfam);

genotype <- tped[,5:ncol(tped)];
genotype <- as.matrix(genotype);

chr_info <- cbind(matrix(tped[,1],ncol=1),matrix(tped[,2],ncol=1),matrix(tped[,4],ncol=1));

vcf <- NULL;
for(i in 1:nrow(genotype)){
  print(i);
  sl =genotype[i,,drop=F];
  num_1 = sum(sl == 2);
  num_2 = sum(sl == 1);
  chr = paste("scaffold","_",chr_info[i,1],sep="");
  
  if(num_1 > num_2) {
    ref_allele = 2;
    alt_allele = 1;
    sl[,which(sl[1,] == ref_allele)] = 0;
  }
  if(num_1 < num_2) {
    ref_allele = 1;
    alt_allele = 2;
    sl[,which(sl[1,] == ref_allele)] = 0;
    sl[,which(sl[1,] != ref_allele && sl[1,] != 0)] = 1;
  }
  
  geno_single <- NULL;
  geno_sample <- NULL;
  for(j in 1:nrow(tfam)){
    sampleID = paste(tfam[j,1],"_",tfam[j,2],sep="");
    geno_sample <- cbind(geno_sample,sampleID);
    sampleID_geno = paste(sl[,j*2-1],"/",sl[,j*2],sep="");
    geno_single <- cbind(geno_single,sampleID_geno);
  }
  geno_single <- cbind(chr,tped[i,4],tped[i,2],ref_allele,alt_allele,geno_single);
  colnames(geno_single) <- c("#CHROM","POS","ID","REF","ATL",geno_sample);
  
  vcf = rbind(vcf,geno_single);
}

write.table(vcf,"20.VCF.CORMAT.vcf",quote = F,sep="\t",col.names = T,row.names = F);

####Set up different threshold of p.values####
setwd("~/Dropbox (UNC Charlotte)/poplar_gwas_data/GS_data-03272018_Wellington/geno_pheno/");
sig_gwas_loci <- NULL;
threshold = c(1.00E-05,1.00E-04,1.00E-03,1.00E-02,5.00E-02);
trait = c("Jmax25","Rdlight25","Resistwp25","WUEref","Av_Diameter_mm");

for(j in 1:length(trait)){
  filename = paste(trait[j],".txt",sep="");
  data <- read.table(filename,header=F);
  for (i in 1:length(threshold)){
    data_sub <- data[which(data[,3] < threshold[i]),];
    newID <- matrix(".",nrow = nrow(data_sub),1);
    newID[which(data_sub[,1] < 10),] <- as.matrix(unlist(paste("Chr0",data_sub[which(data_sub[,1] < 10),1],":",data_sub[which(data_sub[,1] < 10),2],sep="")),ncol=1);
    newID[which(data_sub[,1] >= 10 & data_sub[,1] < 20),] <- as.matrix(unlist(paste("Chr",data_sub[which(data_sub[,1] >= 10 & data_sub[,1] < 20),1],":",data_sub[which(data_sub[,1] >= 10 & data_sub[,1] < 20),2],sep="")),ncol=1);
    newID[which(data_sub[,1] >= 20),] <- as.matrix(unlist(paste("scaffold_",data_sub[which(data_sub[,1] >= 20),1],":",data_sub[which(data_sub[,1] >= 20),2],sep="")),ncol=1);
    
    data_geno_GS <- cbind(newID,data_sub);
    
    write.table(data_geno_GS,paste("./geno_pheno/",trait[j],"_",threshold[i],"_","geno",sep=""),quote=F,sep = "\t",col.names = F,row.names = F)
  }
}


#####Generate genotype and phenotype for GS:
setwd("~/Dropbox (UNC Charlotte)/poplar_gwas_data/GS_data-03272018_Wellington/geno_pheno/");
rm(list=ls());
thresholds = c(1e-05,1e-04,0.001);
traits = c("Jmax25","Rdlight25","Resistwp25","WUEref","Av_Diameter_mm");
for (i in 1:length(thresholds)){
  for(j in 1:length(traits)){
    cat("i is",i,";j is",j,"\n");
    snp <- read.table(paste(traits[j],"_",thresholds[i],"_geno_plink.header",sep=""),header = F);
    snp_sampleid <- unlist(strsplit(t(snp[1,]),"_"));
    snp_sampleid <- matrix(unlist(unique(snp_sampleid)),ncol = 1);
    colnames(snp) = c("CHROM:POS","CHROM","POS","ID","REF","ALT",snp_sampleid[7:nrow(snp_sampleid),]);
    
    pheno <- read.csv("../../Photosynthesis_Water_Traits_01.19.18.csv");
    pheno = as.matrix(pheno);
    sampleid_pheno = pheno[which(!is.na(pheno[,traits[j]])),1];
    inter_sample = intersect(sampleid_pheno,snp_sampleid[7:nrow(snp_sampleid),1]);
    print(length(inter_sample));
    # colnames(snp_convertid) = snp_convertid[1,];
    rownames(pheno) = pheno[,1];
    pheno_specific = pheno[inter_sample,traits[j],drop=FALSE];
    pheno_GT = snp[2:nrow(snp),c("CHROM:POS","CHROM","POS","ID","REF","ALT",inter_sample)];
    write.table(pheno_specific,paste("pheno_",traits[j],"_",thresholds[i],sep=""),quote = F,sep="\t",row.names = T,col.names = T);
    
    
    miss = matrix("NA",nrow(pheno_GT),1);
    for(sl in 1:nrow(pheno_GT)){
      #print(sl);
      slice = pheno_GT[sl,7:ncol(pheno_GT)];
      miss[sl] = sum(slice == "./.");
    }
    fre_miss = matrix(as.numeric(miss)/ncol(pheno_GT),ncol=1);
    pheno_GT_filter_missGT = pheno_GT[which(fre_miss < 0.20),];
    print(nrow(pheno_GT_filter_missGT));
    
    pheno_GT_filter_missGT <- as.matrix(pheno_GT_filter_missGT);
    maf = matrix("NA",nrow = nrow(pheno_GT_filter_missGT),1);
    
    for(sl in 1:nrow(pheno_GT_filter_missGT)){
      #print(sl);
      slice = pheno_GT_filter_missGT[sl,7:ncol(pheno_GT_filter_missGT),drop=FALSE];
      tab = data.frame(table(slice));
      if(dim(tab)[[1]] != 1){
        a1 = tab[which(tab[,1] == "0/0"),2] *2;
        b1 = tab[which(tab[,1] == "0/1" | tab[,1] == "1/0"),2];
        c1 = tab[which(tab[,1] == "1/1"),2] *2;
        if(identical(b1, integer(0))){
          tmp1 = a1/(sum(pheno_GT_filter_missGT[sl,7:ncol(pheno_GT_filter_missGT)] != "./.") *2);
          tmp2 = c1/(sum(pheno_GT_filter_missGT[sl,7:ncol(pheno_GT_filter_missGT)] != "./.") *2);
        }else{
          tmp1 = (a1 + b1) /(sum(pheno_GT_filter_missGT[sl,7:ncol(pheno_GT_filter_missGT)] != "./.") *2);
          tmp2 = (c1 + b1) /(sum(pheno_GT_filter_missGT[sl,7:ncol(pheno_GT_filter_missGT)] != "./.") *2);
        }
      }else{
         tmp1 = 0;
         tmp2 = 1;
      }
      maf[sl] = min(tmp1,tmp2);
    }
    ##########If one SNP only have heteruzyguos for all samples, remove it!
    pheno_GT_filter_missGT_MAF = pheno_GT_filter_missGT[maf > 0.05,];
    print(nrow(pheno_GT_filter_missGT_MAF));
    write.table(pheno_GT_filter_missGT_MAF,paste(traits[j],"_",thresholds[i],".geno",sep=""),quote = F,sep = "\t",col.names = T,row.names = F);
  }
}
