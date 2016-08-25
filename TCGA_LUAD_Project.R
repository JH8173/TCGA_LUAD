
######################
# Read clinical data #
######################

read_clinical <- function(clinical_file = "nationwidechildrens.org_clinical_patient_luad.txt"){
  clinical <- read.csv(clinical_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  clinical <- clinical[3:nrow(clinical),]
  Barcode <- clinical$bcr_patient_barcode
  Gender <- clinical$gender
  Age <-  suppressWarnings(as.numeric(clinical$birth_days_to)*-1)
  Age <- factor(as.integer(Age / 3652.5))
  Status <- factor(clinical$vital_status)
  Contact <- clinical$last_contact_days_to
  Death <- clinical$death_days_to
  Death[which(Status=="Alive")] <- Contact[which(Status=="Alive")]
  Death[which(Death=="[Not Available]")] <- NA
  Time <- suppressWarnings(as.numeric(Death))+1
  
  # Integrate tumor stages into big class level
  clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage=="[Not Available]")] <- NA
  clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage=="Stage IA")] = "Stage I"
  clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage=="Stage IB")] = "Stage I"
  clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage=="Stage IIA")] = "Stage II"
  clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage=="Stage IIB")] = "Stage II"
  clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage=="Stage IIIA")] = "Stage III"
  clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage=="Stage IIIB")] = "Stage III"
  Stage <- factor(clinical$ajcc_pathologic_tumor_stage)
  clinical$tobacco_smoking_history_indicator[which(clinical$tobacco_smoking_history_indicator=="[Not Available]")] <- NA
  clinical$tobacco_smoking_history_indicator[which(clinical$tobacco_smoking_history_indicator=="[Unknown]")] <- NA
  clinical$tobacco_smoking_history_indicator[which(clinical$tobacco_smoking_history_indicator=="Current Reformed Smoker, Duration Not Specified")] <- NA
  Tobacco <- factor(clinical$tobacco_smoking_history_indicator)
  summary(factor(clinical$tobacco_smoking_history_indicator))
  clinical_factor <- data.frame(Barcode, Age, Gender, Status, Time, Stage ,Tobacco)
  return(clinical_factor)
}
LUAD_clinical <- read_clinical("nationwidechildrens.org_clinical_patient_luad.txt")
LUSC_clinical <- read_clinical("nationwidechildrens.org_clinical_patient_lusc.txt")
View(LUAD_clinical)
View(LUSC_clinical)
LUAD_clinical[1:10, 1:10]
LUAD_clinical[1:10,]

####################
# Read RNASeq data #
####################
read_rnaseq <- function(rnaseq_file="LUAD_rnaseq"){
  rnaseq <- read.table(rnaseq_file, sep="\t", header=TRUE, row.names = 1)
  tumor_idx <- which(substr(colnames(rnaseq), 14, 16) == "01A") # Only tumor samples
  colnames(rnaseq) <- substr(gsub("\\.", "-", colnames(rnaseq)), 1, 12)
  return( as.matrix( rnaseq[, tumor_idx]) )
}
LUAD_rnaseq <- read_rnaseq("LUAD_rnaseq")
LUAD_rnaseq <- read_rnaseq("./LUAD_rnaseq/LUAD_rnaseq")
LUAD_rnaseq[1:10, 1:10]
LUSC_rnaseq <- read_rnaseq("./LUSC_rnaseq")
LUSC_rnaseq <- read_rnaseq("./LUAD_rnaseq/LUSC_rnaseq")
View(LUSC_rnaseq)


###############################
# Intersect clinical & RNASeq #
###############################

intersect_data <- function(rnaseq, clinical){
  intersect_Barcode <- intersect( colnames(rnaseq), clinical$Barcode)
  intersect_rnaseq <- c() # clinical Barcode
  intersect_clinical <- c() # clinical Barcode
  for(i in 1:length(intersect_Barcode)){
    intersect_clinical[i] <- which(clinical$Barcode == intersect_Barcode[i])
    intersect_rnaseq[i] <- which(colnames(rnaseq) == intersect_Barcode[i])
  }
  return(list(intersect_rnaseq, intersect_clinical))
} #List 1: rnaseq index, List2: clinical index
LUAD_idx <- intersect_data(LUAD_rnaseq, LUAD_clinical)
LUAD_idx
View(LUAD_clinical)
LUAD_rnaseq <- LUAD_rnaseq[, LUAD_idx[[1]] ]
LUAD_clinical <- LUAD_clinical[LUAD_idx[[2]], ]
LUSC_idx <- intersect_data(LUSC_rnaseq, LUSC_clinical)
LUSC_rnaseq <- LUAD_rnaseq[, LUSC_idx[[1]] ]
LUSC_clinical <- LUAD_clinical[LUSC_idx[[2]], ]


######################
# Survival analysis. #
######################
#Independent variables: One gene expression, Age, Gender, TumorStage, Tobacco, ITGA2
#dependent variable: Survival time


library(survival)
i=1
LUAD_pvalue <- c()
LUSC_pvalue <- c()
dim(LUSC_rnaseq)
for(i in 1:nrow(LUAD_rnaseq)){
  out_LUAD <- coxph(Surv(as.numeric(Time),Status=="Dead")~LUAD_rnaseq[i,]+Age+Gender+Stage+Tobacco, data=LUAD_clinical)
  out_LUSC <- coxph(Surv(as.numeric(Time),Status=="Dead")~LUSC_rnaseq[i,]+Age+Gender+Stage+Tobacco, data=LUSC_clinical)
  LUAD_pvalue[i] <- coef(summary(out_LUAD))[1,5]
  LUSC_pvalue[i] <- coef(summary(out_LUSC))[1,5]
}
LUAD_pvalue<0.05
View(LUAD_clinical)
warnings()
dim(LUSC_rnaseq)
sum(is.na(LUAD_pvalue))
LUAD_pvalue_bonferroni <- p.adjust(LUAD_pvalue, method="bonferroni")
LUSC_pvalue_bonferroni <- p.adjust(LUSC_pvalue, method="bonferroni")
LUAD_genes_associated_with_Survival <- rownames(LUAD_rnaseq)[which(LUAD_pvalue_bonferroni<0.05)]
LUSC_genes_associated_with_Survival <- rownames(LUSC_rnaseq)[which(LUSC_pvalue_bonferroni<0.05)]
LUSC_genes_associated_with_Survival
LUAD_genes_associated_with_Survival

##################################
Gender로 kaplan meier graph 그리기
##################################

install.packages("survival", repos="http://cran.r-project.org" )
require(survival)

f_CXCL17 <- factor(CXCL17 > median(CXCL17),
                   levels = c('FALSE', 'TRUE'),
                   labels = c('Low', 'High'))
LUAD_clinical$CXCL17 <- f_CXCL17
fit <- survfit(Surv(Time, Status=="Dead") ~ CXCL17, 
               data=subset(LUAD_clinical, Stage != "Stage IV")) 
plot(fit, lty=2:3, xlab="CXCL17", ylab="Survival") 
legend("topright", legend=c('Low', 'High'), lty = c(2,3))

#############
#Log-rank test 
#############
survdiff(Surv(Time,Status=="Dead")~f_CXCL17, data=LUAD_clinical)

######################
# Correlation analysis. #
######################

coeffi_cor_SC <- c()
?cor
coeffi_cor_SC <- c()
for (i in 1: dim(LUSC_rnaseq)){
  coeffi_cor_SC[i]  <- cor(LUSC_rnaseq[,"RAET1G"], LUSC[i,])
}
dim(LUSC_rnaseq)
LUSC_rnaseq[,"RAET1G"]
LUSC_rnaseq["RAET1G",]
for (i in 1: dim(LUSC_rnaseq)){
  coeffi_cor_SC[i]  <- cor(LUSC_rnaseq["RAET1G",], LUSC[i,])
}
for (i in 1: dim(LUSC_rnaseq)[1]){
  coeffi_cor_SC[i]  <- cor(LUSC_rnaseq["RAET1G",], LUSC[i,])
}
for (i in 1: dim(LUSC_rnaseq)[1]){
  coeffi_cor_SC[i]  <- cor(LUSC_rnaseq["RAET1G",], LUSC_rnaseq[i,])
}
rownames(LUSC_rnaseq)[coeffi_cor_SC > 0.6]
rownames(LUSC_rnaseq)
coeffi_cor_SC > 0.6
rownames(LUSC_rnaseq)[coeffi_cor_SC > 0.3]

sum(names(CXCL17) == LUAD_clinical$Barcode)
