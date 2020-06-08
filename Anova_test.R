#Rscript

#@Author: Yunlong Ma
#ANOVA test for three group: Tuberculosis, healthy infection, healthy control
#The RNA expression data were downloaded from the GEO database (Accession # GSE133803)

if(!require("psych")) install.packages("psych")
library(psych)

#set the working directory for the Rscript dependent on user their need
setwd("C:\\Users\\Administrator\\Desktop\\anova")
getwd()

#Read data from the working directory
Data <- read.table("GSE139825.txt",header = TRUE)
class(Data)

#Arranged the data for downstream data analysis
anova_matrix<-  Data[,c(-1,-2)]

#Calculate the mean and standard variance for each group
result_mean <-aggregate(Data$CDC16,by = list(Data$Pheno),FUN= mean)
result_sd <-  aggregate(Data$CDC16,by = list(Data$Pheno),FUN= sd)

#Defined a function for calculating Anova P values and Fold changes
Anova_test <- function(x){
  fit <- aov(x~Data$Pheno)
  A_test <- summary(fit)
  Anova_Pvalue <- A_test[[1]][1,5] 
  multi_p <- TukeyHSD(fit)
  p_data <- as.matrix(multi_p$`Data$Pheno`[,4])
  row_names <- c("Anova_Pvalue", rownames(p_data),paste(rownames(p_data),"fold")) 
  Test_summary <- rbind(Anova_Pvalue, p_data)
  #Use the function of describe.by() in psych package to calculate statisitcs for each group
  des_data <- describeBy(x,list(Data$Pheno)) 
  length_des<- length(des_data)
  mean_data <- c()
  
  #Calculate the fold change
  for ( i in 1:length_des){
    b<- i+1
    if (i == length_des ){
      break
    }
    for ( j in b:length_des ){
      mean_data<- c(mean_data,as.data.frame(des_data[j])[3]/as.data.frame(des_data[i])[3])
    }
  }
  mean_data<-t(as.data.frame(mean_data))
  Test_summary <- rbind(Anova_Pvalue,p_data,mean_data)
  return(Test_summary)
}

#Apply the defined function of Anova_test to all the data matrix
test_res <- as.data.frame(apply(anova_matrix, 2, Anova_test))
rownames(test_res) <- row_names
test_res <- t(test_res)

#Write out the result of the analysis
write.csv(test_res, file="anova_test_results.csv")


#Ploting the figures by the function of boxplot in R platform
boxplot(RPS5 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(ZNF197 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(FCHO1 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(CLN8_1 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(SPATA20 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(PDK1 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(ZNF354A ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(TDRKH ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(CDK10_2 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(LIG3 ~ Pheno, data = Data,col=c("green","green","red","red"))
boxplot(CDC16_1 ~ Pheno, data = Data,col=c("green","green","red","red"))


#End

