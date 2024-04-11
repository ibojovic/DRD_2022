suppressMessages(library(minfi))


#####################################1##########################################

#setting the working directory where are the data
setwd("C:/Users/itana/Desktop/DRD project/Report-20230109/Input_data_report")
baseDir <- ("C:/Users/itana/Desktop/DRD project/Report-20230109/Input_data_report")
#read the csv file, returns a dataframe with an additional column (Basename). It contains the path of the input data folder to the Idat file
targets <- read.metharray.sheet(baseDir)
#create an object of class RGChannelSet and save it
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")

###################################2 ###########################################

#Creating data frame of red and green channels with usage of functions getGreen and getRed
Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))

################################### 3 ##########################################


#load the manifest
load("C:/Users/itana/Desktop/DRD project/Report-20230109/Input_data_report/Illumina450Manifest_clean.RData")

#my address: 46801437
#check fluorescences for my address
Red[rownames(Red)=="46801437",]
Green[rownames(Green)=="46801437",]

#Data frame for probes of TypeI and Type II
df <- data.frame(getProbeInfo(RGset))
df_II <- data.frame(getProbeInfo(RGset, type = "II"))

#To find type of probe check column addresses for both dataframes 
df[df$AddressA=="46801437",]
df[df$AddressB=="46801437",]
df_II[df_II$AddressA=="46801437",]



################################# 4 ############################################

#assigne fluorescence to methylation signal according to design of microarray
MSet.raw <- preprocessRaw(RGset)


############################### 5 ##############################################


qc <- getQC(MSet.raw)
qc
plotQC(qc)


#controlStripPlot function will plot the fluorescence of different typesof control probes 

controlStripPlot(RGset, controls="NEGATIVE")



# detectionP function creates matrix where for each probe there is corresponding 
#p-value. With deffined threshold of 0.5 for p-value all probes with
#p-value above given threshold will be marked as failed because fluorescence
#in this case cannot be distinguished from the background
detP <- detectionP(RGset)
failed <- detP>0.05

#summary function report samples that have p.value above the threshold
summary(failed)



################################# 6############################################

#getBeta retrives beta values while getM will retrieve M values
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)


#Knowing that mutants samples are 2,5,6,8 and wt are 1,3,4,7 
#2 different data sets need to be created
beta_wt_df <-beta[,c(1,3,4,7)]
M_wt_df <-M[,c(1,3,4,7)]

beta_mt_df <- beta[,c(2,5,6,8)]
M_mt_df <- M[,c(2,5,6,8)]


#To obtain the mean of each row we can use apply() function. With
#'Margin' argument equal to 1 we specify that function mean should be 
# applied to each row of data frame


mean_beta_wt_df <- apply(beta_wt_df,1,mean)
mean_M_wt_df <- apply(M_wt_df,1,mean)
mean_beta_mt_df <- apply(beta_mt_df,1,mean)
mean_M_mt_df <- apply(M_mt_df,1,mean)


#density distribuition
d_mean_beta_wt_df <- density(mean_beta_wt_df,na.rm=T)
d_mean_M_wt_df <- density(mean_M_wt_df,na.rm=T)
d_mean_beta_mt_df <- density(mean_beta_mt_df,na.rm=T)
d_mean_M_mt_df <- density(mean_M_mt_df,na.rm=T)

#Plot of the distribuition of the density function for 
#Beta and M values

par(mfrow=c(1,2))
plot(d_mean_beta_wt_df,main="Density of Beta Values",col="orange")
lines(d_mean_beta_mt_df,col="purple")
plot(d_mean_M_wt_df,main="Density of M Values",col="orange")
lines(d_mean_M_mt_df,col="purple")




######################### 7 #####################################

#preprocessFunnorm

#retrieve Beta values
beta <- getBeta(MSet.raw)

#separate illumina450Manifest_clean into 2 dataframes depending on type I 
#or type II probes; drop unused levels

df_I <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
df_I <- droplevels(df_I)
df_II <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
df_II <- droplevels(df_II)


#Retrain only the rows which name are in 1st column of df_I and df_II
beta_I <- beta[rownames(beta) %in% df_I$IlmnID,]
beta_II <- beta[rownames(beta) %in% df_II$IlmnID,]

#for each probe in beta_I and beta_II, calculate the mean of beta across 8 samples
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)

#calculate the density of the two vectors of mean values
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

#calculate standard deviations densities
sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,na.rm=T)
d_sd_of_beta_II <- density(sd_of_beta_II,na.rm=T)

#perform the normalization using the preprocessFunnorm function on the input object RGset

preprocessFunnorm_results <- preprocessFunnorm(RGset)

#Retrieve beta values

beta_preprocessFunnorm <- getBeta(preprocessFunnorm_results)

#Obtain rows which name are in the first column of the 2 dataframes (from beta matrix)

beta_preprocessFunnorm_I <- beta_preprocessFunnorm[rownames(beta_preprocessFunnorm) %in% df_I$IlmnID,]
beta_preprocessFunnorm_II <- beta_preprocessFunnorm[rownames(beta_preprocessFunnorm) %in% df_II$IlmnID,]


#calculate mean of beta values
mean_of_beta_preprocessFunnorm_I <- apply(beta_preprocessFunnorm_I,1,mean)
mean_of_beta_preprocessFunnorm_II <- apply(beta_preprocessFunnorm_II,1,mean)

#densities of 2 vectors
d_mean_of_beta_preprocessFunnorm_I <- density(mean_of_beta_preprocessFunnorm_I,na.rm=T)
d_mean_of_beta_preprocessFunnorm_II <- density(mean_of_beta_preprocessFunnorm_II,na.rm=T)

#densities of standard deviations (for normalised data)

sd_of_beta_preprocessFunnorm_I <- apply(beta_preprocessFunnorm_I,1,sd)
sd_of_beta_preprocessFunnorm_II <- apply(beta_preprocessFunnorm_II,1,sd)
d_sd_of_beta_preprocessFunnorm_I <- density(sd_of_beta_preprocessFunnorm_I,na.rm=T)
d_sd_of_beta_preprocessFunnorm_II <- density(sd_of_beta_preprocessFunnorm_II,na.rm=T)

# vector to color wt and mt (orange and purple)
color_code <- c('orange','purple','orange','purple','orange','purple','orange','purple')

#plot with 6 panels
par(mfrow=c(2,3))

#plot non-normalised data, with specified x and y to compare plots
plot(d_mean_of_beta_I,col="blue",main="raw beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_II,col="red")

#box plot of the raw data (samplesof mt and wt color determined by vector)
boxplot(beta, col = color_code, ylim=c(0,1))

#same procedure for normalised data
plot(d_mean_of_beta_preprocessFunnorm_I,col="blue",main="preprocessFunnorm beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_preprocessFunnorm_II,col="red")
plot(d_sd_of_beta_preprocessFunnorme_I,col="blue",main="preprocessFunnorm sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_preprocessFunnorm_II,col="red")
boxplot(beta_preprocessFunnorm,col = color_code, ylim=c(0,1))



###################### 8 ####################################
## The matrix has samples in columns and CpG probes in rows. The the prcomp() 
#function has to be applied to the transposed matrix, which is achieved using 
#the t() function 

pca_results <- prcomp(t(beta_preprocessFunnorm),scale =T)

#plot of PC1 and PC2 and label dots, adjust margins and give the legend #Check if samples cluster according to variables 
#cex-size of dots pch-dot type 
levels(samplesheet$Group) 
palette(c("orange","green")) 
plot(pca_results$x[,1],pca_results$x[,2],cex=2,pch=2,col=samplesheet$Group,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000)) 
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1) legend("bottomright",legend=levels(samplesheet$Group),col=c(1:nlevels(samplesheet$Group)),pch=2)

#visualize if samples cluster based on sex (pink-F; blue-M)


levels(samplesheet$Sex) 
palette(c("pink","blue")) 
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=samplesheet$Sex,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000)) text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1) legend("bottomright",legend=levels(samplesheet$Sex),col=c(1:nlevels(samplesheet$Sex)),pch=2)




########################################### 9 ###########################################à


#load and normalize the data 
load("C:/Users/itana/Desktop/DRD project/Report-20230109/Input_data_report/RGset.RData") 
preprocessFunnorm_results <- preprocessFunnorm(RGset)


#retrieve beta values 
beta_preprocessFunnorm <- getBeta(preprocessFunnorm_results)


#read sample sheet that contains samples 
pheno <-read.csv("C:/Users/itana/Desktop/DRD project/Report-20230109/Input_data_report/Samplesheet_report_2022.csv",header=T, stringsAsFactors=T)



#identifying of the CpG probes that are differentlly methylates 
#between samples of WT and MT groups. Parametric test: Mann-Whitney test 
#Function apply will apply the function, and hoc will retrieve p.values 
MannWhitney_function <- function(x) { 
  wilcox <- wilcox.test(x~ pheno$Group) return(wilcox$p.value) 
  }


#apply function to the all rows of beta_processFunnorm 
p_values_MWtest <- apply(beta_preprocessFunnorm,1, MannWhitney_function)


#data frame of all beta values and p_value column 
df_MW_test <- data.frame(beta_preprocessFunnorm, p_values_MWtest)


#ordering in ascending order of probes based on pValues column 
df_MW_test <- df_MW_test[order(df_MW_test $p_values_MWtest),]


#Check how many probes have p value =< 0.05 (differentially methylated) 
df_MW_test_0.05 <- df_MW_test[df_MW_test$p_values_MWtest<=0.05,]


############################### 10 #####################################


#Benjamini & Hochberg and the Bonferroni corrections 
corrected_pValues_BH <- p.adjust(df_MW_test$p_values_MWtest,"BH") 
corrected_pValues_Bonf<- p.adjust(df_MW_test$p_values_MWtest,"bonferroni") 

#data frame of corrected values 
df_MW_test_corrected <- data.frame(df_MW_test, corrected_pValues_BH, corrected_pValues_Bonf) #Visualize the corrected p-values using the box plot colnames(df_MW_test_corrected) boxplot(df_MW_test_corrected[,9:11])


#Visualize the corrected p-values using the box plot 
colnames(df_MW_test_corrected) 
boxplot(df_MW_test_corrected[,9:11])


#Number of probes that survived the multiple test correction

dim(df_MW_test_corrected[df_MW_test_corrected$p_values_MWtest<=0.05,]) #57261 11 
dim(df_MW_test_corrected[df_MW_test_corrected$corrected_pValues_BH<=0.05,]) #0 11 
dim(df_MW_test_corrected[df_MW_test_corrected$corrected_pValues_Bonf<=0.05,]) # 0 11 

#the cut-off becomes very stringent after the corrections



############################ 11 ##########################################

library(qqman) 
#Volcano plot 

#The matrix with the beta values without p.values, from the matrix with corrected p.values 
beta_test <- df_MW_test_corrected[,1:8] 
#difference between the average of Wild type values and the average of mutant values 
#Two matrices with beta values of wild type group and mutant group samples 
#Calculate mean within each group for each row 

beta_test_groupWT <- beta_test[,pheno$Group=="WT"] 
beta_test_groupMUT <- beta_test[,pheno$Group=="MUT"] 
mean_beta_test_groupWT <- apply(beta_test_groupWT,1,mean) 
mean_beta_test_groupMUT <- apply(beta_test_groupMUT,1,mean) 

#calculate difference between average values 
delta_test <- mean_beta_test_groupMUT-mean_beta_test_groupWT 

#data frame of two columns (one with delta values and other with -log10 of p-values) 
#plot 
toVolcPlot <- data.frame(delta_test, -log10(df_MW_test_corrected$p_values_MWtest)) 
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5) 

#highlight probes which have absolute delta (difference between mutant group 
# and wilt type group), higher than certain threshold, (0.1) 
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.05)),] 
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="red")

#Manhattan plot 

# Annotate the dataframe, adding ggenome annotation information for each CpG probe. 
#Illumina450Manifest_clean 
load("C:/Users/itana/Desktop/DRD project/Report-20230109/Input_data_report/Illumina450Manifest_clean.RData") 

#Merge d_MW_test corrected with illumina450Manifest_clean object on common column 
#CpGs probes are stored in rows of the df_MW_test_corrected, transform in column and add to dataframe 
df_MW_test_corrected <- data.frame(rownames(df_MW_test_corrected),df_MW_test_corrected) 
colnames(df_MW_test_corrected)[1] <- "IlmnID" 

#Merge data df_MW_test_corrected_annotated <- merge(df_MW_test_corrected, Illumina450Manifest_clean,by="IlmnID") 

#input for Manhattan plot analysis (probe, chromosome,position on chromosome and the p-value) 
input_Manhattan <- df_MW_test_corrected_annotated[colnames(df_MW_test_corrected_annotated) %in% c("IlmnID","CHR","MAPINFO","p_values_MWtest")] 

#reorder levels of chromosome 
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y") 
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr ) 

#Convert “CHR” to numbers 
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)

#Manhattan plot 
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="p_values_MWtest",annotatePval = 0.00001,col=rainbow(24) ) -log10(0.00001)


######################################## 12 ##############################################


library(gplots) 
#matrix of beta values using top 100 most significant probes 

input_heatmap=as.matrix(df_MW_test[1:100,2:9]) 

# Bar of colors for Group WT(orange) or Group MUT(purple) 
color_heatmap <- c("red","blue","red","red","blue","blue","red","blue") 
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=color_heatmap,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")





