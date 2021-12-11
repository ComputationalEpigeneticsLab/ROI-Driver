##disorder 
#ROT region of total
library(getopt)
library(tidyverse)
library(dbplyr)
rm(list=ls())
gc()
getROI <- function(x){
  # x=  mut.ROT[1,]
  if(as.numeric(x[7])>=as.numeric(x[9])&as.numeric(x[7])<=as.numeric(x[10])){
    a <- x
  }else{
    x[14] <- "non_ROI"
    a <- x
  }
  return(a)
}
getSumLIDR <- function(x){
  x1 <- unique(x[,c(6,7)])
  b = 0
  for (i in 1:dim(x1)[1]) {
    a = as.numeric(x1[i,2])-as.numeric(x1[i,1])+1
    b = b +a
  }
  return(b)
}
getPaste <-  function(x,n){
  if(n==0){
    x <- unique(x)
    x <- gsub(" ","",x)
    a <- c()
    for (i in 1:length(x)) {
      a <- paste(a,x[i],sep = ";")
    }
    return(a)
  }else if(n==1){
    x <- gsub(" ","",x)
    a <- c()
    for (i in 1:length(x)) {
      a <- paste(a,x[i],sep = "|")
    }
    return(a)
  }else if(n==2){
    x <- gsub(" ","",x)
    a <- c()
    for (i in 1:length(x)) {
      a <- paste(a,x[i],sep = ";")
    }
    return(a)
  }else if(n==3){
    x <- gsub(" ","",x)
    a <- c()
    for (i in 1:length(x)) {
      a <- paste(a,x[i],sep = "-")
    }
    return(a)
  }else if(n==4){
    a <- c()
    for (i in 1:dim(x)[1]) {
      x[i,] <- gsub(" ","",x[i,])
      b <- paste(x[i,1],x[i,2],sep = "-")
      a <- rbind(a,b)
    }
    a <- data.frame(region=a)
    return(a)
  }
}
command=matrix(c("Input","i",1,"character",
                 "type","t",1,"character",
                 "output","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
#帮助信息
if (is.null(args$Input) || is.null(args$type) || is.null(args$output) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}
Can <- readr::read_delim(args$Input,delim = "\t")
input_type <- args$type
#input_type="domain"
if(input_type=="disorder"){
  ROT <- read.table("./data/disorder_region.txt",header=T,stringsAsFactors = F,quote = '',sep = '\t')
  colnames(ROT)[1] <- "ENST"
  
  ROT <- unique(ROT[,c(1:6)])
  ROT$ROI_name <- paste0("ROT:",ROT$from,"-",ROT$to)
  colnames(ROT) <- c("ENST","Gene","from","to","LIDR","LWild","ROI_name")
  ROT$ENST <- matrix(unlist(strsplit(ROT$ENST,split = '[.]')),ncol=2,byrow = T)[,1]
  ROT$Gene <- matrix(unlist(strsplit(ROT$Gene,split = '[.]')),ncol=2,byrow = T)[,1]
}else if(input_type=="domain"){
  ROT <- read.table("./data/Input_mutation.txt",header=T,stringsAsFactors = F,quote = '',sep = '\t')
  ROT <- ROT[which(ROT$type=="Domain"&ROT$E.value<0.0001),]
  ENSTWildNum <- matrix(unlist(strsplit(ROT$seq.id,split = "[|]")),ncol = 3,byrow = T) 
  ROT <- data.frame(ENST=ENSTWildNum[,1],Gene=ENSTWildNum[,2],ROT[,c(4:5)],LIDR=ROT[,5]-ROT[,4]+1,LWild=ENSTWildNum [,3],hmm.name=ROT[,7])#envelope.start-end
  ROT$paste <- paste(ROT$ENST,ROT$Gene,ROT$LWild,sep = ';')
  ROT <- unique(ROT[,c(1:7)])#75806
  ROT$ENST <- matrix(unlist(strsplit(ROT$ENST,split = '[.]')),ncol=2,byrow = T)[,1]
  ROT$Gene <- matrix(unlist(strsplit(ROT$Gene,split = '[.]')),ncol=2,byrow = T)[,1]
  colnames(ROT) <- c("ENST","Gene","from","to","LIDR","LWild","ROI_name")
}
#print(head(ROT))
#print(args$type)
#print(args$output)
##########
colnames(ROT)[1] <- "ENST"
wild.frame <- data.frame(ROT,type="ROI")
colnames(wild.frame)[1] <- c("ENST")
Can.mut.sample <- Can[,c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Tumor_Seq_Allele1','Tumor_Seq_Allele2','Tumor_Sample_Barcode','Transcript_ID','Protein_position','HGVSp_Short','Gene')]
Can.mut.sample <- Can.mut.sample[which(Can.mut.sample$Variant_Classification=="Missense_Mutation"&Can.mut.sample$Variant_Type=="SNP"),]
Mut.sample <- data.frame(Hugo_ENST=paste(Can.mut.sample$Gene,Can.mut.sample$Transcript_ID,sep = '_'),Hugo_Symbol=Can.mut.sample$Hugo_Symbol,Gene=Can.mut.sample$Gene,Mutation=paste(Can.mut.sample$Hugo_Symbol,Can.mut.sample$Chromosome,Can.mut.sample$Start_Position,
                                                                                                                                                                                    Can.mut.sample$End_Position,Can.mut.sample$Tumor_Seq_Allele1,Can.mut.sample$Tumor_Seq_Allele2,
                                                                                                                                                                                    Can.mut.sample$Transcript_ID,Can.mut.sample$Protein_position,Can.mut.sample$HGVSp_Short,sep=';'),
                         Sample=Can.mut.sample$Tumor_Sample_Barcode,ENST=Can.mut.sample$Transcript_ID,POS=Can.mut.sample$Protein_position,HGVSp_Short=Can.mut.sample$HGVSp_Short)
mut.ROT <- merge(Mut.sample,wild.frame,by=c("ENST","Gene"),all.x = T)
mut.ROT <- mut.ROT[which(mut.ROT$type=="ROI"),]
if(length(which(mut.ROT$POS=="."))!=0){
  mut.ROT <- mut.ROT[-which(mut.ROT$POS=="."),]
}
mut.ROT1 <- as.data.frame(t(apply(mut.ROT, 1, getROI)),ncol=14,byrow=T,stringsAsFactors = F)
mut.ROT1 <- mut.ROT1[which(mut.ROT1$type=="ROI"),]
mut.ROT1_Symbol_ENST <- unique(mut.ROT1$Hugo_ENST)
IDRData <- c()
for (jj in 1:length(mut.ROT1_Symbol_ENST)) {
  totMut <- unique(Mut.sample[which(Mut.sample$Hugo_ENST==mut.ROT1_Symbol_ENST[jj]),]$Mutation)
  singleSymbol <- unique(mut.ROT1[which(mut.ROT1$Hugo_ENST==mut.ROT1_Symbol_ENST[jj]),])
  singleSymbol$Region <- paste(singleSymbol$from,singleSymbol$to,sep=';')
  singleSymbolRegion <- unique(singleSymbol$Region)
  IDRData1 <- c()
  for (kk in 1:length(singleSymbolRegion)) {
    singleIDR <- singleSymbol[which(singleSymbol$Region==singleSymbolRegion[kk]),]
    #if(dim(singleIDR)[1]>2)print(paste(jj,kk,sep = ';'))
    singleIDRData <- data.frame(Hugo_ENST=mut.ROT1_Symbol_ENST[jj],Hugo_symbol=unique(singleIDR$Hugo_Symbol),Gene=unique(singleIDR$Gene),
                                ENST=unique(singleIDR$ENST),ROI_name=getPaste(singleIDR$ROI_name,1),HGVSp_Short=getPaste(singleIDR$HGVSp_Short,1),POS=getPaste(singleIDR$POS,1),region=getPaste(unique(cbind(singleIDR$from,singleIDR$to)),4),
                                nIDR=length(unique(singleIDR$Mutation)),ntot=length(totMut),LIDR=unique(singleIDR$LIDR),
                                LGene=unique(as.numeric(singleIDR$LWild)),Type="Scatter")#,from=unique(singleIDR$from),to=unique(singleIDR$to)
    #if(dim(singleIDRData)[1]>1&length(unique(singleIDR$Sample))!=1)print(paste(ii,jj,kk,sep = ';'))
    IDRData1 <- rbind(IDRData1,singleIDRData)
    IDRData1$POS <- gsub("^[|]","",IDRData1$POS)
    IDRData1$HGVSp_Short <- gsub("^[|]","",IDRData1$HGVSp_Short)
    #if(dim(IDRData1)[1]>2)print(paste(jj,kk))
  } 
  IDRData1$E <- (as.numeric(IDRData1$nIDR)/as.numeric(IDRData1$ntot))/(as.numeric(IDRData1$LIDR)/as.numeric(IDRData1$LGene))
  IDRData1$P <- 1-pbinom(as.numeric(IDRData1$nIDR),as.numeric(IDRData1$ntot),as.numeric(IDRData1$LIDR)/as.numeric(IDRData1$LGene))
  IDRData <- rbind(IDRData,IDRData1)
  IDRData$POS <- gsub("^;","",IDRData$POS)
  IDRData$region <- gsub("^;","",IDRData$region)
  IDRData$ROI_name <- gsub("^;","",IDRData$ROI_name)
  IDRData$ROI_name <- gsub("^[|]","",IDRData$ROI_name)
  IDRData$HGVSp_Short <- gsub("^;","",IDRData$HGVSp_Short)
}



IDRData$E <- (as.numeric(IDRData$nIDR)/as.numeric(IDRData$ntot))/(as.numeric(IDRData$LIDR)/as.numeric(IDRData$LGene))
IDRData$P <- 1-pbinom(as.numeric(IDRData$nIDR),as.numeric(IDRData$ntot),as.numeric(IDRData$LIDR)/as.numeric(IDRData$LGene))
IDRData <- IDRData[order(IDRData$P,decreasing = F),]
IDRData$fdr <-  p.adjust(as.numeric(IDRData$P),method="BH")#bonferroni之前是BH
data <- IDRData
data$E <- as.numeric(data$E)
data$P <- as.numeric(data$P)
data$nIDR <- as.numeric(data$nIDR)
data <- data %>% filter(nIDR>=3)
data$fdr <- as.numeric(data$fdr)
data$ntot <- as.numeric(data$ntot)
data1 <- data %>% filter(P<0.01,fdr<0.05,E>2)
write.table(data1,args$output,row.names = F,quote = F,sep = '\t')
