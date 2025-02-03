library(data.table)
library(stringr)
library(RobustRankAggreg)
library(plyr)
############Union Pro_Phos_RNA and cal RRA for each project############
setwd('/hwdata/home/wuyj2/eKPI/spearman/')
sample_Thrd <- 6

filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','PMID','TNType','Omics')] <- str_split_fixed(cancerFrame$fileName,'_',5)
cancerFrame$Omics <- str_replace_all(cancerFrame$Omics,'Kinase.','')
cancerFrame$ID <- paste(cancerFrame$CancerName,cancerFrame$TNType,sep='_')

kinasedf <- read.csv('/hwdata/home/wuyj2/eKPI/CancerOmics/kinase_basic(Sheet1).csv',header = T)
kinaseNames <- union(kinasedf$Offical_gene_symbol,kinasedf$KinBase_name)[-1]

Sample_THRD <- 6
Pvalue_THRD <- 0.05
Rvalue_THRD <- 0.2

print('Union Pro_Phos_RNA and Prediction Tool for each project V2')

exprimentDF <- readRDS('/hwdata/home/wuyj2/eKPI/Prediction_Tool/Expriment.clean.rds')
exprimentDF$ID <- paste(exprimentDF$Kinase_names,str_split_fixed(exprimentDF$Pro.site,'#',2)[,2],sep='#')

ScansiteDF <- readRDS('/hwdata/home/wuyj2/eKPI/Prediction_Tool/Scansite.clean.rds')
ScansiteDF$ID <- paste(ScansiteDF$Kinase_names,str_split_fixed(ScansiteDF$Pro.site,'#',2)[,2],sep='#')

# MusiteDeepDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/MusiteDeep.clean.rds')
# MusiteDeepDF$ID <- paste(MusiteDeepDF$Kinase_names,MusiteDeepDF$Pro.site,sep='#')
# print('MusiteDeepDF')
#
# NetworKINDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/NetworKIN.clean.rds')
# NetworKINDF$ID <- paste(NetworKINDF$Kinase_names,NetworKINDF$Pro.site,sep='#')
# print('NetworKINDF')
#
# GPShighDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_High.rds')
# GPShighDF$ID <- paste(GPShighDF$kinaseName,GPShighDF$Pro.Site,sep='#')
# print('GPShighDF')
#
# GPSmedianDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_Medium.rds')
# GPSmedianDF$ID <- paste(GPSmedianDF$kinaseName,GPSmedianDF$Pro.Site,sep='#')
# print('GPSmedianDF')

ProjectNames <- unique(cancerFrame$ID)
for(PjName in ProjectNames){
  print(PjName)
  myUnionDF <- data.frame()
  fileDataFrame <- cancerFrame[cancerFrame$ID==PjName,]
  write_file <- paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Spearman_RRAV2/',PjName,'.tsv',sep='')
  if(nrow(fileDataFrame)>0 & !file.exists(write_file)){
    if('Phos' %in% fileDataFrame$Omics){
      filedata <- fread(paste(fileDataFrame[fileDataFrame$Omics=='Phos','fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[!is.na(filedata$Rho),]
      filedata <- filedata[filedata$Rho>0&filedata$sample.num>=Sample_THRD,]

      filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
      filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$FDR <- p.adjust(filedata$Rho.pvalue,method='BH')

      max_sample <- max(c(217,max(filedata$sample.num)))
      k1 <- max_sample^3-max_sample
      filedata$p.adjust <- ((k1-(filedata$sample.num^3-filedata$sample.num))*filedata$Rho.pvalue+(filedata$sample.num^3-filedata$sample.num)*filedata$FDR)/k1


      filedata <- filedata[filedata$sample.num>=Sample_THRD&
                             filedata$Rho>=Rvalue_THRD&
                             filedata$p.adjust<=Pvalue_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair
      if(nrow(filedata)>0){
        filedata$ID  <- paste(filedata$KinaseName,str_split_fixed(filedata$GeneSite,'#',2)[,2],sep='#')
        filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
        myUnionDF <- data.frame(ID = filedata$Group.1,
                                Phos.Rho = filedata$x)
        myUnionDF <- myUnionDF[order(myUnionDF$Phos.Rho,decreasing=T),]
        myUnionDF$Phos.Rank <- nrow(myUnionDF) - rank(myUnionDF$Phos.Rho,ties.method = 'max') +1##相同spearman值，取平均秩
      }
      print('Phos')
    }
    if('Pro' %in% fileDataFrame$Omics){
      filedata <- fread(paste(fileDataFrame[fileDataFrame$Omics=='Pro','fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[!is.na(filedata$Rho),]
      filedata <- filedata[filedata$Rho>0&filedata$sample.num>=Sample_THRD,]

      filedata <- filedata[which(filedata$RNA.names %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$RNA.names!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$FDR <- p.adjust(filedata$Rho.pvalue,method='BH')

      max_sample <- max(c(217,max(filedata$sample.num)))
      k1 <- max_sample^3-max_sample
      filedata$p.adjust <- ((k1-(filedata$sample.num^3-filedata$sample.num))*filedata$Rho.pvalue+(filedata$sample.num^3-filedata$sample.num)*filedata$FDR)/k1


      filedata <- filedata[filedata$sample.num>=Sample_THRD&
                             filedata$Rho>=Rvalue_THRD&
                             filedata$p.adjust<=Pvalue_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair
      if(nrow(filedata)>0){
        filedata$ID  <- paste(filedata$RNA.names,str_split_fixed(filedata$GeneSite,'#',2)[,2],sep='#')
        length(unique(filedata$ID)) == nrow(filedata)
        tmpDF <- data.frame(ID = filedata$ID,
                            Pro.Rho = filedata$Rho)
        tmpDF <- tmpDF[order(tmpDF$Pro.Rho,decreasing=T),]
        tmpDF$Pro.Rank <- nrow(tmpDF) - rank(tmpDF$Pro.Rho,ties.method = 'max')+1
        myUnionDF <- merge(myUnionDF,tmpDF,by='ID',all=T)
      }
      print('Pro')
    }
    if('RNA' %in% fileDataFrame$Omics){
      filedata <- fread(paste(fileDataFrame[fileDataFrame$Omics=='RNA','fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[!is.na(filedata$Rho),]
      filedata <- filedata[filedata$Rho>0&filedata$sample.num>=Sample_THRD,]

      filedata <- filedata[which(filedata$RNA.names %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$RNA.names!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$FDR <- p.adjust(filedata$Rho.pvalue,method='BH')

      max_sample <- max(c(217,max(filedata$sample.num)))
      k1 <- max_sample^3-max_sample
      filedata$p.adjust <- ((k1-(filedata$sample.num^3-filedata$sample.num))*filedata$Rho.pvalue+(filedata$sample.num^3-filedata$sample.num)*filedata$FDR)/k1


      filedata <- filedata[filedata$sample.num>=Sample_THRD&
                             filedata$Rho>=Rvalue_THRD&
                             filedata$p.adjust<=Pvalue_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair
      if(nrow(filedata)>0){
        filedata$ID  <- paste(filedata$RNA.names,str_split_fixed(filedata$GeneSite,'#',2)[,2],sep='#')
        length(unique(filedata$ID)) == nrow(filedata)
        tmpDF <- data.frame(ID = filedata$ID,
                            RNA.Rho = filedata$Rho)
        tmpDF <- tmpDF[order(tmpDF$RNA.Rho,decreasing=T),]
        tmpDF$RNA.Rank <- nrow(tmpDF) - rank(tmpDF$RNA.Rho,ties.method = 'max')+1
        myUnionDF <- merge(myUnionDF,tmpDF,by='ID',all=T)
      }
      print('RNA')
    }
    if(nrow(myUnionDF)>0){
      rownames(myUnionDF) <- myUnionDF$ID
      if(ncol(myUnionDF)==5){
        # myUnionDF[is.na(myUnionDF[,3]),3] <- nrow(myUnionDF)
        # myUnionDF[is.na(myUnionDF[,5]),5] <- nrow(myUnionDF)
        # RRA_rst <- aggregateRanks(rmat = as.matrix(myUnionDF[,c('Phos.Rank','Pro.Rank')]), N = length(myUnionDF$ID),method = 'RRA')
        colIndex <- c(3,5)
        project_list <- list()
        for (i in colIndex){
          tmpdata <- myUnionDF[!is.na(myUnionDF[,i]),]
          tmpdata <- tmpdata[order(tmpdata[,colnames(myUnionDF)[i]]),]
          project_list[[colnames(myUnionDF)[i]]] <- rownames(tmpdata)#as.character(unlist(plyr::arrange(na.omit(myUnionDF[!is.na(myUnionDF[,i]),c('ID',colnames(myUnionDF)[i])]),colnames(myUnionDF)[i])[1]))
        }
        RRA_rst <- aggregateRanks(project_list,method = 'RRA')
      }
      if(ncol(myUnionDF)==7){
        # myUnionDF[is.na(myUnionDF$Phos.Rank),]$Phos.Rank <- nrow(myUnionDF)
        # myUnionDF[is.na(myUnionDF$Pro.Rank),]$Pro.Rank <- nrow(myUnionDF)
        # myUnionDF[is.na(myUnionDF$RNA.Rank),]$RNA.Rank <- nrow(myUnionDF)
        # RRA_rst <- aggregateRanks(rmat = as.matrix(myUnionDF[,c('Phos.Rank','Pro.Rank','RNA.Rank')]), N = length(myUnionDF$ID),method = 'RRA')
        colIndex <- c(3,5,7)
        project_list <- list()
        for (i in colIndex){
          print(nrow(myUnionDF))
          print(nrow(myUnionDF[!is.na(myUnionDF[,i]),]))
          tmpdata <- myUnionDF[!is.na(myUnionDF[,i]),]
          tmpdata <- tmpdata[order(tmpdata[,colnames(myUnionDF)[i]]),]
          project_list[[colnames(myUnionDF)[i]]] <- rownames(tmpdata)#as.character(unlist(plyr::arrange(myUnionDF[!is.na(myUnionDF[,i]),c('ID',colnames(myUnionDF)[i])],colnames(myUnionDF)[i])[1]))
        }
        RRA_rst <- aggregateRanks(project_list,method = 'RRA')

      }
      # print('order before')
      myUnionDF <- merge(myUnionDF,RRA_rst,by.x='ID',by.y='Name')
      myUnionDF <- myUnionDF[order(myUnionDF$Score,decreasing=T),]
      myUnionDF$RRA.Rank <- rank(myUnionDF$Score,ties.method = 'max')
      # print('order after')
      myUnionDF$expriment <- FALSE
      myUnionDF[which(myUnionDF$ID %in% exprimentDF$ID),]$expriment <- TRUE

      myUnionDF$Scansite <- FALSE
      myUnionDF[which(myUnionDF$ID %in% ScansiteDF$ID),]$Scansite <- TRUE
      if(!file.exists(write_file)){
        write.table(myUnionDF,write_file,quote=F,row.names=F,sep='\t')
      }
    }
  }
}


############plot Figure2A#################
setwd('/hwdata//home/wuyj2/eKPI/ArticleCancerStat/Spearman_RRAV2/')
filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','TNType')] <- str_split_fixed(str_replace_all(cancerFrame$fileName,'.tsv',''),'_',2)

tmpdata <- data.frame(table(cancerFrame$CancerName))
cancerFrame <- cancerFrame[which(cancerFrame$CancerName %in% tmpdata[tmpdata$Freq==2,]$Var1),]

finaldata <- data.frame()
for(i in 1:nrow(cancerFrame)){
  filedata <- fread(cancerFrame[i,'fileName'],header = T,sep='\t',data.table = F)
  if(ncol(filedata)==11){
    filedata$Phos.Rank <- 1 - filedata$Phos.Rank/length(na.omit(filedata$Phos.Rank))
    filedata$Pro.Rank <- 1 - filedata$Pro.Rank/length(na.omit(filedata$Pro.Rank))
    filedata$RNA.Rank <- 1 - filedata$RNA.Rank/length(na.omit(filedata$RNA.Rank))
    filedata$RRA.Rank <- rank(filedata$Score,ties.method = 'random' )
    filedata$RRA.Rank <- 1 - filedata$RRA.Rank/length(na.omit(filedata$RRA.Rank))
    filedata$Freq <- rowSums(!is.na(filedata[,c('Phos.Rank','Pro.Rank','RNA.Rank')]))

    filedata <- filedata[filedata$expriment==TRUE,]
    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            RANKvalue = na.omit(filedata$Phos.Rank),
                                            OmicName = 'Phos'))
    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            RANKvalue = na.omit(filedata$Pro.Rank),
                                            OmicName = 'Pro'))
    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            RANKvalue = na.omit(filedata$RNA.Rank),
                                            OmicName = 'RNA'))
    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            RANKvalue = na.omit(filedata$RRA.Rank),
                                            OmicName = 'RRA'))
  }
  if(ncol(filedata)==9){
    filedata$Phos.Rank <- 1 - filedata$Phos.Rank/length(na.omit(filedata$Phos.Rank))
    filedata$Pro.Rank <- 1 - filedata$Pro.Rank/length(na.omit(filedata$Pro.Rank))
    filedata$RRA.Rank <- rank(filedata$Score,ties.method = 'random')
    filedata$RRA.Rank <- 1 - filedata$RRA.Rank/length(na.omit(filedata$RRA.Rank))
    filedata <- filedata[filedata$expriment==TRUE,]

    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            RANKvalue = na.omit(filedata$Phos.Rank),
                                            OmicName = 'Phos'))
    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            RANKvalue = na.omit(filedata$Pro.Rank),
                                            OmicName = 'Pro'))
    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            RANKvalue = na.omit(filedata$RRA.Rank),
                                            OmicName = 'RRA'))
  }
}
write.table(finaldata,'/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/Figure2A_pvalue_RRA_rank.tsv',quote=F,row.names=F,sep='\t')

finaldata$TNType <- factor(finaldata$TNType,levels =c('Normal','Tumor'))
finaldata$OmicName <- factor(finaldata$OmicName,levels =c('Pro','Phos','RNA','RRA'))

p <- ggplot(finaldata,aes(x=OmicName,y=RANKvalue,fill=TNType))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif",method="wilcox.test",size=6)+
  labs(x='',y='Normalized Rank')+
  scale_y_continuous(expand = c(0.01,0.1))+
  scale_fill_manual(values=c('#455a82','#7c3031'))+
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        legend.box = element_blank(),
        legend.background = element_blank(),
        legend.title=element_blank(),
        legend.position = 'none',
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.text.x = element_text(size=18,face="plain",color="black"),#,angle = 45,vjust = 1,hjust = 1
        legend.text = element_text(size=18,face="plain",color="black"))

# p
finaldata$ID <- paste(finaldata$cancerName,finaldata$TNType,finaldata$OmicName,sep='_')
plotdata <- data.frame(table(finaldata$ID))
plotdata[,c('cancerName','TNType','OmicName')] <- str_split_fixed(plotdata$Var1,'_',3)
plotdata$TNType <- factor(plotdata$TNType,levels =c('Normal','Tumor'))
plotdata$OmicName <- factor(plotdata$OmicName,levels =c('Pro','Phos','RNA','RRA'))

p1 <- ggplot(plotdata,aes(x=OmicName,y=Freq,fill=TNType))+
  geom_bar(stat='summary',fun='median',position = position_dodge(),width = 0.8)+
  # stat_compare_means(label = "p.format",method="wilcox.test")+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge( .9))+
  scale_fill_manual(values=c('#455a82','#7c3031'))+
  labs(x="",y="Number of Exp-confirmed\nCorrelated Pairs")+
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        legend.box = element_blank(),
        legend.background = element_blank(),
        legend.title=element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.25,0.85),
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.text.x = element_text(size=18,face="plain",color="black"),#,angle = 45,vjust = 1,hjust = 1
        legend.text = element_text(size=18,face="plain",color="black"))

# p1
ggsave('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure2A_exp_confirmed_Number_spearman.pdf',p1|p,width=10,height = 4)

############Figure2B write GeneList for 5 tools###############
# ToolPreList <- list()
# GPShighDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_High.rds')
# GPShighDF$ID <- paste(GPShighDF$kinaseName,str_split_fixed(GPShighDF$Pro.Site,'#',2)[,2],sep='#')
# print('GPShighDF')
# ToolPreList[['GPS_high']] <- GPShighDF$ID
#
# GPSmedianDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_Medium.rds')
# GPSmedianDF$ID <- paste(GPSmedianDF$kinaseName,str_split_fixed(GPSmedianDF$Pro.Site,'#',2)[,2],sep='#')
# print('GPSmedianDF')
# ToolPreList[['GPS_Medium']] <- GPSmedianDF$ID
#
# exprimentDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/Expriment.clean.rds')
# exprimentDF$ID <- paste(exprimentDF$Kinase_names,str_split_fixed(exprimentDF$Pro.site,'#',2)[,2],sep='#')
# ToolPreList[['Expriment']] <- exprimentDF$ID
#
# ScansiteDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/Scansite.clean.rds')
# ScansiteDF$ID <- paste(ScansiteDF$Kinase_names,str_split_fixed(ScansiteDF$Pro.site,'#',2)[,2],sep='#')
# ToolPreList[['Scansite']] <- ScansiteDF$ID
#
# MusiteDeepDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/MusiteDeep.clean.rds')
# MusiteDeepDF$ID <- paste(MusiteDeepDF$Kinase_names,str_split_fixed(MusiteDeepDF$Pro.site,'#',2)[,2],sep='#')
# print('MusiteDeepDF')
# ToolPreList[['MusiteDeep']] <- MusiteDeepDF$ID
#
# NetworKINDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/NetworKIN.clean.rds')
# NetworKINDF$ID <- paste(NetworKINDF$Kinase_names,str_split_fixed(NetworKINDF$Pro.site,'#',2)[,2],sep='#')
# print('NetworKINDF')
# ToolPreList[['NetworKIN']] <- NetworKINDF$ID
#
# save(ToolPreList, file = "/home/wuyj2/eKPI/Prediction_Tool/PredictToolList_GeneName.RData")


setwd('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Spearman_RRAV2/')
filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','TNType')] <- str_split_fixed(str_replace_all(cancerFrame$fileName,'.tsv',''),'_',2)

kinasedf <- read.csv('/home/wuyj2/eKPI/CancerOmics/kinase_basic(Sheet1).csv',header = T)
kinaseNames <- union(kinasedf$Offical_gene_symbol,kinasedf$KinBase_name)[-1]
Sample_THRD <- 6
Pvalue_THRD <- 0.05
Rvalue_THRD <- 0.2
RRA_THRD <- 1
print('GSEA and Prediction Tool for each project')
write_dir <- '/hwdata/home/wuyj2/eKPI/ArticleCancerStat/fGSEA_rst/'
load("/home/wuyj2/eKPI/Prediction_Tool/PredictToolList_GeneName.RData")

finaldata <- data.frame()
for(i in 1:nrow(cancerFrame)){
  filedata <- fread(cancerFrame[i,'fileName'],header = T,sep='\t',data.table = F)
  filedata <- filedata[order(filedata$Score,decreasing = F),]
  rownames(filedata) <- filedata$ID
  filedata2 <- filedata#[filedata$Score<RRA_THRD,]
  # filedata2 <- filedata[1:round(nrow(filedata)/5,0),]
  if(ncol(filedata)==11){
    omicNames <- c('Phos','Pro','RNA')
  }
  if(ncol(filedata)==9){
    omicNames <- c('Phos','Pro')
  }
  for(omicName in omicNames){
    write_file <- paste(write_dir,cancerFrame[i,'CancerName'],'_',cancerFrame[i,'TNType'],'_',omicName,'.tsv',sep='')
    print(write_file)
    if(!file.exists(write_file)){
      genelist <- filedata2[!is.na(filedata2[,paste(omicName,'Rho',sep = '.')]),c('ID',paste(omicName,'Rho',sep = '.'))]
      genelist <- genelist[order(genelist[,2],decreasing = T),]
      
      genelistvalue <- genelist[,2]
      names(genelistvalue) <- genelist$ID
      gsea.re2 <- fgsea(pathways = ToolPreList,#基因集列表
                        stats = genelistvalue,#排序后的基因level,这里是logFC
                        # nperm=1000,#置换检验的次数
                        minSize=1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
                        maxSize=100000000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
                        #                 nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
                        #                 gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
                        #                 BPPARAM  = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
      )
      gesa.re <- data.frame(pathway = gsea.re2[[1]],
                            pval = gsea.re2[[2]],
                            gesa.re2 = gsea.re2[[3]],
                            ES = gsea.re2[[4]],
                            NES = gsea.re2[[5]],
                            nMoreExtreme = gsea.re2[[6]],
                            size =  gsea.re2[[7]],
                            leadingEdge = sapply(gsea.re2[[8]],paste,collapse=';'))
      write.table(gesa.re,write_file,quote = F,sep='\t')
    }
  }
}

setwd('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Spearman_RRAV2/')
filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','TNType')] <- str_split_fixed(str_replace_all(cancerFrame$fileName,'.tsv',''),'_',2)

kinasedf <- read.csv('/home/wuyj2/eKPI/CancerOmics/kinase_basic(Sheet1).csv',header = T)
kinaseNames <- union(kinasedf$Offical_gene_symbol,kinasedf$KinBase_name)[-1]
Sample_THRD <- 6
Pvalue_THRD <- 0.05
Rvalue_THRD <- 0.2
RRA_THRD <- 0.05
print('GSEA and Prediction Tool for each project')
write_dir <- '/hwdata/home/wuyj2/eKPI/ArticleCancerStat/fGSEA_rst/'
load("/home/wuyj2/eKPI/Prediction_Tool/PredictToolList_GeneName.RData")

finaldata <- data.frame()
for(i in 1:nrow(cancerFrame)){
  filedata <- fread(cancerFrame[i,'fileName'],header = T,sep='\t',data.table = F)
  filedata <- filedata[order(filedata$Score,decreasing = F),]
  rownames(filedata) <- filedata$ID
  filedata2 <- filedata#[filedata$Score<RRA_THRD,]
  # filedata2 <- filedata[1:round(nrow(filedata)/5,0),]
  
  omicNames <- c('Score')
  
  for(omicName in omicNames){
    write_file <- paste(write_dir,cancerFrame[i,'CancerName'],'_',cancerFrame[i,'TNType'],'_','RRA.tsv',sep='')
    print(write_file)
    if(!file.exists(write_file)){
      genelist <- filedata2[!is.na(filedata2[,omicName]),c('ID',omicName)]
      genelist <- genelist[genelist$Score<=RRA_THRD,]
      genelist$log10RRA <- -log10(genelist$Score)
      genelist <- genelist[order(genelist[,'log10RRA'],decreasing = T),]
      
      genelistvalue <- genelist[,'log10RRA']
      names(genelistvalue) <- genelist$ID
      gsea.re2 <- fgsea(pathways = ToolPreList,#基因集列表
                        stats = genelistvalue,#排序后的基因level,这里是logFC
                        # nperm=1000,#置换检验的次数
                        minSize=1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
                        maxSize=100000000,#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
                        #                 nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
                        #                 gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
                        #                 BPPARAM  = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用
      )
      gesa.re <- data.frame(pathway = gsea.re2[[1]],
                            pval = gsea.re2[[2]],
                            gesa.re2 = gsea.re2[[3]],
                            ES = gsea.re2[[4]],
                            NES = gsea.re2[[5]],
                            nMoreExtreme = gsea.re2[[6]],
                            size =  gsea.re2[[7]],
                            leadingEdge = sapply(gsea.re2[[8]],paste,collapse=';'))
      write.table(gesa.re,write_file,quote = F,sep='\t')
    }
  }
}

############Figure2B GSEA Dotplot############
finaldata <- data.frame()
setwd('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/fGSEA_rst/')
fileNames <- list.files('.')
for(filename in fileNames){
  tmpdata <- fread(filename,data.table = F,sep='\t')
  tmpdata$fileName <- str_replace(filename,'\\.tsv','')
  finaldata <- rbind(finaldata,tmpdata[,c(2:8,10)])
}
finaldata[is.na(finaldata$pval),'pval'] <- 1
finaldata[is.na(finaldata$gesa.re2),'gesa.re2'] <- 1
finaldata$FDR <- p.adjust(finaldata$pval,method = 'BH')
# colnames(finaldata) <- c("Project","Omics","TNType","Proof.Cor","Proof.UnCor",
#                          "NonProof.Cor","NonProof.Uncor","Oddsr" ,"pvalue", "OddInter1","OddInter2","ExternalProof")

write.table(finaldata,'/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/Figure2B_Predict_Spearman_fgsea_result.tsv',sep = '\t',quote = F,row.names = F)

GSEA_rst <- finaldata[!is.na(finaldata$NES),]
GSEA_rst[,c('CancerName','TNtype','Omics')] <- str_split_fixed(GSEA_rst$fileName,'_',3)
GSEA_rst$xName <- ifelse(GSEA_rst$TNtype=='Tumor',
                         paste(GSEA_rst$CancerName,'T',sep ='_'),
                         paste(GSEA_rst$CancerName,'N',sep ='_'))
GSEA_rst$Omics <- factor(GSEA_rst$Omics,levels = c('RRA','RNA','Phos','Pro'))
GSEA_rst$pathway <- factor(GSEA_rst$pathway,levels = c('Expriment',
                                                       'Scansite',
                                                       'NetworKIN',
                                                       'GPS_high',
                                                       'GPS_Medium',
                                                       'MusiteDeep'))
GSEA_rst$pvalueLabel <- ifelse(GSEA_rst$FDR<0.001,0.001,
                               ifelse(GSEA_rst$FDR<0.01,0.01,
                                      ifelse(GSEA_rst$FDR<0.05,0.05,
                                             ifelse(GSEA_rst$FDR<0.5,0.5,0.7))))
GSEA_rst$logpvalue <- -log10(GSEA_rst$pvalueLabel)
GSEA_rst$NESLabel <- ifelse(GSEA_rst$NES>=0.25&GSEA_rst$FDR<0.05,'NES>=0.25&FDR<0.05','others')
GSEA_rst$NESLabel <- factor(GSEA_rst$NESLabel,levels = c('NES>=0.25&FDR<0.05','others'))
GSEA_rst$xName <- factor(GSEA_rst$xName,levels=c('BRCA.Cell_T',
                                                 'CCRCC_T',
                                                 'COC_T',
                                                 'DTGC_T',
                                                 'ECA_T',
                                                 'ENCA_T',
                                                 'EOGC_T',
                                                 'ESCC_T',
                                                 'ESHCC_T',
                                                 'GBM_T',
                                                 'HBVHCC_T',
                                                 'HGOV_T',
                                                 'HNSCC_T',
                                                 'ICCA_T',
                                                 'LSCC_T',
                                                 'LUAD.CN_T',
                                                 'LUAD.US_T',
                                                 'mCRC_T',
                                                 'MEBL_T',
                                                 'NSLC_T',
                                                 'PBC_T',
                                                 'PDAC_T',
                                                 'TNBC_T',
                                                 'CCRCC_N',
                                                 'COC_N',
                                                 'DTGC_N',
                                                 'ECA_N',
                                                 'ENCA_N',
                                                 'EOGC_N',
                                                 'ESCC_N',
                                                 'ESHCC_N',
                                                 'GBM_N',
                                                 'HBVHCC_N',
                                                 'HNSCC_N',
                                                 'LSCC_N',
                                                 'LUAD.CN_N',
                                                 'LUAD.US_N',
                                                 'NSLC_N',
                                                 'PDAC_N'))

# GSEA_rst1 <- GSEA_rst[GSEA_rst$pathway=='Expriment',]
g1 <- ggplot(GSEA_rst,aes(y = xName, x = Omics)) +
  coord_flip()+
  geom_point(mapping = aes(size=logpvalue,fill = NES,colour = NESLabel),shape=21) +
  scale_size_continuous(name='',breaks=c(0.3,1.3,2,3),range = c(2,7),labels = c('<0.5','<0.05','<0.01','<0.001'))+
  facet_wrap(~pathway,nrow=6,strip.position='right')+
  scale_color_manual(values = c('black','grey'))+
  scale_fill_gradient2(low="#01508b",high="#8f323a",mid='white',midpoint = 0) + #自定义配色
  # scale_y_continuous(position = "left",
  #                    breaks = 1:nrow(df),
  #                    labels = Hmisc::capitalize(rev(df$geneSet))) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(title = "FDR",order=1),
         fill = guide_colorbar(title = "Normalized ES",order=2),
         color = guide_legend(title = "",order=3)) +
  # ggtitle('Expriment') +
  theme_bw(base_size = 15) +
  theme(panel.grid =element_blank(),
        panel.grid.major = element_line(color = "gray"), #, linetype = "dashed"
        legend.position = 'bottom',
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.text = element_text(size = 15, family = 'sans'),
        plot.title = element_text(size = 15, family = 'sans',hjust=0.5),
        axis.text.y = element_text(size = 15, family = 'sans'),
        axis.text.x = element_text(angle = 45,vjust=1,hjust = 1,size = 15)) #去除网格线
g1
ggsave(paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/','Figure2B_','fgsea_NES_DotPoints.pdf',sep = ''),g1,width = 12.5,height = 9)


############FigureS2AB GSEA plot####################
library(msigdbr)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
RRA_THRD <- 0.05
load("/home/wuyj2/eKPI/Prediction_Tool/PredictToolList_GeneName.RData")
myDB <- data.frame()
for(ListName in names(ToolPreList)){
  KPSNames <- ToolPreList[[ListName]]
  myDB <- rbind(myDB,data.frame(gs_name = ListName,
                                gene_symbol = KPSNames))
}
setwd('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Spearman_RRAV2/')
filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','TNType')] <- str_split_fixed(str_replace_all(cancerFrame$fileName,'.tsv',''),'_',2)

finalRst <- data.frame()
ROWINdexs <- c(1,20,23,24,26,31,34)
for(i in ROWINdexs){
  rowIndex <- i
  filedata <- fread(cancerFrame[rowIndex,'fileName'],header = T,sep='\t',data.table = F)
  filedata <- filedata[order(filedata$Score,decreasing = F),]
  rownames(filedata) <- filedata$ID
  filedata2 <- filedata[filedata$Score<RRA_THRD,]
  
  #Annotation Prediction FALSE/True
  for(predictName in names(ToolPreList)){
    filedata2[,predictName] <- FALSE
    KPSNames <- ToolPreList[[predictName]]
    filedata2[which(filedata2$ID %in% KPSNames),predictName] <- TRUE
  }
  
  plotdata <- data.frame()
  for(omicName in c('RNA','Pro','Phos')){
    tmpName <- paste(omicName,'Rho',sep='.')
    tmpdata <- filedata2[!is.na(filedata2[,tmpName]),]
    
    for(predictName in names(ToolPreList)){
      # print(predictName)
      tmpdata1 <- tmpdata[tmpdata[,predictName]==FALSE,]
      plotdata <- rbind(plotdata,data.frame(ID  = tmpdata1$ID,
                                            omicNames = omicName,
                                            ExternalProof = paste(predictName,'All',sep = '_'),
                                            Rho = tmpdata1[,tmpName]))
      
      tmpdata2 <- tmpdata[tmpdata[,predictName]==TRUE,]
      print(paste(omicName,predictName,nrow(tmpdata2),sep='_'))
      plotdata <- rbind(plotdata,data.frame(ID  = tmpdata2$ID,
                                            omicNames = omicName,
                                            ExternalProof = predictName,
                                            Rho = tmpdata2[,tmpName]))
    }
  }
  plotdata2 <- plotdata[!str_detect(plotdata$ExternalProof,'_All'),]
  
  
  for(omicName in c('RNA','Pro','Phos')){
    for(predictName in names(ToolPreList)){
      tmpdata <- plotdata2[plotdata2$omicNames==omicName&plotdata2$ExternalProof==predictName,]
      finalRst <- rbind(finalRst,data.frame(OmicName = omicName,
                                            Project = paste(cancerFrame[rowIndex,'CancerName'],cancerFrame[rowIndex,'TNType'],sep='_'),
                                            ExternalProof = predictName,
                                            MedianValue = round(median(tmpdata$Rho),3),
                                            upquantile = quantile(tmpdata$Rho,0.75)))
    }
  }
  
  tt <- ttheme_minimal(base_size = 14,
                       core=list(
                         #bg_params = list(fill = NA, col=NA),
                         fg_params=list(col=c('#F1B543','#1963B3','#025939','#A5405E','#7C1A97','#019092')))
  )
  gglist <- list()
  for(omicName in c('Pro','Phos','RNA')){
    tmpdata <- finalRst[finalRst$OmicName==omicName&finalRst$Project==paste(cancerFrame[rowIndex,'CancerName'],cancerFrame[rowIndex,'TNType'],sep='_'),]
    tmpdata$ExternalProof <- factor(tmpdata$ExternalProof,levels = c('Expriment','Scansite','NetworKIN','GPS_high','GPS_Medium','MusiteDeep'))
    tmpdata <- tmpdata[order(tmpdata$ExternalProof),]
    tp <- tableGrob(tmpdata[,'MedianValue'],rows = NULL,theme = tt)
    tp$heights <- unit(rep(0.5,nrow(tp)),"cm") # cell height
    
    # 修改表格每个格子的宽度和高度
    #tp$widths <- unit(rep(1.2,ncol(tp)), "cm")
    
    tp$heights <- unit(rep(0.5,nrow(tp)),"cm") # cell height
    plotdata3 <- plotdata2[plotdata2$omicNames==omicName,]
    plotdata3$ExternalProof <- factor(plotdata3$ExternalProof,levels = c('Expriment','Scansite','NetworKIN','GPS_high','GPS_Medium','MusiteDeep'))
    g1 <- ggplot(plotdata3, aes(x=Rho,Group=ExternalProof,fill = ExternalProof))+
      geom_density(alpha=0.4)+
      # geom_rug()+
      facet_grid(~omicNames)+
      # scale_color_manual(values=c('#0C4E9B','#F98F34','#C72228'))+
      scale_fill_manual(values=c('#F1B543','#1963B3','#025939','#A5405E','#7C1A97','#019092'))+
      theme_bw()+
      # geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
      # geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
      # geom_vline(xintercept = quantile(phosdf_1$Rho,0.75),linetype=2) + # 添加base mean的水平线
      scale_x_continuous(name = "Spearman' Rho")+
      scale_y_continuous(name = "Density")+
      # scale_fill_manual(values=c('#B54764','#9D9EA3'),
      #                   breaks = c(TRUE,FALSE),
      #                   labels = c('Exp-confirmed','Non Exp-confirmed'))+
      annotation_custom(tp,
                        xmin = 0.7,
                        xmax = 1
      )+
      theme_classic(base_size = 18)+
      theme(panel.grid=element_blank(),
            legend.box = element_blank(),
            # legend.position= "top",
            panel.border = element_rect(fill=NA,linewidth = 1),
            # text=element_text(size=14,face="plain",color="black"),
            # axis.title=element_text(size=16,face="plain",color="black"),
            # axis.text = element_text(size=14,face="plain",color="black"),
            legend.title = element_blank(),
            # legend.text = element_text(size=14,face="plain",color="black"),
            legend.background = element_blank())
    # g1
    g1
    gglist[[omicName]] <- g1
  }
  
  p <- wrap_plots(gglist)
  p
  ggsave(paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/','Figure2C_',cancerFrame[rowIndex,'CancerName'],'_',cancerFrame[rowIndex,'TNType'],'_Omics_Density.pdf',sep = ''),p,width = 10,height = 3.5)
}

# g1 <- ggplot(plotdata2, aes(x=Rho,Group=ExternalProof,fill = ExternalProof))+
#   geom_density(alpha=0.8)+
#   # geom_rug()+
#   facet_grid(~omicNames)+
#   # scale_color_manual(values=c('#0C4E9B','#F98F34','#C72228'))+
#   scale_fill_manual(breaks = c('Expriment','Expriment_All'),
#                     values=c('#B54764','#9D9EA3'),
#                     labels = c('Expriment-confirmed','All'))+
#   theme_bw()+
#   # geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
#   # geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
#   # geom_vline(xintercept = quantile(phosdf_1$Rho,0.75),linetype=2) + # 添加base mean的水平线
#   scale_x_continuous(name = "Spearman' Rho")+
#   scale_y_continuous(name = "Density")+
#   # scale_fill_manual(values=c('#B54764','#9D9EA3'),
#   #                   breaks = c(TRUE,FALSE),
#   #                   labels = c('Exp-confirmed','Non Exp-confirmed'))+
#   theme_classic(base_size = 18)+
#   theme(panel.grid=element_blank(),
#         legend.box = element_blank(),
#         # legend.position= "inside",
#         # legend.position.inside = c(0.75,0.9),
#         panel.border = element_rect(fill=NA,linewidth = 1),
#         # text=element_text(size=14,face="plain",color="black"),
#         # axis.title=element_text(size=16,face="plain",color="black"),
#         # axis.text = element_text(size=14,face="plain",color="black"),
#         legend.title = element_blank(),
#         # legend.text = element_text(size=14,face="plain",color="black"),
#         legend.background = element_blank())
# # g1
# g1

rowIndex <- 31
filedata <- fread(cancerFrame[rowIndex,'fileName'],header = T,sep='\t',data.table = F)
filedata <- filedata[order(filedata$Score,decreasing = F),]
rownames(filedata) <- filedata$ID

filedata2 <- filedata[filedata$Score<RRA_THRD,]
filedata2$log10RRA <- -log10(filedata2$Score)
genelist <- filedata2[,c('ID','log10RRA')]
genelist <- genelist[order(genelist[,2],decreasing = T),]

genelistvalue <- genelist[,2]
names(genelistvalue) <- genelist$ID


gseahallmark_KEGG2 <- GSEA(genelistvalue, TERM2GENE = myDB, pvalueCutoff = 1,minGSSize = 1,eps = 0,
                           maxGSSize = 500000000,pAdjustMethod='BH',by='fgsea')
gsearst <- gseahallmark_KEGG2@result
gsearst <- gsearst[c('Expriment','Scansite','NetworKIN','GPS_high','GPS_Medium','MusiteDeep'),]
gsearst$color <-  c('#F1B543','#1963B3','#025939','#A5405E','#7C1A97','#019092')
# gsearst <- gsearst[order(gsearst$NES,decreasing = T),]
gsearst <- gsearst[gsearst$p.adjust<0.05,]
geneSetID <- gsearst$ID#c('Expriment','GPS_high','GPS_Medium','MusiteDeep')#'Scansite','NetworKIN',,
pd <- gseahallmark_KEGG2[geneSetID, c( "NES", "p.adjust")]# 提取NES，P值等信息
for (i in seq_len(ncol(pd))) {pd[, i] <- format(pd[, i], digits = 3)}

# 通过修改table的主题来修改表格细节
tt <- ttheme_minimal(base_size = 14,
                     core=list(
                       #bg_params = list(fill = NA, col=NA),
                       fg_params=list(col=gsearst$color))
)
tp <- tableGrob(pd,rows = NULL,theme = tt)

# 修改表格每个格子的宽度和高度
#tp$widths <- unit(rep(1.2,ncol(tp)), "cm")
tp$heights <- unit(rep(0.5,nrow(tp)),"cm") # cell height


#L.US_Tumor
pdf(paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure2C',cancerFrame[rowIndex,'CancerName'],cancerFrame[rowIndex,'TNType'],omicName,'GSEAplot_RRA.pdf',sep='_'),width =5,height = 5)
p <- gseaplot2(gseahallmark_KEGG2,
               geneSetID = geneSetID,#展示的通路序号
               # title=gse_rst$Description[c(5,8,2,23)],
               subplots = 1,
               base_size = 20,
               ES_geom ='line',
               rel_heights = c(1.5, 0.5, 1),
               color = c('#F1B543','#1963B3','#025939','#A5405E','#7C1A97','#019092')#'#A5405E','#025939',
               # pvalue_table = T
)
p <- p+
  ggtitle(paste(cancerFrame[rowIndex,'CancerName'],cancerFrame[rowIndex,'TNType'],sep='_'))+
  annotation_custom(tp,
                    xmin = 64000,
                    xmax = 75000,
                    ymin = 0.30,
                    ymax = 0.6
  )+
  guides(color=guide_legend(nrow=2))+
  theme_bw(base_size = 18)+
  theme(plot.title = element_text(size = 18),
        legend.position = "top",
        legend.title = element_blank(),
        legend.direction = "vertical"
  )

p
dev.off()

############TableS3 and Prediction Tool for each project############
setwd('/hwdata/home/wuyj2/eKPI/spearman/')
sample_Thrd <- 6

filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','PMID','TNType','Omics')] <- str_split_fixed(cancerFrame$fileName,'_',5)
cancerFrame$Omics <- str_replace_all(cancerFrame$Omics,'Kinase.','')
cancerFrame$ID <- paste(cancerFrame$CancerName,cancerFrame$TNType,sep='_')

kinasedf <- read.csv('/hwdata/home/wuyj2/eKPI/CancerOmics/kinase_basic(Sheet1).csv',header = T)
kinaseNames <- union(kinasedf$Offical_gene_symbol,kinasedf$KinBase_name)[-1]

Sample_THRD <- 6
Pvalue_THRD <- 0.05
Rvalue_THRD <- 0.2

print('Union Pro_Phos_RNA and Prediction Tool for each project V2')

exprimentDF <- readRDS('/hwdata/home/wuyj2/eKPI/Prediction_Tool/Expriment.clean.rds')
exprimentDF$ID <- paste(exprimentDF$Kinase_names,str_split_fixed(exprimentDF$Pro.site,'#',2)[,2],sep='#')

ScansiteDF <- readRDS('/hwdata/home/wuyj2/eKPI/Prediction_Tool/Scansite.clean.rds')
ScansiteDF$ID <- paste(ScansiteDF$Kinase_names,str_split_fixed(ScansiteDF$Pro.site,'#',2)[,2],sep='#')

MusiteDeepDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/MusiteDeep.clean.rds')
MusiteDeepDF$ID <- paste(MusiteDeepDF$Kinase_names,str_split_fixed(MusiteDeepDF$Pro.site,'#',2)[,2],sep='#')
print('MusiteDeepDF')

NetworKINDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/NetworKIN.clean.rds')
NetworKINDF$ID <- paste(NetworKINDF$Kinase_names,str_split_fixed(NetworKINDF$Pro.site,'#',2)[,2],sep='#')
print('NetworKINDF')

GPShighDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_High.rds')
GPShighDF$ID <- paste(GPShighDF$kinaseName,str_split_fixed(GPShighDF$Pro.Site,'#',2)[,2],sep='#')
print('GPShighDF')

GPSmedianDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_Medium.rds')
GPSmedianDF$ID <- paste(GPSmedianDF$kinaseName,str_split_fixed(GPSmedianDF$Pro.Site,'#',2)[,2],sep='#')
print('GPSmedianDF')

exprimentUnion <- c()
ScansiteUnion <- c()
MusiteDeepUnion <- c()
NetworKINUnion <- c()
GPShighUnion <- c()
GPSmedianUnion <- c()

ProjectNames <- unique(cancerFrame$ID)
for(PjName in ProjectNames){
  print(PjName)
  myUnionDF <- data.frame()
  fileDataFrame <- cancerFrame[cancerFrame$ID==PjName,]
  # write_file <- paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Spearman_RRAV2/',PjName,'.tsv',sep='')
  if(nrow(fileDataFrame)>0){
    if('Phos' %in% fileDataFrame$Omics){
      filedata <- fread(paste(fileDataFrame[fileDataFrame$Omics=='Phos','fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[!is.na(filedata$Rho),]
      filedata <- filedata[filedata$Rho>0&filedata$sample.num>=Sample_THRD,]
      
      filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
      filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$FDR <- p.adjust(filedata$Rho.pvalue,method='BH')
      
      max_sample <- max(c(217,max(filedata$sample.num)))
      k1 <- max_sample^3-max_sample
      filedata$p.adjust <- ((k1-(filedata$sample.num^3-filedata$sample.num))*filedata$Rho.pvalue+(filedata$sample.num^3-filedata$sample.num)*filedata$FDR)/k1
      
      
      filedata <- filedata[filedata$sample.num>=Sample_THRD&
                             filedata$Rho>=Rvalue_THRD&
                             filedata$p.adjust<=Pvalue_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair
      if(nrow(filedata)>0){
        filedata$ID  <- paste(filedata$KinaseName,str_split_fixed(filedata$GeneSite,'#',2)[,2],sep='#')
        
        filedata$expriment <- FALSE
        filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
        filedata$expriment <- factor(filedata$expriment,levels = c(TRUE,FALSE))
        exprimentUnion <- unique(c(exprimentUnion,filedata[filedata$expriment==TRUE,'ID']))
        print(paste('experiment',length(exprimentUnion),sep = ' '))
        
        filedata$Scansite <- FALSE
        filedata[which(filedata$ID %in% ScansiteDF$ID),]$Scansite <- TRUE
        filedata$Scansite <- factor(filedata$Scansite,levels = c(TRUE,FALSE))
        ScansiteUnion <- unique(c(ScansiteUnion,filedata[filedata$Scansite==TRUE,'ID']))
        print(paste('Scansite',length(ScansiteUnion),sep = ' '))
        
        filedata$MusiteDeep <- FALSE
        filedata[which(filedata$ID %in% MusiteDeepDF$ID),]$MusiteDeep <- TRUE
        filedata$MusiteDeep <- factor(filedata$MusiteDeep,levels = c(TRUE,FALSE))
        MusiteDeepUnion <- unique(c(MusiteDeepUnion,filedata[filedata$MusiteDeep==TRUE,'ID']))
        print(paste('MusiteDeep',length(MusiteDeepUnion),sep = ' '))

        filedata$NetworKIN <- FALSE
        filedata[which(filedata$ID %in% NetworKINDF$ID),]$NetworKIN <- TRUE
        filedata$NetworKIN <- factor(filedata$NetworKIN,levels = c(TRUE,FALSE))
        NetworKINUnion <- unique(c(NetworKINUnion,filedata[filedata$NetworKIN==TRUE,'ID']))
        print(paste('NetworKIN',length(NetworKINUnion),sep = ' '))

        filedata$GPShigh <- FALSE
        filedata[which(filedata$ID %in% GPShighDF$ID),]$GPShigh <- TRUE
        filedata$GPShigh <- factor(filedata$GPShigh,levels = c(TRUE,FALSE))
        GPShighUnion <- unique(c(GPShighUnion,filedata[filedata$GPShigh==TRUE,'ID']))
        print(paste('GPShigh',length(GPShighUnion),sep = ' '))

        filedata$GPSmedian <- FALSE
        filedata[which(filedata$ID %in% GPSmedianDF$ID),]$GPSmedian <- TRUE
        filedata$GPSmedian <- factor(filedata$GPSmedian,levels = c(TRUE,FALSE))
        GPSmedianUnion <- unique(c(GPSmedianUnion,filedata[filedata$GPSmedian==TRUE,'ID']))
        print(paste('GPSmedian',length(GPSmedianUnion),sep = ' '))

      }
      print('Phos')
    }
    if('Pro' %in% fileDataFrame$Omics){
      filedata <- fread(paste(fileDataFrame[fileDataFrame$Omics=='Pro','fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[!is.na(filedata$Rho),]
      filedata <- filedata[filedata$Rho>0&filedata$sample.num>=Sample_THRD,]
      
      filedata <- filedata[which(filedata$RNA.names %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$RNA.names!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$FDR <- p.adjust(filedata$Rho.pvalue,method='BH')
      
      max_sample <- max(c(217,max(filedata$sample.num)))
      k1 <- max_sample^3-max_sample
      filedata$p.adjust <- ((k1-(filedata$sample.num^3-filedata$sample.num))*filedata$Rho.pvalue+(filedata$sample.num^3-filedata$sample.num)*filedata$FDR)/k1
      
      
      filedata <- filedata[filedata$sample.num>=Sample_THRD&
                             filedata$Rho>=Rvalue_THRD&
                             filedata$p.adjust<=Pvalue_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair
      if(nrow(filedata)>0){
        filedata$ID  <- paste(filedata$RNA.names,str_split_fixed(filedata$GeneSite,'#',2)[,2],sep='#')
        filedata$expriment <- FALSE
        filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
        filedata$expriment <- factor(filedata$expriment,levels = c(TRUE,FALSE))
        exprimentUnion <- unique(c(exprimentUnion,filedata[filedata$expriment==TRUE,'ID']))
        print(paste('experiment',length(exprimentUnion),sep = ' '))
        
        filedata$Scansite <- FALSE
        filedata[which(filedata$ID %in% ScansiteDF$ID),]$Scansite <- TRUE
        filedata$Scansite <- factor(filedata$Scansite,levels = c(TRUE,FALSE))
        ScansiteUnion <- unique(c(ScansiteUnion,filedata[filedata$Scansite==TRUE,'ID']))
        print(paste('Scansite',length(ScansiteUnion),sep = ' '))
        
        filedata$MusiteDeep <- FALSE
        filedata[which(filedata$ID %in% MusiteDeepDF$ID),]$MusiteDeep <- TRUE
        filedata$MusiteDeep <- factor(filedata$MusiteDeep,levels = c(TRUE,FALSE))
        MusiteDeepUnion <- unique(c(MusiteDeepUnion,filedata[filedata$MusiteDeep==TRUE,'ID']))
        print(paste('MusiteDeep',length(MusiteDeepUnion),sep = ' '))
        
        filedata$NetworKIN <- FALSE
        filedata[which(filedata$ID %in% NetworKINDF$ID),]$NetworKIN <- TRUE
        filedata$NetworKIN <- factor(filedata$NetworKIN,levels = c(TRUE,FALSE))
        NetworKINUnion <- unique(c(NetworKINUnion,filedata[filedata$NetworKIN==TRUE,'ID']))
        print(paste('NetworKIN',length(NetworKINUnion),sep = ' '))
        
        filedata$GPShigh <- FALSE
        filedata[which(filedata$ID %in% GPShighDF$ID),]$GPShigh <- TRUE
        filedata$GPShigh <- factor(filedata$GPShigh,levels = c(TRUE,FALSE))
        GPShighUnion <- unique(c(GPShighUnion,filedata[filedata$GPShigh==TRUE,'ID']))
        print(paste('GPShigh',length(GPShighUnion),sep = ' '))
        
        filedata$GPSmedian <- FALSE
        filedata[which(filedata$ID %in% GPSmedianDF$ID),]$GPSmedian <- TRUE
        filedata$GPSmedian <- factor(filedata$GPSmedian,levels = c(TRUE,FALSE))
        GPSmedianUnion <- unique(c(GPSmedianUnion,filedata[filedata$GPSmedian==TRUE,'ID']))
        print(paste('GPSmedian',length(GPSmedianUnion),sep = ' '))      }
      print('Pro')
    }
    if('RNA' %in% fileDataFrame$Omics){
      filedata <- fread(paste(fileDataFrame[fileDataFrame$Omics=='RNA','fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[!is.na(filedata$Rho),]
      filedata <- filedata[filedata$Rho>0&filedata$sample.num>=Sample_THRD,]
      
      filedata <- filedata[which(filedata$RNA.names %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$RNA.names!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$FDR <- p.adjust(filedata$Rho.pvalue,method='BH')
      
      max_sample <- max(c(217,max(filedata$sample.num)))
      k1 <- max_sample^3-max_sample
      filedata$p.adjust <- ((k1-(filedata$sample.num^3-filedata$sample.num))*filedata$Rho.pvalue+(filedata$sample.num^3-filedata$sample.num)*filedata$FDR)/k1
      
      
      filedata <- filedata[filedata$sample.num>=Sample_THRD&
                             filedata$Rho>=Rvalue_THRD&
                             filedata$p.adjust<=Pvalue_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair
      if(nrow(filedata)>0){
        filedata$ID  <- paste(filedata$RNA.names,str_split_fixed(filedata$GeneSite,'#',2)[,2],sep='#')
        filedata$expriment <- FALSE
        filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
        filedata$expriment <- factor(filedata$expriment,levels = c(TRUE,FALSE))
        exprimentUnion <- unique(c(exprimentUnion,filedata[filedata$expriment==TRUE,'ID']))
        print(paste('experiment',length(exprimentUnion),sep = ' '))
        
        filedata$Scansite <- FALSE
        filedata[which(filedata$ID %in% ScansiteDF$ID),]$Scansite <- TRUE
        filedata$Scansite <- factor(filedata$Scansite,levels = c(TRUE,FALSE))
        ScansiteUnion <- unique(c(ScansiteUnion,filedata[filedata$Scansite==TRUE,'ID']))
        print(paste('Scansite',length(ScansiteUnion),sep = ' '))
        
        filedata$MusiteDeep <- FALSE
        filedata[which(filedata$ID %in% MusiteDeepDF$ID),]$MusiteDeep <- TRUE
        filedata$MusiteDeep <- factor(filedata$MusiteDeep,levels = c(TRUE,FALSE))
        MusiteDeepUnion <- unique(c(MusiteDeepUnion,filedata[filedata$MusiteDeep==TRUE,'ID']))
        print(paste('MusiteDeep',length(MusiteDeepUnion),sep = ' '))
        
        filedata$NetworKIN <- FALSE
        filedata[which(filedata$ID %in% NetworKINDF$ID),]$NetworKIN <- TRUE
        filedata$NetworKIN <- factor(filedata$NetworKIN,levels = c(TRUE,FALSE))
        NetworKINUnion <- unique(c(NetworKINUnion,filedata[filedata$NetworKIN==TRUE,'ID']))
        print(paste('NetworKIN',length(NetworKINUnion),sep = ' '))
        
        filedata$GPShigh <- FALSE
        filedata[which(filedata$ID %in% GPShighDF$ID),]$GPShigh <- TRUE
        filedata$GPShigh <- factor(filedata$GPShigh,levels = c(TRUE,FALSE))
        GPShighUnion <- unique(c(GPShighUnion,filedata[filedata$GPShigh==TRUE,'ID']))
        print(paste('GPShigh',length(GPShighUnion),sep = ' '))
        
        filedata$GPSmedian <- FALSE
        filedata[which(filedata$ID %in% GPSmedianDF$ID),]$GPSmedian <- TRUE
        filedata$GPSmedian <- factor(filedata$GPSmedian,levels = c(TRUE,FALSE))
        GPSmedianUnion <- unique(c(GPSmedianUnion,filedata[filedata$GPSmedian==TRUE,'ID']))
        print(paste('GPSmedian',length(GPSmedianUnion),sep = ' '))      }
      print('RNA')
    }
  }
}

A <- data.frame(correlated.Pairs = c(length(exprimentUnion),
                                length(ScansiteUnion),length(MusiteDeepUnion),
                                length(NetworKINUnion),
                                length(GPShighUnion),
                                length(GPSmedianUnion)),
           row.names = c('Experiment','Scansite','MusiteDeep','NetworKIN','GPS(High)','GPS(Medium)'))
write.csv(A,paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/','Table2_predictor_correlated.csv',sep=''))