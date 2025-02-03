library(data.table)
library(stringr)
library(ggplot2)
##########Figure1 exp-confirmed Rho Number and Density############
setwd('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/spearman_clean/')
kinasedf <- read.csv('/home/wuyj2/eKPI/CancerOmics/kinase_basic(Sheet1).csv',header = T)
kinaseNames <- union(kinasedf$Offical_gene_symbol,kinasedf$KinBase_name)[-1]
Sample_THRD <- 6
Pvalue_THRD <- 0.05
Rvalue_THRD <- 0.2

filename <- list.files('.')
exprimentDF <- readRDS('/hwdata/home/wuyj2/eKPI/Prediction_Tool/Expriment.clean.rds')
exprimentDF$ID <- paste(exprimentDF$Kinase_names,exprimentDF$Pro.site,sep='#')

cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','PMID','TNType','Omics')] <- str_split_fixed(cancerFrame$fileName,'_',5)
cancerFrame$Omics <- str_replace_all(cancerFrame$Omics,'Kinase.','')

Rho_value_1 <- data.frame()
for(omicName in c('RNA')){#'Phos','Pro',
  fileDataFrame <- cancerFrame[cancerFrame$Omics==omicName,]
  kinase_Rho_value <- data.frame()
  Rho_value_1 <- data.frame()
  write_file <- paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/','Figure1F_correlated_Pos_',omicName,'V2.tsv',sep='')
  if(nrow(fileDataFrame)>0 & omicName=='Phos' & !file.exists(write_file)){
    for(i in 1:nrow(fileDataFrame)){
      print(fileDataFrame[i,'fileName'])
      filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                        header = T,data.table = F,fill=TRUE)
      filedata <- filedata[filedata$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$Rho.pvalue<=Pvalue_THRD&abs(filedata$Rho)>=Rvalue_THRD
      if(nrow(filedata)>0){
        filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
        filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
        filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
        filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
        filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
        filedata <- filedata[filedata$Rho>=Rvalue_THRD&filedata$p.adjust<Pvalue_THRD,]
        # filedata <- filedata[filedata$Rho>0,]
        # filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
        # colnames(filedata) <- c('ID','Rho')
        filedata$expriment <- FALSE
        if(length(intersect(filedata$ID,exprimentDF$ID))>0){
          filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
        }
        print(nrow(filedata))
        kinase_Rho_value <- rbind(kinase_Rho_value,data.frame(ID = filedata$KinaseName,
                                                              Rho = filedata$Rho,
                                                              expriment = filedata$expriment,
                                                              canerName = paste(fileDataFrame[i,'CancerName'],str_sub(fileDataFrame[i,'TNType'],1,1),sep='_')
                                                              )

                                  )
        # tmpdata <- filedata[filedata$expriment==TRUE&filedata$Rho==1,]
        # Rho_value_1 <- rbind(Rho_value_1,data.frame(RNA.names = tmpdata$RNA.names,
        #                                             GeneSite = tmpdata$GeneSite,
        #                                             Rho = tmpdata$Rho,
        #                                             # expriment = filedata$expriment,
        #                                             canerName = paste(fileDataFrame[i,'CancerName'],str_sub(fileDataFrame[i,'TNType'],1,1),sep='_'))
        # )
      }
    }
  }
  if(nrow(fileDataFrame)>0 & omicName=='Pro' & !file.exists(write_file)){
    for(i in 1:nrow(fileDataFrame)){
      print(fileDataFrame[i,'fileName'])
      filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[filedata$sample.num>Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$Rho.pvalue<=Pvalue_THRD&abs(filedata$Rho)>=Rvalue_THRD
      if(nrow(filedata)>0){
        filedata$KinaseName <- filedata$RNA.names
        filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
        filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
        filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
        filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
        filedata <- filedata[filedata$Rho>=Rvalue_THRD&filedata$p.adjust<Pvalue_THRD,]
        # filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
        # colnames(filedata) <- c('ID','Rho')
        filedata$expriment <- FALSE
        if(length(intersect(filedata$ID,exprimentDF$ID))>0){
          filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
        }
        print(nrow(filedata))
        kinase_Rho_value <- rbind(kinase_Rho_value,data.frame(ID = filedata$KinaseName,
                                                              Rho = filedata$Rho,
                                                              expriment = filedata$expriment,
                                                              canerName = paste(fileDataFrame[i,'CancerName'],str_sub(fileDataFrame[i,'TNType'],1,1),sep='_')))

      }
    }
  }
  if(nrow(fileDataFrame)>0 & omicName=='RNA' & !file.exists(write_file)){
    for(i in 1:nrow(fileDataFrame)){
      print(fileDataFrame[i,'fileName'])
      filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[filedata$sample.num>Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$Rho.pvalue<=Pvalue_THRD&abs(filedata$Rho)>=Rvalue_THRD
      if(nrow(filedata)>0){
        filedata$KinaseName <- filedata$RNA.names
        filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
        filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
        filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
        filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
        filedata <- filedata[filedata$Rho>=Rvalue_THRD&filedata$p.adjust<Pvalue_THRD,]
        # filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
        # colnames(filedata) <- c('ID','Rho')
        filedata$expriment <- FALSE
        if(length(intersect(filedata$ID,exprimentDF$ID))>0){
          filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
        }
        print(nrow(filedata))
        kinase_Rho_value <- rbind(kinase_Rho_value,data.frame(ID = filedata$KinaseName,
                                                              Rho = filedata$Rho,
                                                              expriment = filedata$expriment,
                                                              canerName = paste(fileDataFrame[i,'CancerName'],str_sub(fileDataFrame[i,'TNType'],1,1),sep='_')))
      }
    }
  }
  write.table(kinase_Rho_value,write_file,sep='\t',quote = F,row.names = F)
}

##########Figure1AB plot Rho Number and boxplot############
RNAdata <- fread(paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/','Figure1F_correlated_Pos_','RNA','V2.tsv',sep=''),
              header = T,sep='\t',data.table = F)
RNAdata <- RNAdata[!is.na(RNAdata$Rho),]
RNAdata$omicName <- 'RNA'

Prodata <- fread(paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/','Figure1F_correlated_Pos_','Pro','V2.tsv',sep=''),
              header = T,sep='\t',data.table = F)
Prodata <- Prodata[!is.na(Prodata$Rho),]
Prodata$omicName <- 'Pro'

Phosdata <- fread(paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/','Figure1F_correlated_Pos_','Phos','V2.tsv',sep=''),
              header = T,sep='\t',data.table = F)
Phosdata <- Phosdata[!is.na(Phosdata$Rho),]
Phosdata$omicName <- 'Phos'

data <- rbind(RNAdata,Prodata)
data <- rbind(data,Phosdata)
data$ID2 <- paste(data$canerName,data$omicName,sep = '_')
data$omicName <- factor(data$omicName,levels = c('Pro','Phos','RNA'))

ProRho <- median(data[data$omicName=='Pro','Rho'])
PhosRho <- median(data[data$omicName=='Phos','Rho'])
RNARho <- median(data[data$omicName=='RNA','Rho'])

g1 <- ggplot(data=data, aes(x=Rho,Group=omicName,fill = omicName))+
  geom_density(alpha=0.4)+
  # geom_rug()+
  # facet_grid(~canerName)+
  # scale_color_manual(values=c('#0C4E9B','#F98F34','#C72228'))+
  # scale_fill_manual(breaks = c('Expriment-confirmed','All'),values=c('#DD7C4F','#9D9EA3'))+
  theme_bw()+
  # geom_vline(xintercept = ProRho,linetype=2,color='#455a82') + # 添加base mean的水平线
  # geom_vline(xintercept = PhosRho,linetype=2,color='#2a6147') + # 添加base mean的水平线
  # geom_vline(xintercept = RNARho,linetype=2,color='#7c3031') + # 添加base mean的水平线
  scale_x_continuous(name = "Spearman' Rho")+
  scale_y_continuous(name = "Density")+
  scale_fill_manual(values=c('#2a6147','#F0A73A','#7c3031'),
                    breaks = c('Phos','Pro','RNA'),
                    labels = c( paste('Pro(Rho=',round(ProRho,3),')',sep=''),
                              paste('Phos(Rho=',round(PhosRho,3),')',sep=''),
                               paste('RNA(Rho=',round(RNARho,3),')',sep='')))+
  theme_bw(base_size = 18)+
  theme(panel.grid=element_blank(),
        legend.box = element_blank(),
        legend.position= "inside",
        legend.position.inside = c(0.7,0.9),
        panel.border = element_rect(fill=NA,linewidth = 1),
        # text=element_text(size=14,face="plain",color="black"),
        # axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=18,face="plain",color="black"),
        legend.text = element_text(size=18,face="plain",color="black"),
        legend.title = element_blank(),
        # legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank())
g1
ggsave(paste('/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure1B_Correlated_Pos_RNA_Pro_Phos_Density.pdf',sep=''),g1,width = 5,height = 5)

plotdata <- data.frame(table(data$ID2))
plotdata[,c('CancerName','TNType','OmicName')] <- str_split_fixed(plotdata$Var1,'_',3)
plotdata$OmicName <- factor(plotdata$OmicName,levels = c('Pro','Phos','RNA'))
mean(plotdata[plotdata$OmicName=='Phos','Freq'])
mean(plotdata[plotdata$OmicName=='Pro','Freq'])
mean(plotdata[plotdata$OmicName=='RNA','Freq'])

p1 <- ggplot(plotdata,aes(x=OmicName,y=Freq,fill=OmicName))+
  geom_bar(stat='summary',fun='mean',position = position_dodge(),width = 0.8)+
  # stat_compare_means(label = "p.format",method="wilcox.test")+
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge(.7))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                     labels = scales::scientific)+
  labs(x='',y='Number of Correlated Pairs')+
  scale_fill_manual(values=c('#455a82','#2a6147','#7c3031'))+
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        legend.box = element_blank(),
        legend.background = element_blank(),
        legend.title=element_blank(),
        legend.position = 'none',
        axis.title = element_text(size=20,face="plain",color="black"),
        axis.text = element_text(size=20,face="plain",color="black"),
        axis.text.x = element_text(size=20,face="plain",color="black"),#,angle = 45,vjust = 1,hjust = 1
        legend.text = element_text(size=20,face="plain",color="black"))
p1
ggsave(paste('/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure1A_Correlated_Pos_RNA_Pro_Phos_barplot.pdf',sep=''),p1,width = 5,height = 5)

# data2 <- data[data$expriment==TRUE,]
# plotdata <- data.frame(table(data2$ID2))
# plotdata[,c('CancerName','TNType','OmicName')] <- str_split_fixed(plotdata$Var1,'_',3)
# plotdata$OmicName <- factor(plotdata$OmicName,levels = c('Pro','Phos','RNA'))
# 
# p1 <- ggplot(plotdata,aes(x=OmicName,y=Freq,fill=OmicName))+
#   geom_bar(stat='summary',fun='median',position = position_dodge(),width = 0.8)+
#   # stat_compare_means(label = "p.format",method="wilcox.test")+
#   stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", colour = "black",
#                width = 0.15,position = position_dodge( .9))+
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))
#                      # labels = scales::scientific
#                      )+
#   labs(x='',y='Number of Exp-confirmed\nCorrelated Pairs')+
#   scale_fill_manual(values=c('#455a82','#2a6147','#7c3031'))+
#   theme_bw(base_size = 18) +
#   theme(panel.grid = element_blank(),
#         legend.box = element_blank(),
#         legend.background = element_blank(),
#         legend.title=element_blank(),
#         legend.position = 'none',
#         axis.text = element_text(size=18,face="plain",color="black"),
#         axis.text.x = element_text(size=18,face="plain",color="black"),#,angle = 45,vjust = 1,hjust = 1
#         legend.text = element_text(size=18,face="plain",color="black"))
# p1
# ggsave(paste('/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure1F_Exp_Correlated_Pos_RNA_Pro_Phos_barplot.pdf',sep=''),p1,width = 5,height = 5)

##########Figure1D Pro sig.pos.neg number#########
data_RNA <- fread('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/NewAll_RNA_Sig_Positive_Rho_KPS_allproject.csv',header = T,sep=',',data.table = F,quote = F)
data_RNA$OmicName <- 'RNA'
data_Pro <- fread('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/NewAll_Pro_Sig_Positive_Rho_KPS_allproject.csv',header = T,sep=',',data.table = F,quote = F)
data_Pro$OmicName <- 'Pro'
data_Phos <- fread('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/NewAll_Phos_Sig_Positive_Rho_KPS_allproject.csv',header = T,sep=',',data.table = F,quote = F)
data_Phos$OmicName <- 'Phos'

data <- rbind(data_RNA,data_Pro)
data <- rbind(data,data_Phos)
data$TNType <- factor(data$TNType,levels=c('Normal','Tumor'))
data$ID <- paste(data$CancerName,data$TNType,data$OmicName,sep='_')

plotdata <- data.frame(table(data$ID))
plotdata[,c('CancerName','TNType','OmicName')] <- str_split_fixed(plotdata$Var,'_',3)
plotdata$TNType <- factor(plotdata$TNType,levels=c('Normal','Tumor'))
plotdata$OmicName <- factor(plotdata$OmicName,levels=c('Pro','Phos','RNA'))
p1 <- ggplot(data=plotdata,aes(x=OmicName,y=Freq,fill=TNType))+
  # geom_bar(stat = 'identity',width = 0.8,position = 'dodge')+
  geom_bar(stat='summary',fun='mean',position = position_dodge(),width = 0.8)+
  # stat_compare_means(label = "p.format",method="wilcox.test")+
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge(.7))+
  stat_compare_means(label = "p.signif",
                     method="wilcox.test",
                     # comparisons = list(c('Normal','Tumor')),
                     # method.args = list(alternative = "less"),
                     size=6)+
  # geom_text(aes(label=FreqLabel),size=4,
  #           position = position_dodge(width = 0.8), #相应的注释宽度也调整
  #           vjust=-0.3)+    #调节注释高度
  # scale_y_break(c(18000,45000),#截断位置及范围
  #               space = 0.3,#间距大小
  #               scales = 1.2)+#上下显示比例，大于1上面比例大，小于1下面比例大
  scale_fill_manual(values=c('#455a82','#7c3031'))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                     labels = scales::scientific)+
  labs(x="",y="Number of Correlated Pairs")+
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        legend.position = 'nonr',
        # legend.position.inside = c(0.8,0.8),
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=18,face="plain",color="black"))
p1

omicsnames <- unique(data$OmicName)
testrst <- data.frame()
for(omics in omicsnames){
  tmpdata <- data[data$OmicName==omics,]
  tmptest <- wilcox.test(tmpdata[tmpdata$TNType=='Normal','Rho'],tmpdata[tmpdata$TNType=='Tumor','Rho'],alternative='greater')
  testrst <- rbind(testrst,data.frame(OmicName=omics,pValue=tmptest[['p.value']]))
}
# tmpdata$CancerName <- factor(tmpdata$CancerName)
testrst$siglabel <- ifelse(testrst$pValue<0.001,'***',
                           ifelse(testrst$pValue<0.01,'**',
                                  ifelse(testrst$pValue<0.05,'*','ns')))
data$OmicName <- factor(data$OmicName,levels=c('Pro','Phos','RNA'))
testrst$OmicName <- factor(testrst$OmicName,levels=c('Pro','Phos','RNA'))


p2 <- ggplot(data=data,aes(x=OmicName,y=Rho,fill=TNType))+
  geom_boxplot(outlier.shape = NA)+
  # stat_compare_means(label = "p.signif",
  #                    method="wilcox.test",
  #                    # comparisons = list(c('Normal','Tumor')),
  #                    method.args = list(alternative = "less"),
  #                    size=6)+
  geom_text(data = testrst,
            # max.overlaps = getOption("ggrepel.max.overlaps", default = 3),
            aes(x =OmicName,y=0.7, label = siglabel),
            size = 6,
            color = 'black')+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  scale_fill_manual(values=c('#455a82','#7c3031'))+
  labs(x="",y="Spearman' Rho")+
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        legend.box = element_blank(),
        legend.background = element_blank(),
        legend.title=element_blank(),
        legend.position = 'none',
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.text.x = element_text(size=18,face="plain",color="black",angle = 45,vjust = 1,hjust = 1),
        legend.text = element_text(size=18,face="plain",color="black"))

# p1
design <- "111
           ###
           222"
p_combined <- wrap_plots(list(p1,p2),
                         heights = c(4.5,0,4.5), design = design,nrow = 2)
ggsave('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure1D_Liu_levels_Sig_Pos_project_boxplot_barplot.pdf',p_combined,width=4,height=8)




##########Figure1E Pro sig.pos KPS pairs################
data_agg <- fread('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/KPS_exprimentV2/NewAll_Pro_Sig_Positive_Rho_KPS_allproject.csv',header = T,sep=',',data.table = F,quote = F)
data_agg <- data_agg[!is.na(data_agg$Rho),]
data_agg$TNType <- factor(data_agg$TNType,levels=c('Normal','Tumor'))
data_agg$CancerName <- factor(data_agg$CancerName,levels = c('CCRCC','COC','DTGC','ECA','ENCA','EOGC','ESCC','ESHCC','HBVHCC','HNSCC','GBM','LSCC','LUAD.CN','LUAD.US','PDAC'))
data_agg$ID <- paste(data_agg$CancerName,data_agg$TNType,sep='_')

data_count <- data.frame(table(data_agg$ID))
data_count$FreqLabel <- ifelse(data_count$Freq<10000,as.numeric(data_count$Freq),format(as.numeric(data_count$Freq),scientific=TRUE, digits = 2L))
data_count[,c('CancerName','TNType')] <- str_split_fixed(data_count$Var1,'_',2)
data_count$CancerName <- factor(data_count$CancerName,levels = c('CCRCC','COC','DTGC','ECA','ENCA','EOGC','ESCC','ESHCC','HBVHCC','HNSCC','GBM','LSCC','LUAD.CN','LUAD.US','PDAC'))
data_count$TNType <- factor(data_count$TNType,levels=c('Normal','Tumor'))

p1 <- ggplot(data=data_count,aes(x=CancerName,y=Freq,fill=TNType))+
  geom_bar(stat = 'identity',width = 0.8,position = 'dodge')+
  # geom_text(aes(label=FreqLabel),size=4,
  #           position = position_dodge(width = 0.8), #相应的注释宽度也调整
  #           vjust=-0.3)+    #调节注释高度
  # scale_y_break(c(18000,45000),#截断位置及范围
  #               space = 0.3,#间距大小
  #               scales = 1.2)+#上下显示比例，大于1上面比例大，小于1下面比例大
  scale_fill_manual(values=c('#455a82','#7c3031'))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                     labels = scales::scientific)+
  labs(x="",y="Number of Correlated Pairs")+
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        legend.box = element_blank(),
        legend.background = element_blank(),
        legend.title=element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=18,face="plain",color="black"))
# p1
# ggsave('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure1D_Liu_Pro_Sig_Pos_project_number_barplot.pdf',p1,width=10,height=3.5)

p2 <- ggplot(data=data_agg,aes(x=CancerName,y=Rho,fill=TNType))+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(label = "p.signif",
                     method="wilcox.test",
                     # comparisons = list(c('Normal','Tumor')),
                     method.args = list(alternative = "less"),
                     size=6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  scale_fill_manual(values=c('#455a82','#7c3031'))+
  labs(x="",y="Spearman' Rho")+
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        legend.box = element_blank(),
        legend.background = element_blank(),
        legend.title=element_blank(),
        legend.position = 'none',
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.text.x = element_text(size=18,face="plain",color="black",angle = 45,vjust = 1,hjust = 1),
        legend.text = element_text(size=18,face="plain",color="black"))
design <- "111
           ###
           222"
p_combined <- wrap_plots(list(p1,p2),
                         heights = c(4.5,0,4.5), design = design,nrow = 2)
ggsave('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure1D_Liu_Pro_Sig_Pos_project_boxplot_barplot.pdf',p_combined,width=10,height=8)
print('Pro_Sig_Pos')


##########Figure1F################
setwd('/hwdata//home/wuyj2/eKPI/ArticleCancerStat/Spearman_RRAV2/')
filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','TNType')] <- str_split_fixed(str_replace_all(cancerFrame$fileName,'.tsv',''),'_',2)
finaldata <- data.frame()
for(i in 1:nrow(cancerFrame)){
  filedata <- fread(cancerFrame[i,'fileName'],header = T,sep='\t',data.table = F)
  filedata <- filedata[filedata$expriment==TRUE,]
  if(ncol(filedata)==11){
    Phosnum <- length(na.omit(filedata$Phos.Rank))
    RNAnum <- length(na.omit(filedata$RNA.Rank))
    Pronum <- length(na.omit(filedata$Pro.Rank))

    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            KPSNUM = c(Phosnum,Pronum,RNAnum),
                                            OmicName = c('Phos','Pro','RNA')))
  }
  if(ncol(filedata)==9){
    Phosnum <- length(na.omit(filedata$Phos.Rank))
    # RNAnum <- length(na.omit(filedata$RNA.Rank))
    Pronum <- length(na.omit(filedata$Pro.Rank))

    finaldata <- rbind(finaldata,data.frame(cancerName = cancerFrame[i,'CancerName'],
                                            TNType = cancerFrame[i,'TNType'],
                                            KPSNUM = c(Phosnum,Pronum),
                                            OmicName = c('Phos','Pro')))
  }
}
finaldata$OmicName <- factor(finaldata$OmicName,levels = c('RNA','Phos','Pro'))

mean(finaldata[finaldata$OmicName=='Phos','KPSNUM'])
mean(finaldata[finaldata$OmicName=='Pro','KPSNUM'])
mean(finaldata[finaldata$OmicName=='RNA','KPSNUM'])

p1 <- ggplot(data=finaldata,aes(x=KPSNUM,y=OmicName,fill=OmicName))+
  # geom_bar(stat = 'identity',width = 0.8,position = 'dodge')+
  geom_bar(stat='summary',fun='mean',position = position_dodge(),width = 0.8)+
  # stat_compare_means(label = "p.format",method="wilcox.test")+
  stat_summary(fun.data = 'mean_cl_boot', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge(.8))+
  # stat_compare_means(label = "p.signif",
  #                    method="wilcox.test",
  #                    comparisons = list(c('Pro','Phos'),
  #                                       c('Phos','RNA'),
  #                                       c('Pro','RNA')),
  #                    label.y = c(600,650,700),
  #                    # method.args = list(alternative = "less"),
  #                    size=6)+
  # geom_text(aes(label=FreqLabel),size=4,
  #           position = position_dodge(width = 0.8), #相应的注释宽度也调整
  #           vjust=-0.3)+    #调节注释高度
  # scale_y_break(c(18000,45000),#截断位置及范围
  #               space = 0.3,#间距大小
  #               scales = 1.2)+#上下显示比例，大于1上面比例大，小于1下面比例大
  scale_fill_manual(values=c('#7c3031','#2a6147','#455a82'))+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  labs(y="",x="Number of Exp-confirmed\nCorrelated Pairs")+
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        legend.position = 'nonr',
        # legend.position.inside = c(0.8,0.8),
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.text.x =  element_text(size=18,face="plain",color="black"),#,angle = 45,vjust = 1,hjust = 1
        # axis.ticks.x = element_blank(),
        legend.text = element_text(size=18,face="plain",color="black"))
p1
ggsave('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/FigureV2/Figure1F_Liu_exp_confirmed_KPS_barplot.pdf',p1,width=4,height=2.5)

##########Figure1H Fisher test for cor.pairs and expriment/Tools############
cal_fishertest <- function(filedata,Project,type,Tool,omics,f){
  tmptab <- xtabs(f,data=filedata)
  tmpfisherst <- fisher.test(tmptab)
  exp.cor <- tmptab[1,1]
  exp.Uncor <- tmptab[2,1]
  Noexp.cor <- tmptab[1,2]
  Noexp.Uncor <- tmptab[2,2]
  Oddr <- tmpfisherst$estimate
  pvalue <- tmpfisherst$p.value
  OddrStr <- paste(round(tmpfisherst$estimate,3),'(',round(tmpfisherst$conf.int[1],3),',',round(tmpfisherst$conf.int[2],3),')',sep='')

  tmpfishertest <- data.frame(
    Project = Project,
    Omics = omics,
    TNType = type,
    Proof.Cor=exp.cor,
    Proof.UnCor=exp.Uncor,
    NonProof.Cor = Noexp.cor,
    NonProof.Uncor = Noexp.Uncor,
    Oddsr = Oddr,
    pvalue = pvalue,
    OddInter = OddrStr,
    ExternalProof = Tool)

  return(tmpfishertest)
}

setwd('/home/liangzr/eKPI/0.data/spearman/')
filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','PMID','TNType','Omics')] <- str_split_fixed(cancerFrame$fileName,'_',5)
cancerFrame$Omics <- str_replace_all(cancerFrame$Omics,'Kinase.','')

kinasedf <- read.csv('/home/wuyj2/eKPI/CancerOmics/kinase_basic(Sheet1).csv',header = T)
kinaseNames <- union(kinasedf$Offical_gene_symbol,kinasedf$KinBase_name)[-1]

Sample_THRD <- 6
Pvalue_THRD <- 0.05
Rvalue_THRD <- 0.2

print('Fisher test for correlated and expriment/Tools')
GPShighDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_High.rds')
GPShighDF$ID <- paste(GPShighDF$kinaseName,GPShighDF$Pro.Site,sep=',')
print('GPShighDF')

GPSmedianDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/GPS_Medium.rds')
GPSmedianDF$ID <- paste(GPSmedianDF$kinaseName,GPSmedianDF$Pro.Site,sep=',')
print('GPSmedianDF')

exprimentDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/Expriment.clean.rds')
exprimentDF$ID <- paste(exprimentDF$Kinase_names,exprimentDF$Pro.site,sep=',')

ScansiteDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/Scansite.clean.rds')
ScansiteDF$ID <- paste(ScansiteDF$Kinase_names,ScansiteDF$Pro.site,sep=',')

MusiteDeepDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/MusiteDeep.clean.rds')
MusiteDeepDF$ID <- paste(MusiteDeepDF$Kinase_names,MusiteDeepDF$Pro.site,sep=',')
print('MusiteDeepDF')

NetworKINDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/NetworKIN.clean.rds')
NetworKINDF$ID <- paste(NetworKINDF$Kinase_names,NetworKINDF$Pro.site,sep=',')
print('NetworKINDF')

for(omicName in c('Phos','Pro','RNA')){
  fileDataFrame <- cancerFrame[cancerFrame$Omics==omicName,]
  fishertestDF <- data.frame()
  write_file <- paste('/home/wuyj2/eKPI/ArticleCancerStat/KPS_expriment/','Figure1H_Uncor_Correlated_Pos_fishertest',omicName,'.tsv',sep='')
  if(nrow(fileDataFrame)>0 & omicName=='Phos' & !file.exists(write_file)){
    for(i in 1:nrow(fileDataFrame)){
      filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[filedata$sample.num>=Sample_THRD&filedata$Rho>0,]#根据样本数量过滤满足条件的KPS pair
      filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
      filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep=',')

      filedata$CorLabel <- ifelse(filedata$Rho>=Rvalue_THRD&filedata$Rho.pvalue<Pvalue_THRD,'Cor.','Uncor.')#根据Rho、pvalue值筛选显著相关的KPS pair
      filedata$CorLabel <- factor(filedata$CorLabel,levels = c('Cor.','Uncor.'))

      filedata$expriment <- FALSE
      filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
      filedata$expriment <- factor(filedata$expriment,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+expriment)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'expriment','phos',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$Scansite <- FALSE
      filedata[which(filedata$ID %in% ScansiteDF$ID),]$Scansite <- TRUE
      filedata$Scansite <- factor(filedata$Scansite,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+Scansite)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'Scansite','Phos',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$MusiteDeep <- FALSE
      filedata[which(filedata$ID %in% MusiteDeepDF$ID),]$MusiteDeep <- TRUE
      filedata$MusiteDeep <- factor(filedata$MusiteDeep,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+MusiteDeep)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'MusiteDeep','Phos',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$NetworKIN <- FALSE
      filedata[which(filedata$ID %in% NetworKINDF$ID),]$NetworKIN <- TRUE
      filedata$NetworKIN <- factor(filedata$NetworKIN,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+NetworKIN)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'NetworKIN','Phos',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$GPShigh <- FALSE
      filedata[which(filedata$ID %in% GPShighDF$ID),]$GPShigh <- TRUE
      filedata$GPShigh <- factor(filedata$GPShigh,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+GPShigh)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'GPShigh','Phos',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$GPSmedian <- FALSE
      filedata[which(filedata$ID %in% GPSmedianDF$ID),]$GPSmedian <- TRUE
      filedata$GPSmedian <- factor(filedata$GPSmedian,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+GPSmedian)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'GPSmedian','Phos',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

    }
  }
  if(nrow(fileDataFrame)>0 & omicName=='Pro' & !file.exists(write_file)){
    for(i in 1:nrow(fileDataFrame)){
      filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[filedata$sample.num>=2,]#根据样本数量过滤满足条件的KPS pair
      filedata$KinaseName <- filedata$RNA.names
      filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep=',')

      filedata$CorLabel <- ifelse(filedata$Rho>=Rvalue_THRD&filedata$Rho.pvalue<Pvalue_THRD,'Cor.','Uncor.')#根据Rho、pvalue值筛选显著相关的KPS pair
      filedata$CorLabel <- factor(filedata$CorLabel,levels = c('Cor.','Uncor.'))

      filedata$expriment <- FALSE
      filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
      filedata$expriment <- factor(filedata$expriment,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+expriment)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'expriment','Pro',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$Scansite <- FALSE
      filedata[which(filedata$ID %in% ScansiteDF$ID),]$Scansite <- TRUE
      filedata$Scansite <- factor(filedata$Scansite,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+Scansite)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'Scansite','Pro',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$MusiteDeep <- FALSE
      filedata[which(filedata$ID %in% MusiteDeepDF$ID),]$MusiteDeep <- TRUE
      filedata$MusiteDeep <- factor(filedata$MusiteDeep,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+MusiteDeep)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'MusiteDeep','Pro',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$NetworKIN <- FALSE
      filedata[which(filedata$ID %in% NetworKINDF$ID),]$NetworKIN <- TRUE
      filedata$NetworKIN <- factor(filedata$NetworKIN,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+NetworKIN)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'NetworKIN','Pro',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$GPShigh <- FALSE
      filedata[which(filedata$ID %in% GPShighDF$ID),]$GPShigh <- TRUE
      filedata$GPShigh <- factor(filedata$GPShigh,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+GPShigh)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'GPShigh','Pro',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$GPSmedian <- FALSE
      filedata[which(filedata$ID %in% GPSmedianDF$ID),]$GPSmedian <- TRUE
      filedata$GPSmedian <- factor(filedata$GPSmedian,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+GPSmedian)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'GPSmedian','Pro',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

    }
  }
  if(nrow(fileDataFrame)>0 & omicName=='RNA' & !file.exists(write_file)){
    for(i in 1:nrow(fileDataFrame)){
      filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                        header = T,data.table = F)
      filedata <- filedata[filedata$sample.num>=2,]#根据样本数量过滤满足条件的KPS pair
      filedata$KinaseName <- filedata$RNA.names
      filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
      filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
      filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
      filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep=',')

      filedata$CorLabel <- ifelse(filedata$Rho>=Rvalue_THRD&filedata$Rho.pvalue<Pvalue_THRD,'Cor.','Uncor.')#根据Rho、pvalue值筛选显著相关的KPS pair
      filedata$CorLabel <- factor(filedata$CorLabel,levels = c('Cor.','Uncor.'))

      filedata$expriment <- FALSE
      filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
      filedata$expriment <- factor(filedata$expriment,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+expriment)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'expriment','RNA',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$Scansite <- FALSE
      filedata[which(filedata$ID %in% ScansiteDF$ID),]$Scansite <- TRUE
      filedata$Scansite <- factor(filedata$Scansite,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+Scansite)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'Scansite','RNA',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$MusiteDeep <- FALSE
      filedata[which(filedata$ID %in% MusiteDeepDF$ID),]$MusiteDeep <- TRUE
      filedata$MusiteDeep <- factor(filedata$MusiteDeep,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+MusiteDeep)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'MusiteDeep','RNA',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$NetworKIN <- FALSE
      filedata[which(filedata$ID %in% NetworKINDF$ID),]$NetworKIN <- TRUE
      filedata$NetworKIN <- factor(filedata$NetworKIN,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+NetworKIN)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'NetworKIN','RNA',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$GPShigh <- FALSE
      filedata[which(filedata$ID %in% GPShighDF$ID),]$GPShigh <- TRUE
      filedata$GPShigh <- factor(filedata$GPShigh,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+GPShigh)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'GPShigh','RNA',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)

      filedata$GPSmedian <- FALSE
      filedata[which(filedata$ID %in% GPSmedianDF$ID),]$GPSmedian <- TRUE
      filedata$GPSmedian <- factor(filedata$GPSmedian,levels = c(TRUE,FALSE))

      f <- as.formula(~CorLabel+GPSmedian)
      tmpfisherDF <- cal_fishertest(filedata,fileDataFrame[i,'CancerName'],fileDataFrame[i,'TNType'],'GPSmedian','RNA',f)
      fishertestDF <- rbind(fishertestDF,tmpfisherDF)    }
  }

  write.table(fishertestDF,write_file,sep='\t',quote = F,row.names = F)
}



