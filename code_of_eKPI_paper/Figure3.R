library(data.table)
library(stringr)
library(patchwork)
library(cowplot)
##########Figure3A Kinase trueback comp FDR LiuData#############
setwd('/hwdata/home/wuyj2/eKPI/spearman/')
sample_Thrd <- 6

filename <- list.files('.')
cancerFrame <- data.frame(fileName=filename)
cancerFrame[,c('CancerName','PMID','TNType','Omics')] <- str_split_fixed(cancerFrame$fileName,'_',5)
cancerFrame$Omics <- str_replace_all(cancerFrame$Omics,'Kinase.','')

exprimentDF <- readRDS('/home/wuyj2/eKPI/Prediction_Tool/Expriment.clean.rds')
exprimentDF$ID <- paste(exprimentDF$Kinase_names,exprimentDF$Pro.site,sep='#')

kinasedf <- read.csv('/home/wuyj2/eKPI/CancerOmics/kinase_basic(Sheet1).csv',header = T)
kinaseNames <- union(kinasedf$Offical_gene_symbol,kinasedf$KinBase_name)[-1]
Sample_THRD <- 6
Pvalue_THRD <- 0.05
Rvalue_THRD <- 0.2
print('Kinase and Expriment-proofed for each project')

for(omicName in c('RNA','Pro','Phos')){#'Phos',
    fileDataFrame <- cancerFrame[cancerFrame$Omics==omicName,]
    kinase_Rho_value <- data.frame()
    write_file <- paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Kinase_rst/Kinase_Project/','Correlated_pvalue_FDR_exp_Pos_wilcox_',omicName,'.tsv',sep='')
    if(nrow(fileDataFrame)>0 & omicName=='Phos' & !file.exists(write_file)){
      for(i in 1:nrow(fileDataFrame)){
        filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                          header = T,data.table = F)
        filedata <- filedata[!is.na(filedata$Rho),]
        print(fileDataFrame[i,'fileName'])
        max_sample <- 217
        # if(max_sampleNUM>samp_cutoff){
        #   Pvalue_THRD <- 0.01
        # }else{
        #   Pvalue_THRD <- 0.05
        # }
        filedata <- filedata[filedata$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$p.adjust<=Pvalue_THRD&filedata$Rho>=Rvalue_THRD
        if(nrow(filedata)>0){
          filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
          filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
          filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
          filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
          filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
          # filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
          # colnames(filedata) <- c('ID','Rho')
          filedata$expriment <- FALSE
          filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
          filedata2 <- filedata[filedata$expriment==TRUE,]
          tmpkinases <- unique(filedata2$KinaseName)

          for(kiName in tmpkinases){
            print(kiName)
            tmpdata <- filedata[filedata$KinaseName==kiName,]
            tmpdata$FDR <- p.adjust(tmpdata$Rho.pvalue,method = 'BH')

            k1 <- max_sample^3-max_sample
            tmpdata$p.adjust <- ((k1-(tmpdata$sample.num^3-tmpdata$sample.num))*tmpdata$Rho.pvalue+(tmpdata$sample.num^3-tmpdata$sample.num)*tmpdata$FDR)/k1

            tmpdata2 <- tmpdata[tmpdata$sample.num>=Sample_THRD&tmpdata$p.adjust<=Pvalue_THRD&tmpdata$Rho>=Rvalue_THRD,]

            proved_data <- tmpdata2[tmpdata2$expriment==TRUE,]
            if(nrow(proved_data)>0){
              pvalue <- wilcox.test(proved_data$Rho,tmpdata2$Rho,alternative='greater')$p.value

              proved_num <- nrow(proved_data)
              proved_Median_Rho <- median(proved_data$Rho)
              back_Median_Rho <- median(tmpdata$Rho)
              back_num <- nrow(tmpdata)
              kinase_Rho_value <- rbind(kinase_Rho_value,
                                        data.frame(Project = fileDataFrame[i,'CancerName'],
                                                   Kinase = kiName,
                                                   Level = 'Phos',
                                                   proved_num = proved_num,
                                                   proved_Median_Rho = proved_Median_Rho,
                                                   back_num = back_num,
                                                   back_Median_Rho = back_Median_Rho,
                                                   Pvalue = pvalue,
                                                   TNtype = fileDataFrame[i,'TNType']))
            }
          }
        }
      }
    }
    if(nrow(fileDataFrame)>0 & omicName=='Pro' & !file.exists(write_file)){
      for(i in 1:nrow(fileDataFrame)){
        filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                          header = T,data.table = F)
        filedata <- filedata[!is.na(filedata$Rho),]
        filedata <- filedata[filedata$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#
        # filedata <- filedata[filedata$sample.num>=Sample_THRD&filedata$p.adjust<=0.01&filedata$Rho>=Rvalue_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#
        max_sample <- max(c(217,max(filedata$sample.num)))
        print(fileDataFrame[i,'fileName'])
        if(nrow(filedata)>0){
          filedata$KinaseName <- filedata$RNA.names
          filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
          filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
          filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
          filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
          # filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
          # colnames(filedata) <- c('ID','Rho')
          filedata$expriment <- FALSE
          if(length(intersect(filedata$ID,exprimentDF$ID))>0){
            filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
          }
          filedata2 <- filedata[filedata$expriment==TRUE,]
          tmpkinases <- unique(filedata2$KinaseName)

          for(kiName in tmpkinases){
            print(kiName)
            tmpdata <- filedata[filedata$KinaseName==kiName,]
            tmpdata$FDR <- p.adjust(tmpdata$Rho.pvalue,method = 'BH')
            k1 <- max_sample^3-max_sample
            tmpdata$p.adjust <- ((k1-(tmpdata$sample.num^3-tmpdata$sample.num))*tmpdata$Rho.pvalue+(tmpdata$sample.num^3-tmpdata$sample.num)*tmpdata$FDR)/k1


            tmpdata2 <- tmpdata[tmpdata$sample.num>=Sample_THRD&tmpdata$p.adjust<=Pvalue_THRD&tmpdata$Rho>=Rvalue_THRD,]

            proved_data <- tmpdata2[tmpdata2$expriment==TRUE,]

            if(nrow(proved_data)>0){
              pvalue <- wilcox.test(proved_data$Rho,tmpdata2$Rho,alternative='greater')$p.value
              proved_num <- nrow(proved_data)
              proved_Median_Rho <- median(proved_data$Rho)
              back_Median_Rho <- median(tmpdata$Rho)
              back_num <- nrow(tmpdata)
              kinase_Rho_value <- rbind(kinase_Rho_value,
                                        data.frame(Project = fileDataFrame[i,'CancerName'],
                                                   Kinase = kiName,
                                                   Level = 'Pro',
                                                   proved_num = proved_num,
                                                   proved_Median_Rho = proved_Median_Rho,
                                                   back_num = back_num,
                                                   back_Median_Rho = back_Median_Rho,
                                                   Pvalue = pvalue,
                                                   TNtype = fileDataFrame[i,'TNType']))
            }
          }
        }
      }
    }
    if(nrow(fileDataFrame)>0 & omicName=='RNA' & !file.exists(write_file)){
      for(i in 1:nrow(fileDataFrame)){
        print(fileDataFrame[i,'fileName'])
        filedata <- fread(paste(fileDataFrame[i,'fileName'],sep=''),
                          header = T,data.table = F)
        filedata <- filedata[!is.na(filedata$Rho),]
        filedata <- filedata[filedata$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$p.adjust<=Pvalue_THRD&filedata$Rho>=Rvalue_THRD
        max_sample <- max(c(217,max(filedata$sample.num)))
        if(nrow(filedata)>0){
          filedata$KinaseName <- filedata$RNA.names
          filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
          filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
          filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
          filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
          # filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
          # colnames(filedata) <- c('ID','Rho')
          filedata$expriment <- FALSE
          if(length(intersect(filedata$ID,exprimentDF$ID))>0){
            filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
          }
          filedata2 <- filedata[filedata$expriment==TRUE,]
          tmpkinases <- unique(filedata2$KinaseName)

          for(kiName in tmpkinases){
            tmpdata <- filedata[filedata$KinaseName==kiName,]
            tmpdata$FDR <- p.adjust(tmpdata$Rho.pvalue,method = 'BH')
            k1 <- max_sample^3-max_sample
            tmpdata$p.adjust <- ((k1-(tmpdata$sample.num^3-tmpdata$sample.num))*tmpdata$Rho.pvalue+(tmpdata$sample.num^3-tmpdata$sample.num)*tmpdata$FDR)/k1


            tmpdata2 <- tmpdata[tmpdata$sample.num>=Sample_THRD&tmpdata$p.adjust<=Pvalue_THRD&tmpdata$Rho>=Rvalue_THRD,]

            proved_data <- tmpdata2[tmpdata2$expriment==TRUE,]
            if(nrow(proved_data)>0){
              print('Expriment data >0')
              pvalue <- wilcox.test(proved_data$Rho,tmpdata2$Rho,alternative='greater')$p.value
              proved_num <- nrow(proved_data)
              proved_Median_Rho <- median(proved_data$Rho)
              back_Median_Rho <- median(tmpdata$Rho)
              back_num <- nrow(tmpdata)
              kinase_Rho_value <- rbind(kinase_Rho_value,
                                        data.frame(Project = fileDataFrame[i,'CancerName'],
                                                   Kinase = kiName,
                                                   Level = 'RNA',
                                                   proved_num = proved_num,
                                                   proved_Median_Rho = proved_Median_Rho,
                                                   back_num = back_num,
                                                   back_Median_Rho = back_Median_Rho,
                                                   Pvalue = pvalue,
                                                   TNtype = fileDataFrame[i,'TNType']))
            }
          }
        }
      }
    }

    write.table(kinase_Rho_value,write_file,sep='\t',quote = F,row.names = F)
  }


##########Figure3C CDC7###############
fileDataFrame <- cancerFrame[cancerFrame$CancerName=='HNSCC',]
filedata <- fread(paste(fileDataFrame[3,'fileName'],sep=''),
                  header = T,data.table = F)
filedata <- filedata[filedata$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$Rho.pvalue<=Pvalue_THRD&filedata$Rho>=Rvalue_THRD
filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
# filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
# colnames(filedata) <- c('ID','Rho')
filedata$expriment <- FALSE
filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
filedata$Project <- 'HNSCC_N'

filedata2 <- fread(paste(fileDataFrame[6,'fileName'],sep=''),
                   header = T,data.table = F)
filedata2 <- filedata2[filedata2$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata2$Rho.pvalue<=Pvalue_THRD&filedata2$Rho>=Rvalue_THRD
filedata2$KinaseName <- str_split_fixed(filedata2$RNA.names,'_',2)[,1]
filedata2 <- filedata2[which(filedata2$KinaseName %in% kinaseNames),]#筛选激酶
filedata2$SubstrateName <- str_split_fixed(filedata2$GeneSite,'#',3)[,2]
filedata2 <- filedata2[filedata2$KinaseName!=filedata2$SubstrateName,]#过滤自磷酸化的KPS pairs
filedata2$ID  <- paste(filedata2$KinaseName,filedata2$GeneSite,sep='#')
# filedata2 <- aggregate(filedata2[,c('Rho')],list(filedata2$ID),max,na.rm=T)
# colnames(filedata2) <- c('ID','Rho')
filedata2$expriment <- FALSE
filedata2[which(filedata2$ID %in% exprimentDF$ID),]$expriment <- TRUE
filedata2$Project <- 'HNSCC_T'

filedata3 <- rbind(filedata,filedata2)
filedata4 <- filedata3[filedata3$KinaseName=='CDC7'&filedata3$Rho.pvalue<=Pvalue_THRD&filedata3$Rho>=Rvalue_THRD,]

tmpdata <- filedata4[filedata4$Project=='HNSCC_T',]
all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
all_median
exp_median
wilcoxtest <- wilcox.test(tmpdata[tmpdata$expriment==TRUE,]$Rho,tmpdata[tmpdata$expriment==FALSE,]$Rho,alternative = 'greater')
pvalue <- format(wilcoxtest[["p.value"]], scientific = T,digits = 3)

p1 <- ggplot(tmpdata,aes(x=Rho,Group=expriment,fill=expriment))+
  geom_density(alpha=0.3)+
  theme_bw()+
  facet_wrap(~Project)+
  geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
  geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
  scale_x_continuous(name = "")+
  scale_y_continuous(name = "Density")+
  scale_fill_manual(values=c('#B54764','#9D9EA3'),
                    breaks = c(TRUE,FALSE),
                    labels = c('Exp-confirmed','All'))+
  annotate("text",label=paste('p = ',pvalue,sep=''),x=0.6,y=4,size=6)+
  theme_classic(base_size = 18)+
  theme(panel.grid=element_blank(),
        legend.box = element_blank(),
        legend.position= "inside",
        legend.position.inside = c(0.6,0.87),
        panel.border = element_rect(fill=NA,linewidth = 1),
        # text=element_text(size=14,face="plain",color="black"),
        # axis.title=element_text(size=16,face="plain",color="black"),
        # axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_blank(),
        # legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank())
# p1

tmpdata <- filedata4[filedata4$Project=='HNSCC_N',]
all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
all_median
exp_median
wilcoxtest <- wilcox.test(tmpdata[tmpdata$expriment==TRUE,]$Rho,tmpdata[tmpdata$expriment==FALSE,]$Rho,alternative = 'greater')
pvalue <- format(wilcoxtest[["p.value"]], scientific = T,digits = 3)

p2 <- ggplot(tmpdata,aes(x=Rho,Group=expriment,fill=expriment))+
  geom_density(alpha=0.3)+
  theme_bw()+
  facet_wrap(~Project)+
  geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
  geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
  scale_x_continuous(name = "Spearman' Rho")+
  scale_y_continuous(name = "Density")+
  scale_fill_manual(values=c('#B54764','#9D9EA3'),
                    breaks = c(TRUE,FALSE),
                    labels = c('Exp-confirmed','All'))+
  annotate("text",label=paste('p = ',pvalue,sep=''),x=0.73,y=4.5,size=6)+
  theme_classic(base_size = 18)+
  theme(panel.grid=element_blank(),
        legend.box = element_blank(),
        legend.position= "none",
        panel.border = element_rect(fill=NA,linewidth = 1),
        # text=element_text(size=14,face="plain",color="black"),
        # axis.title=element_text(size=16,face="plain",color="black"),
        # axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_blank(),
        # legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank())

# p <- p1 %>% insert_bottom(p2)

p <- plot_grid(plotlist = list(p2,NULL,p1),ncol=1,
               rel_heights = c(1,0.0000001,1),
               rel_widths = c(1))
# p
ggsave(paste('/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Figure/Figure3B_CDC7_HNSCC_RNA_density','.pdf',sep=''),p,width =4,height = 8)


#write fastq
tmpdata <- filedata4[filedata4$Project=='HNSCC_T',]
# all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
# all_median
exp_median

subsiteDF <- aggregate(tmpdata[,c('Rho','sample.num')],list(tmpdata$GeneSite),max,na.rm=T)
colnames(subsiteDF)[1] <- 'GeneSite'

subsiteDF[,c('UniprotID','GeneName','PTM')] <- str_split_fixed(subsiteDF$GeneSite,'#',3)
subsiteDF$PTMtype <- str_sub(subsiteDF$PTM,1,1)
subsiteDF$PTMPos <- as.numeric(str_sub(subsiteDF$PTM,2,30))

subsiteDF <- merge(subsiteDF,uniprot_fa_DF,by.x='UniprotID',by.y='UniprotID',all.x=T)
subsiteDF$Length <- sapply(subsiteDF$Sequence,str_length)
subsiteDF$Peptide <- paste(str_sub(subsiteDF$Sequence,ifelse((subsiteDF$PTMPos-7)<1,1,subsiteDF$PTMPos-7),subsiteDF$PTMPos-1),
                           str_sub(subsiteDF$Sequence,subsiteDF$PTMPos,subsiteDF$PTMPos),
                           str_sub(subsiteDF$Sequence,subsiteDF$PTMPos+1,ifelse((subsiteDF$PTMPos+7)>subsiteDF$Length,subsiteDF$Length,subsiteDF$PTMPos+7)),sep='')
subsiteDF <- subsiteDF[order(subsiteDF$Rho,decreasing=T),]
subsiteDF2 <- subsiteDF[subsiteDF$Rho>=exp_median,]
subsiteDF2$PeptideLen <- sapply(subsiteDF2$Peptide,str_length)
subsiteDF2 <- subsiteDF2[subsiteDF2$PeptideLen==15,]
subsiteDF2$Project <- 'HNSCC_T'


tmpdata <- filedata4[filedata4$Project=='HNSCC_N',]
# all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
# all_median
exp_median

subsiteDF <- aggregate(tmpdata[,c('Rho','sample.num')],list(tmpdata$GeneSite),max,na.rm=T)
colnames(subsiteDF)[1] <- 'GeneSite'

subsiteDF[,c('UniprotID','GeneName','PTM')] <- str_split_fixed(subsiteDF$GeneSite,'#',3)
subsiteDF$PTMtype <- str_sub(subsiteDF$PTM,1,1)
subsiteDF$PTMPos <- as.numeric(str_sub(subsiteDF$PTM,2,30))

subsiteDF <- merge(subsiteDF,uniprot_fa_DF,by.x='UniprotID',by.y='UniprotID',all.x=T)
subsiteDF$Length <- sapply(subsiteDF$Sequence,str_length)
subsiteDF$Peptide <- paste(str_sub(subsiteDF$Sequence,ifelse((subsiteDF$PTMPos-7)<1,1,subsiteDF$PTMPos-7),subsiteDF$PTMPos-1),
                           str_sub(subsiteDF$Sequence,subsiteDF$PTMPos,subsiteDF$PTMPos),
                           str_sub(subsiteDF$Sequence,subsiteDF$PTMPos+1,ifelse((subsiteDF$PTMPos+7)>subsiteDF$Length,subsiteDF$Length,subsiteDF$PTMPos+7)),sep='')
subsiteDF <- subsiteDF[order(subsiteDF$Rho,decreasing=T),]
subsiteDF3 <- subsiteDF[subsiteDF$Rho>=exp_median,]
subsiteDF3$PeptideLen <- sapply(subsiteDF3$Peptide,str_length)
subsiteDF3 <- subsiteDF3[subsiteDF3$PeptideLen==15,]
subsiteDF3$Project <- 'HNSCC_N'

subsiteDF4 <- rbind(subsiteDF2,subsiteDF3)
write.table(subsiteDF4[,c(1:9,15,17)],'/hwdata/home/wuyj2/eKPI/ArticleCancerStat/Kinase_rst/EEF2K_HBVHCC/CDC7_HNSCC_RNA.tsv',sep='\t',row.names = F,quote = F)

length(intersect(subsiteDF4[subsiteDF4$Project=='HNSCC_N','GeneSite'],subsiteDF4[subsiteDF4$Project=='HNSCC_T','GeneSite']))
data <- subsiteDF4
data$PeptideLen <- sapply(data$Peptide,str_length)
data <- data[data$PeptideLen==15,]

data_N <- data[data$Project=='HNSCC_N',c('GeneSite','Peptide')]
data_N[,paste('Pos',seq(1,15),sep='')] <- str_split_fixed(data_N$Peptide,'',15)
data_table_N <- data.frame(table(data_N[,'Pos1']))
data_table_N$Freq <- round(data_table_N$Freq/sum(data_table_N$Freq,na.rm = T),3)
colnames(data_table_N) <- c('ProName','Pos1')
for(nm in paste('Pos',seq(2,15),sep='')){
  tmp_data_table <- data.frame(table(data_N[,nm]))
  tmp_data_table$Freq <- round(tmp_data_table$Freq/sum(tmp_data_table$Freq,na.rm = T),3)
  colnames(tmp_data_table) <- c('ProName',nm)
  data_table_N <- merge(data_table_N,tmp_data_table,by='ProName',all=T)
}
# data_table_N$ProNum <- rowSums(data_table_N[,2:16],na.rm = T)
rownames(data_table_N) <- paste('N_',data_table_N$ProName,sep='')
data_table_N <- data_table_N[,2:16]
data_table_N[is.na(data_table_N)] <- 0
data_table_N <- as.matrix(data_table_N)

# data_table_N <- data_table_N[c('N_S','N_P','N_E','N_K','N_R'),4:12]*100
# anno_data_N <- data_table_N[c('N_S','N_P','N_E','N_K','N_R'),]

data_T <- data[data$Project=='HNSCC_T',c('GeneSite','Peptide')]
data_T[,paste('Pos',seq(1,15),sep='')] <- str_split_fixed(data_T$Peptide,'',15)
data_table_T <- data.frame(table(data_T[,'Pos1']))
data_table_T$Freq <- round(data_table_T$Freq/sum(data_table_T$Freq,na.rm = T),4)
colnames(data_table_T) <- c('ProName','Pos1')
for(nm in paste('Pos',seq(2,15),sep='')){
  tmp_data_table <- data.frame(table(data_T[,nm]))
  tmp_data_table$Freq <- round(tmp_data_table$Freq/sum(tmp_data_table$Freq,na.rm = T),3)
  colnames(tmp_data_table) <- c('ProName',nm)
  data_table_T <- merge(data_table_T,tmp_data_table,by='ProName',all=T)
}
rownames(data_table_T) <- paste('T_',data_table_T$ProName,sep='')
data_table_T <- data_table_T[,2:16]
data_table_T[is.na(data_table_T)] <- 0
data_table_T <- as.matrix(data_table_T)

ttestDF <- data.frame()
for(i in 1:nrow(data_table_T)){
  for(j in 1:nrow(data_table_N)){
    ttest <- t.test(data_table_T[i,3:13],data_table_N[i,3:13],paired = T)
    ttestDF <- rbind(ttestDF,data.frame(Tname = rownames(data_table_T)[i],
                                        Nname = rownames(data_table_N)[j],
                                        TMean = median(data_table_T[i,]),
                                        NMean = median(data_table_N[j,]),
                                        pvalue = ttest[['p.value']]
    ))
  }
}

##########FigureS3A BCKDK#########
fileDataFrame <- cancerFrame[cancerFrame$CancerName=='HBVHCC',]
filedata <- fread(paste(fileDataFrame[1,'fileName'],sep=''),
                  header = T,data.table = F)
filedata <- filedata[filedata$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$Rho.pvalue<=Pvalue_THRD&filedata$Rho>=Rvalue_THRD
filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
# filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
# colnames(filedata) <- c('ID','Rho')
filedata$expriment <- FALSE
filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
filedata$Project <- 'HBVHCC_N'

filedata2 <- fread(paste(fileDataFrame[4,'fileName'],sep=''),
                   header = T,data.table = F)
filedata2 <- filedata2[filedata2$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata2$Rho.pvalue<=Pvalue_THRD&filedata2$Rho>=Rvalue_THRD
filedata2$KinaseName <- str_split_fixed(filedata2$RNA.names,'_',2)[,1]
filedata2 <- filedata2[which(filedata2$KinaseName %in% kinaseNames),]#筛选激酶
filedata2$SubstrateName <- str_split_fixed(filedata2$GeneSite,'#',3)[,2]
filedata2 <- filedata2[filedata2$KinaseName!=filedata2$SubstrateName,]#过滤自磷酸化的KPS pairs
filedata2$ID  <- paste(filedata2$KinaseName,filedata2$GeneSite,sep='#')
# filedata2 <- aggregate(filedata2[,c('Rho')],list(filedata2$ID),max,na.rm=T)
# colnames(filedata2) <- c('ID','Rho')
filedata2$expriment <- FALSE
filedata2[which(filedata2$ID %in% exprimentDF$ID),]$expriment <- TRUE
filedata2$Project <- 'HBVHCC_T'

filedata3 <- rbind(filedata,filedata2)

filedata4 <- filedata3[filedata3$KinaseName=='EEF2K'&filedata3$Rho.pvalue<=Pvalue_THRD&filedata3$Rho>=Rvalue_THRD,]
tmpdata <- filedata4[filedata4$Project=='HBVHCC_T',]
all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
all_median
exp_median

wilcoxtest <- wilcox.test(tmpdata[tmpdata$expriment==TRUE,]$Rho,tmpdata[tmpdata$expriment==FALSE,]$Rho,alternative = 'greater')
pvalue <- format(wilcoxtest[["p.value"]], scientific = T,digits = 3)

p1 <- ggplot(tmpdata,aes(x=Rho,Group=expriment,fill=expriment))+
  geom_density(alpha=0.3)+
  theme_bw()+
  facet_wrap(~Project)+
  geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
  geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
  scale_x_continuous(name = "")+
  scale_y_continuous(name = "Density")+
  scale_fill_manual(values=c('#B54764','#9D9EA3'),
                    breaks = c(TRUE,FALSE),
                    labels = c('Exp-confirmed','All'))+
  annotate("text",label=paste('p = ',pvalue,sep=''),x=0.6,y=5.3,size=6)+
  theme_classic(base_size = 18)+
  theme(panel.grid=element_blank(),
        legend.box = element_blank(),
        legend.position= "inside",
        legend.position.inside = c(0.7,0.87),
        panel.border = element_rect(fill=NA,linewidth = 1),
        # text=element_text(size=14,face="plain",color="black"),
        # axis.title=element_text(size=16,face="plain",color="black"),
        # axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_blank(),
        # legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank())
p1

tmpdata <- filedata4[filedata4$Project=='HBVHCC_N',]
all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
all_median
exp_median
wilcoxtest <- wilcox.test(tmpdata[tmpdata$expriment==TRUE,]$Rho,tmpdata[tmpdata$expriment==FALSE,]$Rho,alternative = 'greater')
pvalue <- format(wilcoxtest[["p.value"]], scientific = T,digits = 3)

p2 <- ggplot(tmpdata,aes(x=Rho,Group=expriment,fill=expriment))+
  geom_density(alpha=0.3)+
  theme_bw()+
  facet_wrap(~Project)+
  geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
  geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
  scale_x_continuous(name = "Spearman' Rho")+
  scale_y_continuous(name = "Density")+
  scale_fill_manual(values=c('#B54764','#9D9EA3'),
                    breaks = c(TRUE,FALSE),
                    labels = c('Exp-confirmed','All'))+
  annotate("text",label=paste('p = ',pvalue,sep=''),x=0.65,y=5,size=6)+
  theme_classic(base_size = 18)+
  theme(panel.grid=element_blank(),
        legend.box = element_blank(),
        legend.position= "none",
        panel.border = element_rect(fill=NA,linewidth = 1),
        # text=element_text(size=14,face="plain",color="black"),
        # axis.title=element_text(size=16,face="plain",color="black"),
        # axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_blank(),
        # legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank())

# p <- p1 %>% insert_bottom(p2)

p <- plot_grid(plotlist = list(p2,NULL,p1),ncol=1,
          rel_heights = c(1,0.0000001,1),
          rel_widths = c(1))
p
ggsave(paste('/home/wuyj2/eKPI/ArticleCancerStat/Figure/Figure3B_EEF2K_HBVHCC_Phos_density','.pdf',sep=''),p,width =4,height = 6.5)


#BCKDK
filedata4 <- filedata3[filedata3$KinaseName=='BCKDK'&filedata3$Rho.pvalue<=Pvalue_THRD&filedata3$Rho>=Rvalue_THRD,]
tmpdata <- filedata4[filedata4$Project=='HBVHCC_T',]
all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
all_median
exp_median

wilcoxtest <- wilcox.test(tmpdata[tmpdata$expriment==TRUE,]$Rho,tmpdata[tmpdata$expriment==FALSE,]$Rho,alternative = 'greater')
pvalue <- format(wilcoxtest[["p.value"]], scientific = T,digits = 3)

p1 <- ggplot(tmpdata,aes(x=Rho,Group=expriment,fill=expriment))+
  geom_density(alpha=0.3)+
  theme_bw()+
  facet_wrap(~Project)+
  geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
  geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
  scale_x_continuous(name = "")+
  scale_y_continuous(name = "Density")+
  scale_fill_manual(values=c('#B54764','#9D9EA3'),
                    breaks = c(TRUE,FALSE),
                    labels = c('Exp-confirmed','All'))+
  annotate("text",label=paste('p = ',pvalue,sep=''),x=0.6,y=4.5,size=6)+
  theme_classic(base_size = 18)+
  theme(panel.grid=element_blank(),
        legend.box = element_blank(),
        legend.position= "inside",
        legend.position.inside = c(0.5,0.87),
        panel.border = element_rect(fill=NA,linewidth = 1),
        # text=element_text(size=14,face="plain",color="black"),
        # axis.title=element_text(size=16,face="plain",color="black"),
        # axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_blank(),
        # legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank())
p1

tmpdata <- filedata4[filedata4$Project=='HBVHCC_N',]
all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
all_median
exp_median
wilcoxtest <- wilcox.test(tmpdata[tmpdata$expriment==TRUE,]$Rho,tmpdata[tmpdata$expriment==FALSE,]$Rho,alternative = 'greater')
pvalue <- format(wilcoxtest[["p.value"]], scientific = T,digits = 3)

p2 <- ggplot(tmpdata,aes(x=Rho,Group=expriment,fill=expriment))+
  geom_density(alpha=0.3)+
  theme_bw()+
  facet_wrap(~Project)+
  geom_vline(xintercept = exp_median,linetype=2,color='#B54764') + # 添加base mean的水平线
  geom_vline(xintercept = all_median,linetype=2,color='#9D9EA3') + # 添加base mean的水平线
  scale_x_continuous(name = "Spearman' Rho")+
  scale_y_continuous(name = "Density")+
  scale_fill_manual(values=c('#B54764','#9D9EA3'),
                    breaks = c(TRUE,FALSE),
                    labels = c('Exp-confirmed','All'))+
  annotate("text",label=paste('p = ',pvalue,sep=''),x=0.58,y=6,size=6)+
  theme_classic(base_size = 18)+
  theme(panel.grid=element_blank(),
        legend.box = element_blank(),
        legend.position= "none",
        panel.border = element_rect(fill=NA,linewidth = 1),
        # text=element_text(size=14,face="plain",color="black"),
        # axis.title=element_text(size=16,face="plain",color="black"),
        # axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_blank(),
        # legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank())

# p <- p1 %>% insert_bottom(p2)

p <- plot_grid(plotlist = list(p2,NULL,p1),ncol=1,
               rel_heights = c(1,0.0000001,1),
               rel_widths = c(1))
p
ggsave(paste('/home/wuyj2/eKPI/ArticleCancerStat/Figure/Figure3B_BCKDK_HBVHCC_Phos_density','.pdf',sep=''),p,width =4,height = 8)



##########FigureS3B readseq BCKDK#############
# uniprot_fa_DF <- read.csv('/home/wuyj2/eKPI/CancerOmics/Uniprot_ID_sequence2.csv',header = T,row.names = 1)
# fileDataFrame <- cancerFrame[cancerFrame$CancerName=='HBVHCC',]
# filedata <- fread(paste(fileDataFrame[1,'fileName'],sep=''),
#                   header = T,data.table = F)
# filedata <- filedata[filedata$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata$Rho.pvalue<=Pvalue_THRD&filedata$Rho>=Rvalue_THRD
# filedata$KinaseName <- str_split_fixed(filedata$RNA.names,'_',2)[,1]
# filedata <- filedata[which(filedata$KinaseName %in% kinaseNames),]#筛选激酶
# filedata$SubstrateName <- str_split_fixed(filedata$GeneSite,'#',3)[,2]
# filedata <- filedata[filedata$KinaseName!=filedata$SubstrateName,]#过滤自磷酸化的KPS pairs
# filedata$ID  <- paste(filedata$KinaseName,filedata$GeneSite,sep='#')
# # filedata <- aggregate(filedata[,c('Rho')],list(filedata$ID),max,na.rm=T)
# # colnames(filedata) <- c('ID','Rho')
# filedata$expriment <- FALSE
# filedata[which(filedata$ID %in% exprimentDF$ID),]$expriment <- TRUE
# filedata$Project <- 'HBVHCC_N'
# 
# filedata2 <- fread(paste(fileDataFrame[4,'fileName'],sep=''),
#                    header = T,data.table = F)
# filedata2 <- filedata2[filedata2$sample.num>=Sample_THRD,]#根据样本数量、Rho、pvalue值过滤满足条件的KPS pair#&filedata2$Rho.pvalue<=Pvalue_THRD&filedata2$Rho>=Rvalue_THRD
# filedata2$KinaseName <- str_split_fixed(filedata2$RNA.names,'_',2)[,1]
# filedata2 <- filedata2[which(filedata2$KinaseName %in% kinaseNames),]#筛选激酶
# filedata2$SubstrateName <- str_split_fixed(filedata2$GeneSite,'#',3)[,2]
# filedata2 <- filedata2[filedata2$KinaseName!=filedata2$SubstrateName,]#过滤自磷酸化的KPS pairs
# filedata2$ID  <- paste(filedata2$KinaseName,filedata2$GeneSite,sep='#')
# # filedata2 <- aggregate(filedata2[,c('Rho')],list(filedata2$ID),max,na.rm=T)
# # colnames(filedata2) <- c('ID','Rho')
# filedata2$expriment <- FALSE
# filedata2[which(filedata2$ID %in% exprimentDF$ID),]$expriment <- TRUE
# filedata2$Project <- 'HBVHCC_T'
# 
# filedata3 <- rbind(filedata,filedata2)
# 
# filedata4 <- filedata3[filedata3$KinaseName=='EEF2K'&filedata3$Rho.pvalue<=Pvalue_THRD&filedata3$Rho>=Rvalue_THRD,]
# tmpdata <- filedata4[filedata4$Project=='HBVHCC_T',]
# # all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
# exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
# # all_median
# exp_median
# 
# subsiteDF <- aggregate(tmpdata[,c('Rho','sample.num')],list(tmpdata$GeneSite),max,na.rm=T)
# colnames(subsiteDF)[1] <- 'GeneSite'
# 
# subsiteDF[,c('UniprotID','GeneName','PTM')] <- str_split_fixed(subsiteDF$GeneSite,'#',3)
# subsiteDF$PTMtype <- str_sub(subsiteDF$PTM,1,1)
# subsiteDF$PTMPos <- as.numeric(str_sub(subsiteDF$PTM,2,30))
# 
# subsiteDF <- merge(subsiteDF,uniprot_fa_DF,by.x='UniprotID',by.y='UniprotID',all.x=T)
# subsiteDF$Length <- sapply(subsiteDF$Sequence,str_length)
# subsiteDF$Peptide <- paste(str_sub(subsiteDF$Sequence,ifelse((subsiteDF$PTMPos-7)<1,1,subsiteDF$PTMPos-7),subsiteDF$PTMPos-1),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos,subsiteDF$PTMPos),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos+1,ifelse((subsiteDF$PTMPos+7)>subsiteDF$Length,subsiteDF$Length,subsiteDF$PTMPos+7)),sep='')
# subsiteDF <- subsiteDF[order(subsiteDF$Rho,decreasing=T),]
# subsiteDF2 <- subsiteDF[subsiteDF$Rho>=exp_median,]
# subsiteDF2$PeptideLen <- sapply(subsiteDF2$Peptide,str_length)
# subsiteDF2 <- subsiteDF2[subsiteDF2$PeptideLen==15,]
# subsiteDF2$Project <- 'HBVHCC_T'
# 
# 
# tmpdata <- filedata4[filedata4$Project=='HBVHCC_N',]
# # all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
# exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
# # all_median
# exp_median
# 
# subsiteDF <- aggregate(tmpdata[,c('Rho','sample.num')],list(tmpdata$GeneSite),max,na.rm=T)
# colnames(subsiteDF)[1] <- 'GeneSite'
# 
# subsiteDF[,c('UniprotID','GeneName','PTM')] <- str_split_fixed(subsiteDF$GeneSite,'#',3)
# subsiteDF$PTMtype <- str_sub(subsiteDF$PTM,1,1)
# subsiteDF$PTMPos <- as.numeric(str_sub(subsiteDF$PTM,2,30))
# 
# subsiteDF <- merge(subsiteDF,uniprot_fa_DF,by.x='UniprotID',by.y='UniprotID',all.x=T)
# subsiteDF$Length <- sapply(subsiteDF$Sequence,str_length)
# subsiteDF$Peptide <- paste(str_sub(subsiteDF$Sequence,ifelse((subsiteDF$PTMPos-7)<1,1,subsiteDF$PTMPos-7),subsiteDF$PTMPos-1),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos,subsiteDF$PTMPos),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos+1,ifelse((subsiteDF$PTMPos+7)>subsiteDF$Length,subsiteDF$Length,subsiteDF$PTMPos+7)),sep='')
# subsiteDF <- subsiteDF[order(subsiteDF$Rho,decreasing=T),]
# subsiteDF3 <- subsiteDF[subsiteDF$Rho>=exp_median,]
# subsiteDF3$PeptideLen <- sapply(subsiteDF3$Peptide,str_length)
# subsiteDF3 <- subsiteDF3[subsiteDF3$PeptideLen==15,]
# subsiteDF3$Project <- 'HBVHCC_N'
# 
# subsiteDF4 <- rbind(subsiteDF2,subsiteDF3)
# write.table(subsiteDF4[,c(1:9,15:16)],'/home/wuyj2/eKPI/ArticleCancerStat/Kinase_rst/EEF2K_HBVHCC/EEF2K_HBVHCC_Phos.tsv',sep='\t',row.names = F,quote = F)
# 
# ###BCKDK
# filedata4 <- filedata3[filedata3$KinaseName=='BCKDK'&filedata3$Rho.pvalue<=Pvalue_THRD&filedata3$Rho>=Rvalue_THRD,]
# tmpdata <- filedata4[filedata4$Project=='HBVHCC_T',]
# # all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
# exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
# # all_median
# exp_median
# 
# subsiteDF <- aggregate(tmpdata[,c('Rho','sample.num')],list(tmpdata$GeneSite),max,na.rm=T)
# colnames(subsiteDF)[1] <- 'GeneSite'
# subsiteDF[,c('UniprotID','GeneName','PTM')] <- str_split_fixed(subsiteDF$GeneSite,'#',3)
# subsiteDF$PTMtype <- str_sub(subsiteDF$PTM,1,1)
# subsiteDF$PTMPos <- as.numeric(str_sub(subsiteDF$PTM,2,30))
# 
# subsiteDF <- merge(subsiteDF,uniprot_fa_DF,by.x='UniprotID',by.y='UniprotID',all.x=T)
# subsiteDF$Length <- sapply(subsiteDF$Sequence,str_length)
# subsiteDF$Peptide <- paste(str_sub(subsiteDF$Sequence,ifelse((subsiteDF$PTMPos-7)<1,1,subsiteDF$PTMPos-7),subsiteDF$PTMPos-1),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos,subsiteDF$PTMPos),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos+1,ifelse((subsiteDF$PTMPos+7)>subsiteDF$Length,subsiteDF$Length,subsiteDF$PTMPos+7)),sep='')
# subsiteDF <- subsiteDF[order(subsiteDF$Rho,decreasing=T),]
# # subsiteDF2 <- subsiteDF[subsiteDF$Rho>=exp_median,]
# subsiteDF2 <- subsiteDF[1:round(nrow(subsiteDF)/10,0),]
# subsiteDF2$PeptideLen <- sapply(subsiteDF2$Peptide,str_length)
# subsiteDF2 <- subsiteDF2[subsiteDF2$PeptideLen==15,]
# subsiteDF2$Project <- 'HBVHCC_T'
# 
# 
# tmpdata <- filedata4[filedata4$Project=='HBVHCC_N',]
# # all_median <- quantile(tmpdata[tmpdata$expriment==FALSE,]$Rho,0.5)
# exp_median <- quantile(tmpdata[tmpdata$expriment==TRUE,]$Rho,0.5)
# # all_median
# exp_median
# 
# subsiteDF <- aggregate(tmpdata[,c('Rho','sample.num')],list(tmpdata$GeneSite),max,na.rm=T)
# colnames(subsiteDF)[1] <- 'GeneSite'
# 
# subsiteDF[,c('UniprotID','GeneName','PTM')] <- str_split_fixed(subsiteDF$GeneSite,'#',3)
# subsiteDF$PTMtype <- str_sub(subsiteDF$PTM,1,1)
# subsiteDF$PTMPos <- as.numeric(str_sub(subsiteDF$PTM,2,30))
# 
# subsiteDF <- merge(subsiteDF,uniprot_fa_DF,by.x='UniprotID',by.y='UniprotID',all.x=T)
# subsiteDF$Length <- sapply(subsiteDF$Sequence,str_length)
# subsiteDF$Peptide <- paste(str_sub(subsiteDF$Sequence,ifelse((subsiteDF$PTMPos-7)<1,1,subsiteDF$PTMPos-7),subsiteDF$PTMPos-1),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos,subsiteDF$PTMPos),
#                            str_sub(subsiteDF$Sequence,subsiteDF$PTMPos+1,ifelse((subsiteDF$PTMPos+7)>subsiteDF$Length,subsiteDF$Length,subsiteDF$PTMPos+7)),sep='')
# subsiteDF <- subsiteDF[order(subsiteDF$Rho,decreasing=T),]
# # subsiteDF3 <- subsiteDF[subsiteDF$Rho>=exp_median,]
# subsiteDF3 <- subsiteDF[1:round(nrow(subsiteDF)/10,0),]
# subsiteDF3$PeptideLen <- sapply(subsiteDF3$Peptide,str_length)
# subsiteDF3 <- subsiteDF3[subsiteDF3$PeptideLen==15,]
# subsiteDF3$Project <- 'HBVHCC_N'
# 
# subsiteDF4 <- rbind(subsiteDF2,subsiteDF3)
# write.table(subsiteDF4[,c(1:9,15,17)],'/home/wuyj2/eKPI/ArticleCancerStat/Kinase_rst/EEF2K_HBVHCC/BCKDK_HBVHCC_Phos.tsv',sep='\t',row.names = F,quote = F)
# 
# length(intersect(subsiteDF4[subsiteDF4$Project=='HBVHCC_N','GeneSite'],subsiteDF4[subsiteDF4$Project=='HBVHCC_T','GeneSite']))
# data <- subsiteDF4
# data$PeptideLen <- sapply(data$Peptide,str_length)
# data <- data[data$PeptideLen==15,]
# 
# data_N <- data[data$Project=='HBVHCC_N',c('GeneSite','Peptide')]
# data_N[,paste('Pos',seq(1,15),sep='')] <- str_split_fixed(data_N$Peptide,'',15)
# data_table_N <- data.frame(table(data_N[,'Pos1']))
# data_table_N$Freq <- round(data_table_N$Freq/sum(data_table_N$Freq,na.rm = T),3)
# colnames(data_table_N) <- c('ProName','Pos1')
# for(nm in paste('Pos',seq(2,15),sep='')){
#   tmp_data_table <- data.frame(table(data_N[,nm]))
#   tmp_data_table$Freq <- round(tmp_data_table$Freq/sum(tmp_data_table$Freq,na.rm = T),3)
#   colnames(tmp_data_table) <- c('ProName',nm)
#   data_table_N <- merge(data_table_N,tmp_data_table,by='ProName',all=T)
# }
# # data_table_N$ProNum <- rowSums(data_table_N[,2:16],na.rm = T)
# rownames(data_table_N) <- paste('N_',data_table_N$ProName,sep='')
# data_table_N <- data_table_N[,2:16]
# data_table_N[is.na(data_table_N)] <- 0
# data_table_N <- as.matrix(data_table_N)
# 
# # data_table_N <- data_table_N[c('N_S','N_P','N_E','N_K','N_R'),4:12]*100
# # anno_data_N <- data_table_N[c('N_S','N_P','N_E','N_K','N_R'),]
# 
# data_T <- data[data$Project=='HBVHCC_T',c('GeneSite','Peptide')]
# data_T[,paste('Pos',seq(1,15),sep='')] <- str_split_fixed(data_T$Peptide,'',15)
# data_table_T <- data.frame(table(data_T[,'Pos1']))
# data_table_T$Freq <- round(data_table_T$Freq/sum(data_table_T$Freq,na.rm = T),4)
# colnames(data_table_T) <- c('ProName','Pos1')
# for(nm in paste('Pos',seq(2,15),sep='')){
#   tmp_data_table <- data.frame(table(data_T[,nm]))
#   tmp_data_table$Freq <- round(tmp_data_table$Freq/sum(tmp_data_table$Freq,na.rm = T),3)
#   colnames(tmp_data_table) <- c('ProName',nm)
#   data_table_T <- merge(data_table_T,tmp_data_table,by='ProName',all=T)
# }
# rownames(data_table_T) <- paste('T_',data_table_T$ProName,sep='')
# data_table_T <- data_table_T[,2:16]
# data_table_T[is.na(data_table_T)] <- 0
# data_table_T <- as.matrix(data_table_T)
# 
# ttestDF <- data.frame()
# for(i in 1:nrow(data_table_T)){
#   for(j in 1:nrow(data_table_N)){
#     ttest <- t.test(data_table_T[i,3:13],data_table_N[i,3:13],paired = T)
#     ttestDF <- rbind(ttestDF,data.frame(Tname = rownames(data_table_T)[i],
#                                         Nname = rownames(data_table_N)[j],
#                                         TMean = median(data_table_T[i,]),
#                                         NMean = median(data_table_N[j,]),           
#                                         pvalue = ttest[['p.value']]
#     ))
#   }
# }
# 