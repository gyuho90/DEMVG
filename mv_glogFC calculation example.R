######## mv-glogFC and p-value estimation (practical coding) #########
### input - RLE or PQN normalized and then VSN(glog) treated matrix , treatment group idx , control group idx
### output - mv-glogFC (and scale factor from DEMVG score) and p-value => make a DESeq2/apeglm/ashr like result table ###

### PPARa agonist liver data
bulk_raw_count<-read.csv("/home/gyuho/SRA_PPARA/countM/PPARA_rawdata.csv", row.names = 1)
bulk_meta_data<-read.csv("/home/gyuho/SRA_PPARA/countM/PPARA_metadata.csv" , row.names = 1)
rownames(bulk_raw_count)<-sub("\\.[0-9]+","",rownames(bulk_raw_count))  #remove dot num in ensembl_ID

bulk_meta_data$mixed<-rep(c("WT (Control)","WT (Wy14643)","WT (Fenofibrate)","KO (Control)"),each=3) 

### AAN
bulk_raw_count<-read.csv("/home/gyuho/kidney_AAN/countM/AAN_rawcount.csv" , row.names = 1)
bulk_meta_data<-data.frame("tissue"=rep("Kidney",10), "mixed"=c(rep("AAN",6),rep("Control",4)))
rownames(bulk_meta_data)<-colnames(bulk_raw_count)
rownames(bulk_raw_count)<-sub("\\.[0-9]+","",rownames(bulk_raw_count))

## HFpEF
###
heart_p53_rawcount<-read.csv("/home/gyuho/heart_p53/countM/heart_p53_rawcount.csv",row.names = 1)
rownames(heart_p53_rawcount)<-sub("\\.[0-9]+","",rownames(heart_p53_rawcount)) 
bulk_raw_count<-heart_p53_rawcount

bulk_meta_data$heart<-rep(rep(c("HFpEF","Control"),each=3),times=2)  
bulk_meta_data$sex<-rep(c("Female","Male"),each=6) 
bulk_meta_data$mixed<-bulk_meta_data$heart 
bulk_meta_data

### C57BL6J and DBA2J
bulk_raw_count<-read.csv("/home/gyuho/C57BL6J_DBA2J/countM/C57BL6J_DBA2J_rawcount.csv", row.names = 1)
rownames(bulk_raw_count)<-sub("\\.[0-9]+","",rownames(bulk_raw_count))  #remove dot num in ensembl_ID

bulk_meta_data<-read.csv("/home/gyuho/C57BL6J_DBA2J/countM/C57BL6J_DBA2J_metadata.csv" , row.names = 1)
bulk_meta_data
# two diff batch
bulk_meta_data$AvgSpotLen==76
bulk_meta_data$AvgSpotLen==70

#for AvgSpotLen==70
bulk_raw_count<-bulk_raw_count[,bulk_meta_data$AvgSpotLen==70]   #length 70(evaluation) or 76(validation)
bulk_meta_data<-bulk_meta_data[bulk_meta_data$AvgSpotLen==70,]
#for AvgSpotLen==76
bulk_raw_count<-bulk_raw_count[,bulk_meta_data$AvgSpotLen==76]   #length 70(evaluation) or 76(validation)
bulk_meta_data<-bulk_meta_data[bulk_meta_data$AvgSpotLen==76,]

########### Assume a column named mixed exists in the metadata and is located in the last column##
data.table(bulk_meta_data)[,.I,by=mixed][,lapply(.SD,list),by=mixed]
data.table(bulk_meta_data)[,.I,by=mixed]
bulk_meta_data_group_idx<-split(data.table(bulk_meta_data)[,.I,by=mixed]$I,data.table(bulk_meta_data)[,.I,by=mixed]$mixed)
bulk_meta_data_group_idx

bulk_meta_data_contrast_name<-names(bulk_meta_data_group_idx)[c(1,2)]  #treatment, control
bulk_meta_data_contrast_name

bulk_raw_count_t_idx <-bulk_meta_data_group_idx[[1]]
bulk_raw_count_c_idx<-bulk_meta_data_group_idx[[2]]
bulk_raw_count_t_idx
bulk_raw_count_c_idx

library(DESeq2)
library(ashr)
library(apeglm)
dds_bulkRNA <<- DESeqDataSetFromMatrix(
  countData =  bulk_raw_count,            
  colData = bulk_meta_data,               
  design =  ~mixed  )

dds_bulkRNA<<-DESeq(dds_bulkRNA)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing


vsd_bulkRNA <- DESeq2::varianceStabilizingTransformation(dds_bulkRNA, blind=TRUE)	
dds_bulkRNA_res<- results(dds_bulkRNA, contrast = c(colnames(bulk_meta_data)[ncol(bulk_meta_data)],bulk_meta_data_contrast_name),
          alpha = 0.05, # p-value
          lfcThreshold = 0) #see below.

dds_bulkRNA_res # log2 fold change (MLE): mixed DBA2J vs C57BL6J 
dds_bulkRNA_res_lfc_ashr<- lfcShrink(dds_bulkRNA,  contrast = c(colnames(bulk_meta_data)[ncol(bulk_meta_data)],bulk_meta_data_contrast_name),
            type = "ashr",
            res= dds_bulkRNA_res)
 dds_bulkRNA_res_lfc_ashr # log2 fold change (MMSE): mixed DBA2J vs C57BL6J 

dds_bulkRNA$mixed<-factor(dds_bulkRNA$mixed)
dds_bulkRNA$mixed<-relevel(dds_bulkRNA$mixed, ref = bulk_meta_data_contrast_name[2])## ref= the name of control condition !

dds_bulkRNA <- nbinomWaldTest(dds_bulkRNA)
dds_bulkRNA_res_lfc_apeglm<- lfcShrink(dds_bulkRNA,  coef = 2,type = "apeglm")     
dds_bulkRNA_res_lfc_apeglm # log2 fold change (MAP): mixed DBA2J vs C57BL6J 

##@## baseMean, pvalue/padj are the same while log2FC value are different! (after treating shrinkage)
##@## ensembl_gene_id orders are conserved!! (all the same) 
# dds_bulkRNA_res[[1]]@listData$log2FoldChange
# dds_bulkRNA_res_lfc[[1]]@listData$log2FoldChange

summary(dds_bulkRNA_res) 
dds_bulkRNA_res_lfc

# cf.
# bulk_raw_count_rownames_tag<-unique(grcm39_all_annot[,c(1,3)])[data.table(ensembl=rownames(bulk_raw_count)),on=.(ensembl_gene_id=ensembl)]
# bulk_raw_count_rownames_tag

### shrinked log2FC and original log2FC
dds_bulkRNA_table<-data.frame(dds_bulkRNA_res@listData,  ## add original log2FC (without considering baseMean exp)
                               "log2FoldChange_manual"=ntd_bulkRNA_log2FC,
                               "log2FoldChange_ashr"=dds_bulkRNA_res_lfc_ashr@listData$log2FoldChange, 
                                "log2FoldChange_apeglm"=dds_bulkRNA_res_lfc_apeglm@listData$log2FoldChange)
dds_bulkRNA_table[dds_bulkRNA_table$baseMean==0,"log2FoldChange"]<-0
dds_bulkRNA_table[dds_bulkRNA_table$baseMean==0,"log2FoldChange_apeglm"]<-0
dds_bulkRNA_table
quantile(dds_bulkRNA_table$log2FoldChange_manual/dds_bulkRNA_table$log2FoldChange,na.rm=T)
quantile(dds_bulkRNA_table$log2FoldChange_ashr/dds_bulkRNA_table$log2FoldChange,na.rm=T)
quantile(dds_bulkRNA_table$log2FoldChange_apeglm/dds_bulkRNA_table$log2FoldChange,na.rm=T)

vsd_bulkRNA_assay<<-assay(vsd_bulkRNA)
DEMVG_score_original<<-DEMVG::DEMVG(vsd_bulkRNA_assay,bulk_raw_count_t_idx, bulk_raw_count_c_idx)
DEMVG_score_original

##############################
# combn_matrix_treatment_idx<<-which(bulk_meta_data_sub$mixed== "DBA2J") 
# combn_matrix_control_idx<<-which(bulk_meta_data_sub$mixed=="C57BL6J")   
bulk_raw_count_t_idx
bulk_raw_count_c_idx
any(bulk_raw_count_t_idx%in% bulk_raw_count_c_idx) || any(bulk_raw_count_c_idx%in% bulk_raw_count_t_idx) # FALSE; mutually exclusive!

combn_matrix_t_and_c_length<- length(bulk_raw_count_t_idx)+length(bulk_raw_count_c_idx)
combn_matrix_t_and_c_length

#@# 12C6 = 924 ;  6C3 = 20, 10C4 = 210
combn_matrix_treatment<-combn(combn_matrix_t_and_c_length, length(bulk_raw_count_t_idx) ) # nCr [(treatment+control) sample size, treatment sample size]
combn_matrix_treatment

combn_matrix_control<-apply(combn_matrix_treatment,2,function(t){seq_len(combn_matrix_t_and_c_length)[! seq_len(combn_matrix_t_and_c_length)%in%t] })
combn_matrix_both<-rbind(combn_matrix_treatment,combn_matrix_control) #@# it works as a general index [background all combination cases]
combn_matrix_both

combn_matrix<-matrix(c(bulk_raw_count_t_idx,bulk_raw_count_c_idx)[combn_matrix_both],nrow=combn_matrix_t_and_c_length) ## put treatment/control vector
combn_matrix

#### limit ncol 3000 ####
#@# if the ncol of combn_matrix exceeds 3000, then make a subset with choose without replacement !
if( ncol(combn_matrix) > 3000){
combn_matrix_sub<<-combn_matrix[,sample(seq_len(ncol(combn_matrix)),3000,replace = FALSE)]  
}else{
combn_matrix_sub<<-combn_matrix
}

combn_matrix_t<-combn_matrix_sub[seq_len(length(bulk_raw_count_t_idx)),]
combn_matrix_c<-combn_matrix_sub[length(bulk_raw_count_t_idx)+seq_len(length(bulk_raw_count_c_idx)),] 

##### parallel computing.  only use core 10-14 d/t less matching memory space !
library(parallel)
detectCores()  # got the number of core in CPU
numCores <<- 10  # allocate the number of core
myCluster <<- parallel::makeCluster(numCores)

vsd_bulkRNA_assay<<-assay(vsd_bulkRNA)
parallel::clusterExport(myCluster, c("vsd_bulkRNA_assay","combn_matrix_t","combn_matrix_c")) #dataset will be used in parallel
clusterEvalQ(cl=myCluster, library(DEMVG))
DEMVG_score_list<<-parLapply(cl=myCluster,seq_len(ncol(combn_matrix_t)),function(t){
  DEMVG::DEMVG(vsd_bulkRNA_assay,combn_matrix_t[,t], combn_matrix_c[,t]) ##row[index for each treatment/control],  column for nth combi case 
})
length(DEMVG_score_list)

##
stopCluster(myCluster)  # Will free up the memory

###
DEMVG_score_matrix<<-data.matrix(data.frame(DEMVG_score_list))  
DEMVG_score_matrix[1:2,]

DEMVG_score_matrix_idx_delta_max<<-apply(DEMVG_score_matrix,1,function(t){ 
  sort_t=sort(t,decreasing=T)
  # sort_t[1] - sort_t[2]
  max(sort_t[1] - sort_t[2],sort_t[ncol(DEMVG_score_matrix)-1] - sort_t[ncol(DEMVG_score_matrix)])
  # abs diff between lowest and second lowest or highest and the second highest
})

DEMVG_score_matrix_others_mean<<-rowMeans(DEMVG_score_matrix)
DEMVG_score_matrix_others_sd<<-apply(DEMVG_score_matrix[,],1,sd) #sd calculated by /n-1

### calculate normal distribution and z-score / pvalue estimation
if( ncol(combn_matrix) > 3000){
  ## not all cases are permutated, so we need to calculate the original DEMVG score vector first
  DEMVG_score_original<<-DEMVG::DEMVG(vsd_bulkRNA_assay,bulk_raw_count_t_idx, bulk_raw_count_c_idx) 

}else{
  # the first case is the original DEMVG score vector
  DEMVG_score_original<<- DEMVG_score_matrix[,1]
}
DEMVG_score_z_values<-(DEMVG_score_original - DEMVG_score_matrix_others_mean)/DEMVG_score_matrix_others_sd 
DEMVG_score_z_values_pvalue<-pnorm(DEMVG_score_z_values, mean =0, sd = 1) #"Cumulative Distribution Function".

DEMVG_score_z_values_pvalue[DEMVG_score_z_values_pvalue>0.5] <- 1-DEMVG_score_z_values_pvalue[DEMVG_score_z_values_pvalue>0.5]  #one tail
DEMVG_score_z_values_pvalue2<<- 2*DEMVG_score_z_values_pvalue  #two tail
DEMVG_score_z_values_pvalue2

DEMVG_score_table<<-data.table("DEMVG_score"=DEMVG_score_original,"DEMVG_Zscore"=DEMVG_score_z_values,"DEMVG_pval"=DEMVG_score_z_values_pvalue2)
DEMVG_score_table[DEMVG_score_matrix_idx_delta_max <= median(DEMVG_score_matrix_idx_delta_max),DEMVG_pval:=NA] 
DEMVG_score_table$DEMVG_padj<-p.adjust(DEMVG_score_table$DEMVG_pval,method ="BH" ) #two tail pvalue treated with BH adjust

#@# mv_glogFC=DEMVG* median(log2FC/DEMVG)
DEMVG_score_table$mv_glogFC<- DEMVG_score_table$DEMVG * median(dds_bulkRNA_table$log2FoldChange/ DEMVG_score_table$DEMVG) ## log2FC from DESeq2

## merge DEMVG/DESeq2/apeglm/ashr results table ##
#bulk_raw_count_rownames_tag<-unique(grcm39_all_annot[,c(1,3)])[data.table(ensembl=rownames(bulk_raw_count)),on=.(ensembl_gene_id=ensembl)]
unique(grcm39_all_annot[,c(1,3)])[data.table(ensembl=rownames(bulk_raw_count)),on=.(ensembl_gene_id=ensembl)]
mv_glogFC_table_tag<-unique(grcm39_all_annot[,c(1,3)])[data.table(ensembl=rownames(bulk_raw_count)),on=.(ensembl_gene_id=ensembl)][data.table(dds_bulkRNA_table,DEMVG_score_table,keep.rownames = "ensembl_gene_id"),on="ensembl_gene_id"]
mv_glogFC_table_tag  #merging DESeq2, ashr, apeglm, DEMVG,mv_glogFC obj and form in a table !

