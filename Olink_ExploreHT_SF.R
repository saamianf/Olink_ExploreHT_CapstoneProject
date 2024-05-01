## READ IN ALL LIBRARIES FOR CODE
library(OlinkAnalyze)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

start <- Sys.time()
print('start Part 1 (LOD Evaluation) at')
print(start)

## ADJUST FILE NAMES FOR INPUTS + OUTPUTS THROUGHOUT CODE IF NEEDED
ht <- read.csv('Project_HT.csv')

## DETERMINE LOD CUTOFF VALUE FOR EACH ASSAY 
ht$Max_LOD <- NA
for (i in seq(1, dim(ht)[1], 96)){
  neg <- subset(ht[i:(i+95),], SampleType == 'NEGATIVE_CONTROL')
  neg_npx <- neg[,14]
  neg_ct <- neg[,12]
  LOD1 = median(neg_npx) + 3*((sd(neg_npx)*((length(neg_npx)-1)/length(neg_npx))))
  LOD2 = median(neg_npx) + 0.2
  LOD <- ifelse (LOD1 > LOD2, LOD1, LOD2)
  ht[i:(i+95),]$Max_LOD <- LOD
}

## DETERMINE WHETHER NPX OF SAMPLE FOR ASSAY FALLS BELOW THE MAXLOD
ht$LOD_check <- NA
for (i in 289:dim(ht)[1]){ ## start at 289 so it doesn't count the first few control things in the counts
  val <- ht[i,14]
  lod <- ht[i,20]
  if (ht[i,2] == 'SAMPLE') {
    ht[i,21] <- ifelse(val == 'NaN', NA, ifelse (val > lod, 'PASS', 'WARNING'))
  } else { ht[i,21] <- NA
  }}

## EVALUATE THE RATE OF WARNINGS BASED ON ASSAY
protein_check <- data.frame(with(ht,table(Assay, LOD_check)))
protein_check <- protein_check[order(protein_check$Assay),]
protein_check$Percent_Warning <- NA
for (i in seq(1, dim(protein_check)[1], 2)){
  p_val <- ifelse(protein_check[i,2] == 'PASS', protein_check[i,3], protein_check[i+1,3])
  w_val <- ifelse(p_val == protein_check[i,3], protein_check[i+1,3], protein_check[i,3])
  protein_check[i+1,4] <- paste(round((w_val / (p_val + w_val))*100, digits=3), '%', sep='')
}

## EVALUATE THE RATE OF WARNINGS BASED ON SAMPLE
sample_check <- data.frame(with(ht,table(SampleID, LOD_check)))
sample_check <- sample_check[order(sample_check$SampleID),]
sample_check$Percent_Warning <- NA
for (i in seq(1, dim(sample_check)[1], 2)){
  p_val <- ifelse(sample_check[i,2] == 'PASS', sample_check[i,3], sample_check[i+1,3])
  w_val <- ifelse(p_val == sample_check[i,3], sample_check[i+1,3], sample_check[i,3])
  sample_check[i+1,4] <- paste(round((w_val / (p_val + w_val))*100, digits=3), '%', sep='')
}

## EVALUATE THE RATE OF WARNINGS BASED ON BLOCK
block_check <- data.frame(with(ht,table(Block, LOD_check)))
block_check <- block_check[order(block_check$Block),]
block_check$Percent_Warning <- NA
for (i in seq(1, dim(block_check)[1], 2)){
  p_val <- ifelse(block_check[i,2] == 'PASS', block_check[i,3], block_check[i+1,3])
  w_val <- ifelse(p_val == block_check[i,3], block_check[i+1,3], block_check[i,3])
  block_check[i+1,4] <- paste(round((w_val / (p_val + w_val))*100, digits=3), '%', sep='')
}

## CREATE A WARNING SUBSETTED DATAFRAME
warnings <- subset(ht, LOD_check == 'WARNING' & AssayType == 'assay')

write.csv(ht, 'HT_LODcheck_SF.csv')
end <- Sys.time()
print('end Part 1 (LOD Evaluation) at')
print(end)

################################################################################

start <- Sys.time()
print('start Part 2 (General Figures) at')
print(start)

## READ IN HT FILE & ASSOCIATED MANIFEST
ht_npx <- read_NPX('Project_HT_LODcheck_SF.csv')
mani <- read_excel('Project_HT_manifest.xlsx')
ht_m <- merge(ht_npx, mani)

## FILTERING FOR PEDIATRIC ONLY SUBJECTS FROM FULL HT DATASET
ht_m_ped_all <- subset(ht_m, Subject_Group == 'P')
ht_m_ped <- subset(ht_m, Subject_Group == 'P' & Group =='PB') 
ht_m_n <- subset(ht_m, is.na(Subject_Group) == TRUE)
ht_m_ped_all <- rbind(ht_m_ped_all, ht_m_n)
ht_m_p <- rbind(ht_m_ped, ht_m_n)

## GENERATE A SCALED & UNSCALED HEATMAP 
heatmap1 <- olink_heatmap_plot(df=ht_m_ped_all, variable_row_list = 'Subject', cluster_rows=FALSE, 
                cluster_cols=FALSE, colnames='assay', show_colnames = FALSE,
                main = 'Heatmap (centered/scaled)')
#ggsave(path = 'figures', filename='heatmap_scaled.tiff', width = 10, height = 10)
heatmap2 <- olink_heatmap_plot(df=ht_m_ped_all, variable_row_list = 'Subject', center_scale=FALSE, 
                cluster_rows=FALSE, cluster_cols=FALSE, show_colnames = FALSE, 
                colnames='assay', main = 'Heatmap (uncentered/unscaled)')
#ggsave(path = 'figures', filename='heatmap_unscaled.tiff', width = 10, height = 10)

## GENERATE NPX DISTRIBUTION PLOTS BASED ON DIFFERENT QC CHECKS (Select the most applicable)
dist_assay <- olink_dist_plot(ht_m_ped_all, color_g = AssayQC) + 
  labs(title = 'AssayQC Distribution') + theme(axis.text.x = element_text(angle = 90))
#ggsave(path = 'figures', filename='dist_assay.tiff', width = 15, height = 4)
dist_sample <- olink_dist_plot(ht_m_ped_all, color_g = SampleQC) + 
  labs(title = 'SampleQC Distribution') + theme(axis.text.x = element_text(angle = 90))
#ggsave(path = 'figures', filename='dist_sample.png', width = 15, height = 4)
dist_lod <- olink_dist_plot(ht_m_ped_all, color_g = LOD_check) + 
  labs(title = 'LOD Distribution') + theme(axis.text.x = element_text(angle = 90))
#ggsave(path = 'figures', filename='dist_lod.tiff', width = 15, height = 4)

## GENERATE QC PLOTS BASED ON DIFFERENT QC CHECKS (Select the most applicable)
qc_assay <- olink_qc_plot(ht_m_ped_all, color_g = AssayQC) + labs(title = 'AssayQC QC Plot')
#ggsave(path = 'figures', filename='qc_assay.tiff', width = 8, height = 8)
qc_sample <- olink_qc_plot(ht_m_ped_all, color_g = SampleQC) +  labs(title = 'SampleQC QC Plot')
#ggsave(path = 'figures', filename='qc_sample.tiff', width = 8, height = 8)
qc_lod <- olink_qc_plot(ht_m_ped_all, color_g = LOD_check) + labs(title = 'LOD QC Plot')
#ggsave(path = 'figures', filename='qc_lod.tiff', width = 8, height = 8)

## GENERATE UMAP PLOTS WITH DIFFERENT GROUPINGS (Select the most applicable) 
umap_subject <- olink_umap_plot(ht_m_ped_all, color_g='Subject', quiet=TRUE)
#ggsave(path = 'figures', filename='umap_subject.tiff', width = 10, height = 10)
umap_p <- olink_umap_plot(ht_m_ped_all, color_g='Pre_Post', quiet=TRUE)
#ggsave(path = 'figures', filename='umap_p.tiff', width = 10, height = 10)
umap_g <- olink_umap_plot(ht_m_ped_all, color_g='Group', quiet=TRUE)
#ggsave(path = 'figures', filename='umap_g.tiff', width = 10, height = 10)

## RUN T-TESTS OF SUBGROUPS & GENERATE VOLCANO PLOTS
t_ped <- olink_ttest(ht_m_p, variable='Pre_Post')
vol_ped <- olink_volcano_plot(t_ped, x_lab = 'Log2FoldChange (NPX Difference)') + labs(title = 'Pediatric Samples')
#ggsave(path = 'figures', filename='ttest.png', width = 10, height = 10)

## GENERATE LOD-COLORED VOLCANO PLOT
protein_lods <- subset(protein_check, LOD_check == 'WARNING')
protein_lods$Percent_Warning <- gsub('%', '', protein_lods$Percent_Warning)
t_ped$LOD_asPercentages <- NA
for (i in 1:dim(t_ped)[1]){
 lod_c <- subset(protein_lods, Assay == as.character(t_ped[i,1]))
 t_ped[i,17] <- ifelse (as.numeric(lod_c[1,4]) < 70, ifelse(as.numeric(lod_c[1,4]) == 0, 'NONE (0)' , 'PASS (<70%)'), 'WARNING (>70%)')
}
cols <- c('NONE (0)' = 'green', 'PASS (<70%)' = 'black', 'WARNING (>70%)' = 'red')
vol_ped_lod <- olink_volcano_plot(t_ped, x_lab = 'Log2FoldChange (NPX Difference)') +
 geom_point(t_ped, mapping = aes(x=estimate, y=(-log10(p.value)), color = LOD_asPercentages)) +
 scale_color_manual(values=cols) + labs(title = 'All Samples - LOD labeled') +
 labs(color=('LOD_warnings'))
#ggsave(path = 'figures', filename='ttest_lod.png', width = 10, height = 10)

## RUN T-TESTS OF EACH PRE-POST PAIR & GENERATE VOLCANO PLOTS
t_ped_times_0 <- subset(ht_m_p, Time == '0' & Pre_Post == 'pre')
t_v <- function(x){
 this_set <- subset(ht_m_p, Time == x & Pre_Post == 'post')
 this_ttest <- olink_ttest(rbind(t_ped_times_0, this_set), variable='Pre_Post')
 this_vol <- olink_volcano_plot(this_ttest, x_lab='log2FoldChange (NPX Difference)') + labs(title = paste('Pediatric Samples - Day 0 pre / Day ', x, ' Pairs', sep=''))
 #ggsave(path = 'figures', filename=paste0('this_vol_0_', as.character(x), '.png', sep=''), width = 10, height = 10)
}
t_v(1)
t_v(4)
t_v(7)
t_v(10)
t_v(14)
t_v(21)
t_v(28)

end <- Sys.time()
print('end Part 2 (General Figures) at')
print(end)

################################################################################

start <- Sys.time()
print('start Part 3 (LOD Figures) at')
print(start)

assays <- unique(ht_m_p$Assay)

## GENERATE SUMMARY GRAPHS OF LOD DISTRIBUTION AMONG SAMPLES FOR EACH ASSAY & OUTPUT TO FOLDER (5440 total)
for (i in 1:length(assays)){
  assay_group <- na.omit(subset(ht_m_p, Assay == assays[i]))
  colls <- c('WARNING' = 'red', 'PASS' = 'black')
  graph <- ggplot(assay_group, aes(x=SampleID, y=NPX)) + geom_point(aes(color=LOD_check), size=2) + 
    theme_bw() + theme(axis.text.x = element_text(angle=90)) +
    geom_hline(yintercept = assay_group$Max_LOD, color='red', linetype='dashed') + 
    scale_color_manual(values = colls) +
    ggtitle(paste(assays[i], ' NPXs', sep=''))
  named <- paste(assays[i], '_OlinkHT_LODcheck.tiff', sep='')
  #ggsave(path = 'OlinkHT_LODcheck_ped_figures', filename=named, width = 8, height = 5)
}
end <- Sys.time()
print('end Part 3 (LOD Figures) at')
print(end)

################################################################################

start <- Sys.time()
print('start Part 4 (Olink Cross-Panel Comparisons) at')
print(start)

## READ IN & SORT TARGET96 DATA & MANIFEST
t96_npx <- read_excel('Project_Target96_NPX.xlsx')
mani96 <- read_excel('Project_Target96_manifest.xlsx') 
mani_re <- mani96[,4:9]
mani_re <- mani_re |> uncount(288) 
t96_wmani <- bind_cols(t96_npx, mani_re)
t96_data <- t96_wmani %>% arrange(Subject_Group, Subject_ID, Time)

## READ IN & SORT FLEX48 DATA & MANIFEST
f48_npx <- read_excel('Project_Flex48_NPX.xlsx')
mani48 <- read_excel('Project_Flex48_manifest.xlsx') 
f48_data <- merge(f48_npx, mani48)

## FILTERING FOR PEDIATRIC ONLY SUBJECTS FROM FULL TARGET96 & FLEX48 DATASET
t96_m_ped <- subset(t96_data, Subject_Group == 'P' & Group == 'PB') 
t96_m_n <- subset(t96_data, is.na(Subject_Group) == TRUE)
t96_m_p <- rbind(t96_m_ped, t96_m_n)
f48_m_ped <- subset(f48_data, Subject_Group == 'P' & Group == 'PB') 
f48_m_n <- subset(f48_data, is.na(Subject_Group) == TRUE)
f48_m_p <- rbind(f48_m_ped, f48_m_n)

## MAKE FULL T96 + F48 DATATABLE
tf_data <- full_join(t96_data, f48_data) 
tf_data$LOD_check <- ifelse(tf_data[,12] > tf_data[,11], 'PASS', 'WARNING')

## GET ASSAY + SUBJECT DATA
proteins <- read.csv('Project_Olink_Panels_Proteins.csv')
subjects <- unique(t96_m_p$Subject[!is.na(tf_data$Subject)])

collls <- c('Explore_HT' = 'black', 'Olink Target 96 Immuno-Oncology' = 'green', 'Olink Target 96 Oncology III' = 'blue', 'Olink Target 96 Immune Response' = 'purple', 'Flex Panel (FJPL-NVCM)' = 'red', 'Flex Panel (FTDS-ZZQZ)' = 'orange')

for (i in 1:dim(proteins)[1]){
  this_assay <- proteins[i,1]
  ht_subset <- na.omit(subset(ht_m_p, Assay == as.character(this_assay)))
  tf_subset <- subset(tf_data, Assay == as.character(this_assay))
  gene_set <- full_join(ht_subset, tf_subset)
  
  ## CALCULATE ADJUSTED NPX 
  baselines <- subset(gene_set, Time == 0 & Pre_Post == 'pre')
  gene_set$AdjustedNPX <- gene_set$NPX
  for (i in 1:dim(gene_set)[1]){
    baselineval <- baselines[which(baselines$Subject == gene_set[i,24] & baselines$Panel == gene_set[i,11] & baselines$Assay == gene_set[i, 9]),]
    if (dim(baselineval)[1] == 0) {
      gene_set[i,41] <- as.numeric(gene_set[i,15])
    } else {
      gene_set[i,41] <- as.numeric(gene_set[i,15]) - as.numeric(baselineval[1,15])
    }}
  
  for (j in 1:length(subjects)){
    t_subject <- subjects[j]
    sub_set <- subset(gene_set, Subject == as.character(t_subject))
    
    ## GENERATE LINE GRAPHS PER ASSAY + PER SUBJECT COMPARING HT + TARGET96 + FLEX48 DATASETS
    line_plot <- ggplot() + 
      geom_point(subset(sub_set, Time > -1), mapping=aes(x=Time, y=AdjustedNPX, color=Panel), size=2) + geom_line(subset(sub_set, Time > -1), mapping=aes(x=Time, y=AdjustedNPX, color=Panel), linetype='longdash') +
      theme_bw() +  xlab('Time (days)') + scale_color_manual(values = collls) + xlim(c(0, 28)) +
      ggtitle(paste(this_assay, ' for Subject ', t_subject, ' : Olink Panels', sep=''))
    named <- paste(this_assay, '-', t_subject, '.tiff', sep='')
    #ggsave(path = 'Olink_Lines', filename=named, height=8, width = 8)
  }

  ## HISTOGRAMS FOR FULL ASSAY
  hist_plot <- ggplot(subset(gene_set, Time > -1), aes(AdjustedNPX, fill=Panel)) +
    geom_histogram(alpha=0.5, position='dodge2', bins=50) + xlab('Values (AdjustedNPX)') + ylab('Count') +
    scale_fill_manual(values = collls) + 
    ggtitle(paste('Olink Value (AdjustedNPX) Distribution for ', this_assay, sep=''))
  thisname <- paste(this_assay, '_hist.tiff', sep='')
  #ggsave(path = 'Olink_Histograms', filename = thisname, width=8, height=8)
    
  ## BAR GRAPHS
  bar_plot <- ggplot(subset(gene_set, Time > -1), aes(factor(Time), AdjustedNPX, fill=Panel, group=Subject)) + 
    geom_bar(stat='identity', position = 'dodge2') + theme_bw() + 
    scale_fill_manual(values=collls) + 
    xlab('Time (days)') + ylab ('Values(AdjustedNPX)') + 
    ggtitle(paste('Olink Values (AdjustedNPX) for ', this_assay, sep=''))
  #ggsave(path = 'Olink_Bars', filename = paste(this_assay, '_bars.png', sep=''), width=8, height=8)
}

end <- Sys.time()
print('end Part 4 (Olink Cross-Panel Comparisons) at')
print(end)
