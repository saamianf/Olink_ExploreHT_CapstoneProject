## READ IN ALL LIBRARIES FOR CODE
library(OlinkAnalyze)
library(readxl)
library(ggplot2)

start <- Sys.time()
print('start Part 1 (LOD) at')
print(start)

## ADJUST FILE NAMES FOR INPUTS + OUTPUTS THROUGHOUT CODE IF NEEDED
ht <- read.csv('Project_HT.csv')

## DETERMINE LOD VALUE FOR EACH ASSAY 
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
print('end Part 1 (LOD) at')
print(end)

################################################################################

start <- Sys.time()
print('start Part 2 (General Figures) at')
print(start)

## FUNCTION TO OUTPUT FIGURES IN NEW WINDOW - REMOVE THE # FROM INSTANCES IN CODE TO HAVE OUTPUTS
## This is done instead of ggsave in order to allow flexibility with sizing of the outputs for each run
p <- function(r){
  dev.new()
  print(r) }

## READ IN HT FILE & ASSOCIATED MANIFEST
ht_npx <- read_NPX('HT_LODcheck_SF.csv')
mani <- read_excel('Project_Manifest.xlsx')
ht_m <- merge(ht_npx, mani)

## FILTERING FOR PEDIATRIC ONLY SUBJECTS FROM FULL HT DATASET
ht_m_ped <- subset(ht_m, Subject_Group == 'P') 
ht_m_n <- subset(ht_m, is.na(Subject_Group) == TRUE)
ht_m_p <- rbind(ht_m_ped, ht_m_n)

## GENERATE A SCALED & UNSCALED HEATMAP
heatmap1 <- olink_heatmap_plot(df=ht_m_p, variable_row_list = 'Subject', cluster_rows=FALSE, 
                cluster_cols=FALSE, colnames='assay', show_colnames = FALSE,
                main = 'Heatmap (centered/scaled)')
heatmap2 <- olink_heatmap_plot(df=ht_m_p, variable_row_list = 'Subject', center_scale=FALSE, 
                cluster_rows=FALSE, cluster_cols=FALSE, show_colnames = FALSE, 
                colnames='assay', main = 'Heatmap (uncentered/unscaled)')

## GENERATE PCA PLOTS W DIFFERENT GROUPINGS (Select the most applicable) 
pca_subject <- olink_pca_plot(ht_m_ped, color_g='Subject', quiet=TRUE)
pca_p <- olink_pca_plot(ht_m_ped, color_g='Pre_Post', quiet=TRUE)
pca_g <- olink_pca_plot(ht_m_ped, color_g='Group', quiet=TRUE)
pca_subject <- olink_pca_plot(subset(ht_m_ped, Group == 'PB'), color_g='Subject', quiet=TRUE)

## GENERATE NPX DISTRIBUTION PLOTS BASED ON DIFFERENT QC CHECKS (Select the most applicable)
dist_assay <- olink_dist_plot(ht_m_p, color_g = AssayQC) + 
  labs(title = 'AssayQC Distribution') + theme(axis.text.x = element_text(angle = 90))
dist_sample <- olink_dist_plot(ht_m_p, color_g = SampleQC) + 
  labs(title = 'SampleQC Distribution') + theme(axis.text.x = element_text(angle = 90))
dist_lod <- olink_dist_plot(ht_m_p, color_g = LOD_check) + 
  labs(title = 'LOD Distribution') + theme(axis.text.x = element_text(angle = 90))

## GENERATE QC PLOTS BASED ON DIFFERENT QC CHECKS (Select the most applicable)
qc_assay <- olink_qc_plot(ht_m_p, color_g = AssayQC) + 
  labs(title = 'AssayQC QC Plot') + theme(axis.text.x = element_text(angle = 90))
qc_sample <- olink_qc_plot(ht_m_p, color_g = SampleQC) +  
  labs(title = 'SampleQC QC Plot') + theme(axis.text.x = element_text(angle = 90))
qc_lod <- olink_qc_plot(ht_m_p, color_g = LOD_check) +  
  labs(title = 'LOD QC Plot') + theme(axis.text.x = element_text(angle = 90))

## GENERATE UMAP PLOTS WITH DIFFERENT GROUPINGS (Select the most applicable) 
umap_subject <- olink_umap_plot(ht_m_ped, color_g='Subject', quiet=TRUE)
umap_p <- olink_umap_plot(ht_m_ped, color_g='Pre_Post', quiet=TRUE)
umap_g <- olink_umap_plot(ht_m_ped, color_g='Group', quiet=TRUE)

#p(heatmap1)
#p(heatmap2)
#p(pca_subject)
#p(pca_p)
#p(pca_g)
#p(dist_assay)
#p(dist_sample)
#p(qc_assay)
#p(qc_sample)
#p(umap_subject)
#p(umap_p)
#p(umap_g)

## RUN T-TESTS OF SUBGROUPS & GENERATE VOLCANO PLOTS 
t_ped <- olink_ttest(ht_m_p, variable='Pre_Post') 
vol_ped <- olink_volcano_plot(t_ped, x_lab = 'Log2FoldChange (NPX Difference)') + labs(title = 'Pediatric Samples') 
t_ped_pb <- olink_ttest(subset(ht_m_p, Group == 'PB'), variable='Pre_Post')
vol_ped_pb <- olink_volcano_plot(t_ped, x_lab = 'Log2FoldChange (NPX Difference)') + labs(title = 'Pediatric Samples (PB only)')

#p(vol_ped)
#p(vol_ped_pb)

## GENERATE LOD-COLORED VOLCANO PLOT 
protein_lods <- subset(protein_check, LOD_check == 'WARNING')
protein_lods$Percent_Warning <- gsub('%', '', protein_lods$Percent_Warning)
t_all$LOD_asPercentages <- NA
for (i in 1:dim(t_ped)[1]){
  lod_c <- subset(protein_lods, Assay == as.character(t_ped[i,1]))
  t_ped[i,17] <- ifelse (as.numeric(lod_c[1,5]) < 50, ifelse(as.numeric(lod_c[1,5]) == 0, 'NONE (0)' , 'PASS (<50%)'), 'WARNING (>50%)')
}
cols <- c('NONE (0)' = 'green', 'PASS (<50%)' = 'black', 'WARNING (>50%)' = 'red')
vol_ped_lod <- olink_volcano_plot(t_ped, x_lab = 'Log2FoldChange (NPX Difference)') + 
  geom_point(t_ped, mapping = aes(x=estimate, y=(-log10(p.value)), color = LOD_asPercentages)) + 
  scale_color_manual(values=cols) + labs(title = 'All Samples - LOD labeled') +
  labs(color=('LOD_warnings'))

#p(vol_ped_lod)

## RUN T-TESTS OF EACH PRE-POST PAIR & GENERATE VOLCANO PLOTS
t_ped_times_0 <- subset(ht_m_p, Group == 'PB' & Time == '0' & Pre_Post == 'pre') 
t_v <- function(x){
  this_set <- subset(ht_m_p, Group == 'PB' & Time == x & Pre_Post == 'post')
  this_ttest <- olink_ttest(rbind(t_ped_times_0, this_set), variable='Pre_Post')
  assign(paste0('this_ttest_0_', as.character(x), sep=''), this_ttest)
  this_vol <- olink_volcano_plot(this_ttest, x_lab='log2FoldChange (NPX Difference)') + labs(title = paste('Pediatric Samples - Day 0 pre / Day ', x, ' Pairs', sep=''))
  assign(paste0('this_vol_0_', as.character(x), sep=''), this_vol)
  #p(this_vol)
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

assays <- unique(ht_m_p$Assay) ## use ht_m_p generated from part 1 --> can probably use the general graph for this + just remove sampleIDs --> which others could I do this for?? 

## GENERATE SUMMARY GRAPHS OF LOD DISTRIBUTION AMONG SAMPLES FOR EACH ASSAY & OUTPUT TO FOLDER (5440 total)
for (i in 1:length(assays)){
  assay_group <- na.omit(subset(ht_m_p, Assay == assays[i]))
  cols <- c('WARNING' = 'red', 'PASS' = 'black')
  graph <- ggplot(assay_group, aes(x=SampleID, y=NPX)) + geom_point(aes(color=LOD_check), size=2) + 
    theme_bw() + theme(axis.text.x = element_text(angle=90)) +
    geom_hline(yintercept = assay_group$Max_LOD, color='red', linetype='dashed') + 
    scale_color_manual(values = cols) +
    ggtitle(paste(assays[i], ' NPXs', sep=''))
  named <- paste(assays[i], '_OlinkHT_LODcheck.png', sep='')
  #ggsave(path = 'OlinkHT_LODcheck_ped_figures', filename=named)
}
end <- Sys.time()
print('end Part 3 (LOD Figures) at')
print(end)
################################################################################
start <- Sys.time()
print('start Part 4 (ExploreHT vs Target96 Figures) at')
print(start)

## READ IN & SORT TARGET96 DATA & MANIFEST
t96_npx <- read_excel('20240208_CAR123_Target96_NPX.xlsx') # #rename
mani <- read_excel('20240208_CAR123_Target96_manifest.xlsx') # rename
reorder = c(1,9,17,25,33,41,49,57,65,73,81,89,2,10,18,26,34,42,50,58,66,74,82,90,3,11,
            19,27,35,43,51,59,67,75,83,91,4,12,20,28,36,44,52,60,68,76,84,92,5,13,21,
            29,37,45,53,61,69,77,85,93,6,14,22,30,38,46,54,62,70,78,86,94,7,15,23,31,
            39,47,55,63,71,79,87,95,8,16,24,32,40,48,56,64,72,80,88,96) 
mani_re <- bind_cols(mani,reorder)
mani_re <- mani_re %>% arrange(reorder) 
mani_re <- mani_re[,4:9]
mani_re <- mani_re |> uncount(96) 
mani_rep <- rbind(mani_re, mani_re, mani_re)
t96_wmani <- bind_cols(t96_npx, mani_rep)
t96_data <- t96_wmani %>% arrange(Subject_Group, Subject_ID, Time)

## FILTERING FOR PEDIATRIC ONLY SUBJECTS FROM FULL TARGET96 DATASET
t96_m_ped <- subset(t96_data, Subject_Group == 'P') 
t96_m_n <- subset(t96_data, is.na(Subject_Group) == TRUE)
t96_m_p <- rbind(t96_m_ped, t96_m_n)

## GET PROTEIN + SAMPLE SUBJECT DATA
proteins <- read_excel('ExploreHT_Target96_OverlappingFull.xlsx')
subjects <- unique(t96_m_p$Subject[!is.na(t96_m_p$Subject)])

## GENERATE GRAPHS PER PROTEIN + PER SUBJECT COMPARING HT + TARGET96 DATASETS
for (i in 1:dim(proteins)[1]){
  this_assay <- proteins[i,1]
  ht_subset <- na.omit(subset(ht_m_p, Assay == as.character(this_assay)))
  t96_subset <- na.omit(subset(t96_m_p, Assay == as.character(this_assay)))
  for (j in 1:length(subjects)){
    t_subject <- subjects[j]
    ht_subset2 <- subset(ht_subset, Subject == as.character(t_subject))
    t96_subset2 <- subset(t96_subset, Subject == as.character(t_subject))
    full_set <- full_join(ht_subset2, t96_subset2)
    cols <- c('Explore_HT' = 'black', 'Olink Target 96 Immuno-Oncology' = 'green', 'Olink Target 96 Oncology III' = 'blue', 'Olink Target 96 Immune Response' = 'red')
    this_plot <- ggplot() + 
      geom_point(full_set, mapping=aes(x=Time, y=NPX, color=Panel), size=2) + geom_line(full_set, mapping=aes(x=Time, y=NPX, color=Panel), linetype='longdash') +
      theme_bw() +  xlab('Time (days)') + scale_color_manual(values = cols) +
      ggtitle(paste(this_assay, ' for Subject ', t_subject, ' : Target96 vs ExploreHT', sep=''))
    named <- paste(this_assay, '-', t_subject, '.png', sep='')
    #ggsave(path = 'HT_T96__figures', filename=named)
  }}

end <- Sys.time()
print('end Part 4 (ExploreHT vs Target96 Figures) at')
print(end)
