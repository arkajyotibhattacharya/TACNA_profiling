
library(data.table)

#Load two expression datasets
gdsc_expression_dataset = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ GDSC/Data/GDSC__Affy_hgu219_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetset_common_genes_with_TCGA.txt", sep = "\t", header = TRUE)

tcga_expression_dataset = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA__RSEM_genes_RNAseq__duplicate_samples_removed__genes_with_all_zeroes_removed_common_genes_with_GDSC.txt", sep = "\t", header = TRUE)

#Standardize two expression datasets
gdsc_expression_dataset_standardized = t(scale(t(gdsc_expression_dataset)))

tcga_expression_dataset_standardized = t(scale(t(tcga_expression_dataset)))

#export two standardized expression datasets
write.table(gdsc_expression_dataset_standardized, file = "/Users/arkajyotibhattacharya/Projects/Analysis\ on\ gdsc/Data/gdsc__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_common_genes_with_TCGA_genelevel_standardized.txt", sep = "\t", row.names = TRUE)
write.table(tcga_expression_dataset_standardized, file = "/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/TCGA__RSEM_genes_RNAseq__duplicate_samples_removed__genes_with_all_zeroes_removed_common_genes_with_gdsc_genelevel_standardized.txt", sep = "\t", row.names = TRUE)

#Load consensus estimated sources of TCGA dataset
tcga_ces = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/Consensus_Independent_Components_20170906_Duplicate_removed_TCGA_data_common_genes_with_gdsc.txt", sep = "\t", header = TRUE)

#merge two standardized expression datasets
combined_expression_dataset = cbind(tcga_expression_dataset_standardized, gdsc_expression_dataset_standardized)
tcga_ces = as.matrix(tcga_ces)

#generate mixing matrix for merged expression dataset
dot_prod = t(tcga_ces)%*%tcga_ces
V = solve(dot_prod)
mix_matrix = V%*%t(tcga_ces)%*%as.matrix(combined_expression_dataset)
rownames(mix_matrix) = colnames(tcga_ces)

#export mixing matrix for merged expression dataset
write.table(mix_matrix, file = "/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/TCGA\ CES\ on\ GDSC/Mixing_matrix_TCGA_GDSC_combined_using_TCGA_ces.txt", sep = "\t", row.names = TRUE)

#load extreme valued genomic region information of TCGA consensus estimated sources
evr_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ TCGA/Data/Genelevel_20170906_Duplicate_removed_TCGA_data_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_.txt", sep = "\t", header = TRUE))

#Create an indicator matrix with same dimension of the consensus estimated sources matrix having value 1 if a gene falls in an extreme valued genomic region in a consensus estimated source and 0 otherwise
tcga_ica_amp_del_indicator = matrix(0, dim(tcga_ces)[1], dim(tcga_ces)[2])
rownames(tcga_ica_amp_del_indicator) = rownames(tcga_ces)
colnames(tcga_ica_amp_del_indicator) = colnames(tcga_ces)

match_vec = match(evr_info$ENTREZID, rownames(tcga_ica_amp_del_indicator))

for(j in 1:dim(evr_info)[1])
{
  tcga_ica_amp_del_indicator[match_vec[j],evr_info$Component[j]] = abs(evr_info$Amp_or_Del[j])
}


#obtain CNA-CES (copy number alteration consensus estimated sources) using the indicator matrix created above
tcga_cna_ces = tcga_ces*tcga_ica_amp_del_indicator

#obtain TACNA profiles of merged expression datasets using the CNA-CESs of TCGA
combined_tacna = tcga_cna_ces%*%mix_matrix

#obtain subset of TACNA profiles for samples from GDSC using the CNA-CESs of TCGA
gdsc_tacna_by_tcga_ces = combined_tacna[,10818:dim(combined_tacna)[2]]

#load CN dataset (derived from SNP) for samples from GDSC
gdsc_snp_data = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ GDSC/Data/GDSC__median_centered_copy_number_data__genes_ranked_on_genomic_mapping.txt", sep = "\t", header = TRUE))

#convert CN dataset from probeset-level to gene-level
jetset_mapping_file = read.csv("/Users/arkajyotibhattacharya/Projects/Databases/Genomic\ mapping\ files/jetset.scores.hgu219_3.4.0.csv")
jetset_mapping_file = jetset_mapping_file[which(jetset_mapping_file$best==TRUE),]
rownames(jetset_mapping_file) = jetset_mapping_file$probeset
gdsc_snp_data = gdsc_snp_data[which(gdsc_snp_data$PROBESET%in%jetset_mapping_file$probeset),]
rownames(gdsc_snp_data) = gdsc_snp_data$PROBESET
jetset_mapping_file = jetset_mapping_file[which(jetset_mapping_file$probeset%in%gdsc_snp_data$PROBESET),]
gdsc_snp_data = gdsc_snp_data[rownames(jetset_mapping_file),]
rownames(gdsc_snp_data) = jetset_mapping_file$EntrezID

#list of genes common between SNP data and TACNA profiles for GDSC using CNA-CES of TCGA
common_genes_snp_tacna = intersect(rownames(gdsc_snp_data), rownames(gdsc_tacna_by_tcga_ces))

#subsetting both SNP data and TACNA profiles of GDSC using CNA-CES of TCGA with common genes
gdsc_snp_data = gdsc_snp_data[common_genes_snp_tacna,]
gdsc_tacna_by_tcga_ces = gdsc_tacna_by_tcga_ces[common_genes_snp_tacna,]

#subsetting both SNP data and TACNA profiles of GDSC using CNA-CES of TCGA with samples present in both the datasets
common_samples = intersect(colnames(gdsc_snp_data),colnames(gdsc_tacna_by_tcga_ces))
gdsc_snp_data = gdsc_snp_data[,common_samples]
gdsc_tacna_by_tcga_ces = gdsc_tacna_by_tcga_ces[,common_samples]


#Obtain correlation between TACNA profiles using CNA-CES of TCGA and SNP profiles for GDSC
gdsc_snp_data = as.matrix(gdsc_snp_data)
cor_between_snp_and_tacna_by_tcga = NULL

for(i in 1:dim(gdsc_tacna_by_tcga_ces)[2])
{
  cor_between_snp_and_tacna_by_tcga[i] = cor(gdsc_snp_data[,i], gdsc_tacna_by_tcga_ces[,i], use = "pairwise.complete.obs")
}

cor_between_snp_and_tacna_by_tcga = as.data.frame(cbind(common_samples,cor_between_snp_and_tacna_by_tcga))

#Export file containing correlation coefficients
write.table(cor_between_snp_and_tacna_by_tcga, file = "/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Tables/Correlation_of_TACNA_with_SNP_of_GDSC_by_TCGA_CES.txt", sep = "\t", row.names = FALSE)

#Plot histogram of the above correlation coefficients
pdf("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Plots/Histogram of correlation between SNP and regenerated GDSC TACNA using TCGA CNA-CES.pdf", height = 6, width = 12)
hist(as.numeric(as.character(cor_between_snp_and_tacna_by_tcga$cor_between_snp_and_tacna_by_tcga)), breaks = 100,main = "Histogram of correlation between SNP and regenerated GDSC TACNA using TCGA CNA-CES"
     ,xlab = "Pearson correlation coefficient")
     
dev.off()



