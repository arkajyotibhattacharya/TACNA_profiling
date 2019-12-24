setwd("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/")
library(sqldf)

total_time = proc.time()[3]
time1 = proc.time()[3]

#load expression dataset
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
expression_data <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt" , header = TRUE, colClasses = classes)

print(head(rownames(expression_data)))
print(head(colnames(expression_data)))

#Distribute samples into 5 groups using multinomial distribution simulation 
set.seed(12345)
order_of_folds = sample(1:5,dim(expression_data)[2], replace = T)
table(order_of_folds)

#Create 5 folds of input expression datasets by removing samples belonging to one group in each fold
expression_data_fold1 = expression_data[,which(order_of_folds!=1)]
expression_data_fold2 = expression_data[,which(order_of_folds!=2)]
expression_data_fold3 = expression_data[,which(order_of_folds!=3)]
expression_data_fold4 = expression_data[,which(order_of_folds!=4)]
expression_data_fold5 = expression_data[,which(order_of_folds!=5)]


#export all the input datasets
write.table(expression_data_fold1, file = "/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold1.txt", sep = "\t", row.names = TRUE)
write.table(expression_data_fold2, file = "/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold2.txt", sep = "\t", row.names = TRUE)
write.table(expression_data_fold3, file = "/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold3.txt", sep = "\t", row.names = TRUE)
write.table(expression_data_fold4, file = "/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold4.txt", sep = "\t", row.names = TRUE)
write.table(expression_data_fold5, file = "/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold5.txt", sep = "\t", row.names = TRUE)

#TACNA profiling of fold1

#set the working directory
setwd("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Results/CCLE_cross_validation/")
library(sqldf)

total_time = proc.time()[3]
time1 = proc.time()[3]
#load expression data fold1
tab5rows  = read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold1.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
expression_data <- read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold1.txt" , header = TRUE, colClasses = classes)

print(head(rownames(expression_data)))
print(head(colnames(expression_data)))


#load genomic mapping file
data_with_chr_bp_mapping = read.csv("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/Final Non-missing 1-1 Genomic Mapping 20161130 hgu133plus2 Genelevel.csv")

print(paste("time taken to load the data is", round((proc.time()[3]-time1)/60,2),"mins"))

#load source TACNA profiling code
source("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Codes/Source code TACNA.R")

#Run the TACNA profiling source code
combat_consensus_ica_degr_clustering(expression_data = expression_data,probelevel_standardization = TRUE,
                                combat_identifier = FALSE,
                                Title="20191030_CCLE_data_fold1",send_email_indicator = FALSE,
                                email_id="arkajyoti.bhattacharya@gmail.com",
                                CES_clustering_algo1 = TRUE,
                                batch_correction_check =FALSE,n_princomp_check = TRUE,
                                prin_comp_cor = FALSE,choose_fastICA_ncomp = FALSE,
                                var_cutoff_to_choose_ncomp = 0.85,ncomp_without_pca=1460,
                                nsample_fastica = 25,
                                seed = 12345678,fastica_alg_typ = "parallel",
                                fastica_fun = "logcosh", fastica_alpha = 1,fastica_method = "R",
                                fastica_row_norm = FALSE, fastica_maxit = 1000,
                                fastica_tol = 0.0001, fastica_verbose = FALSE,CES_clustering = TRUE,
                                consensus_clustering_of_data = FALSE,
                                distance_function = 'pearson',no_mad_ge = 5000,
                                consensus_maxK=20,consensus_reps=100,consensus_pItem=0.8,
                                consensus_pFeature=1,consensus_clusterAlg="hc",
                                permutation_testing_for_degr_indicator = TRUE,
                                gene_level_mapping_indicator=data_with_chr_bp_mapping,
                                intervals_for_gaus_kernel = seq(10000,2000000,by=10000),
                                probe_no_for_gaus_kernel = 10,
                                FDR = c(0.01,0.05,0.1,0.2),
                                CL = 0.5,
                                state_deciding_cutoff = 0.95,
                                min_probesets = 5,
                                set_seed = 12345678,
                                gene_level_mapping_indicator = TRUE, TACNA_profiling = TRUE)



print(paste("time taken to do the whole analysis is", round((proc.time()[3]-total_time)/60,2),"mins"))

#Load the expression dataset for all samples from CCLE
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
ccle_expression_dataset <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt" , header = TRUE, colClasses = classes)


#Load consensus estimated sources matrix for fold1
ces_fold1 = read.table("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Consensus_Independent_Components_20191030_CCLE_data_fold1_.txt", sep = "\t", header = TRUE)
#Subset the expression dataset with genes from the consensus estimated sources matrix (in this case both are same)
ccle_expression_dataset = ccle_expression_dataset[rownames(ces_fold1),]
#standardize the expression dataset
ccle_expression_dataset_standardized = t(scale(t(ccle_expression_dataset)))

#Obtain mixing matrix weights for samples from CCLE using CES of only the samples from fold1 expression dataset
ces_fold1 = as.matrix(ces_fold1)
dot_prod = t(ces_fold1)%*%ces_fold1
V = solve(dot_prod)
mix_matrix = V%*%t(ces_fold1)%*%as.matrix(ccle_expression_dataset_standardized)
rownames(mix_matrix) = colnames(ces_fold1)

#load extreme valued genomic region information of consensus estimated sources
evr_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/SCNA_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_probe_no_for_gaus_kernel_10/Genelevel_20191030_CCLE_data_fold1_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_.txt", sep = "\t", header = TRUE))

#Create an indicator matrix with same dimension of the consensus estimated sources matrix having value 1 if a gene falls in an extreme valued genomic region in a consensus estimated source and 0 otherwise
ccle_ica_amp_del_indicator = matrix(0, dim(ces_fold1)[1], dim(ces_fold1)[2])
rownames(ccle_ica_amp_del_indicator) = rownames(ces_fold1)
colnames(ccle_ica_amp_del_indicator) = colnames(ces_fold1)
match_vec = match(evr_info$ENTREZID, rownames(ccle_ica_amp_del_indicator))
for(j in 1:dim(evr_info)[1])
{
  ccle_ica_amp_del_indicator[match_vec[j],evr_info$Component[j]] = abs(evr_info$Amp_or_Del[j])
}

#obtain CNA-CES (copy number alteration consensus estimated sources) using the indicator matrix created above
ccle_cna_ces_fold1 = ces_fold1*ccle_ica_amp_del_indicator

#obtain TACNA profiles of samples from CCLE using the CNA-CESs of samples ony from fold1
ccle_tacna_fold1 = ccle_cna_ces_fold1%*%mix_matrix

#load CN dataset (derived from SNP) for samples from CCLE
ccle_snp_data = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__median_centered_copy_number_data__genes_ranked_on_genomic_mapping.txt", sep = "\t", header = TRUE)

#convert CN dataset from probeset-level to gene-level
jetset_mapping = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Databases/Genomic\ mapping\ files/Genomic_Mapping_hgu133plus2_using_jetscore_30032018.txt", sep = "\t", header = TRUE))
rownames(jetset_mapping) = jetset_mapping$PROBESET
jetset_mapping = jetset_mapping[which(jetset_mapping$PROBESET%in%ccle_snp_data$PROBESET),]
jetset_mapping = jetset_mapping[which(jetset_mapping$top_probe_indicator==1),]
rownames(ccle_snp_data) = ccle_snp_data$PROBESET
ccle_snp_data = ccle_snp_data[rownames(jetset_mapping),]
rownames(ccle_snp_data) = jetset_mapping$ENTREZID


#list of genes common between SNP data and TACNA profiles of CCLE using CNA-CES of fold1
common_genes_snp_tacna = intersect(rownames(ccle_snp_data), rownames(ccle_tacna_fold1))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold1 with common genes
ccle_snp_data = ccle_snp_data[common_genes_snp_tacna,]
ccle_tacna_fold1 = ccle_tacna_fold1[common_genes_snp_tacna,]

common_samples = intersect(colnames(ccle_snp_data),colnames(ccle_tacna_fold1))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold1 with common samples
ccle_snp_data = ccle_snp_data[,common_samples]
ccle_tacna_fold1 = ccle_tacna_fold1[,common_samples]

#Obtain correlation between TACNA profiles using CNA-CES of fold1 and SNP profiles CCLE
ccle_snp_data = as.matrix(ccle_snp_data)
cor_between_snp_and_tacna_by_ccle = NULL

for(i in 1:dim(ccle_tacna_fold1)[2])
{
  cor_between_snp_and_tacna_by_ccle[i] = cor(ccle_snp_data[,i], ccle_tacna_fold1[,i], use = "pairwise.complete.obs")
}

#load list of samples used in fold1
samples_in_fold1 = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold1.txt", sep = "\t" , header = TRUE, nrows = 5)

#correlation coefficients for samples not belonging to fold1
cor_between_snp_and_tacna_fold1 = cor_between_snp_and_tacna_by_ccle[which(!(common_samples%in%colnames(samples_in_fold1)))]

samples_in_validation_fold1 = common_samples[which(!(common_samples%in%colnames(samples_in_fold1)))]

correlation_fold1 = as.data.frame(cbind(samples_in_validation_fold1, cor_between_snp_and_tacna_fold1))
correlation_fold1$cor_between_snp_and_tacna_fold1 = as.numeric(as.character(correlation_fold1$cor_between_snp_and_tacna_fold1))
#Export file containing correlation coefficients
write.table(correlation_fold1, file = "/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Correlation between TACNA and SNP fold1.txt", sep = "\t", row.names = FALSE)

#Plot histogram of the above correlation coefficients
pdf("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Plots/Histogram of correlation between SNP and regenerated CCLE TACNA fold1.pdf", height = 6, width = 12)
hist(cor_between_snp_and_tacna_fold1, breaks = 100,main = "Histogram of correlation between SNP and regenerated CCLE TACNA fold1"
     ,xlab = "Pearson correlation coefficient")

dev.off()





#TACNA profiling of fold2

#set the working directory
setwd("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Results/CCLE_cross_validation/")
library(sqldf)

total_time = proc.time()[3]
time1 = proc.time()[3]
#load expression data fold2
tab5rows  = read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold2.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
expression_data <- read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold2.txt" , header = TRUE, colClasses = classes)

print(head(rownames(expression_data)))
print(head(colnames(expression_data)))


#load genomic mapping file
data_with_chr_bp_mapping = read.csv("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/Final Non-missing 1-1 Genomic Mapping 20161130 hgu133plus2 Genelevel.csv")

print(paste("time taken to load the data is", round((proc.time()[3]-time1)/60,2),"mins"))

#load source TACNA profiling code
source("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Codes/Source code TACNA.R")

#Run the TACNA profiling source code
combat_consensus_ica_degr_clustering(expression_data = expression_data,probelevel_standardization = TRUE,
                                     combat_identifier = FALSE,
                                     Title="20191030_CCLE_data_fold2",send_email_indicator = FALSE,
                                     email_id="arkajyoti.bhattacharya@gmail.com",
                                     CES_clustering_algo1 = TRUE,
                                     batch_correction_check =FALSE,n_princomp_check = TRUE,
                                     prin_comp_cor = FALSE,choose_fastICA_ncomp = FALSE,
                                     var_cutoff_to_choose_ncomp = 0.85,ncomp_without_pca=1460,
                                     nsample_fastica = 25,
                                     seed = 12345678,fastica_alg_typ = "parallel",
                                     fastica_fun = "logcosh", fastica_alpha = 1,fastica_method = "R",
                                     fastica_row_norm = FALSE, fastica_maxit = 1000,
                                     fastica_tol = 0.0001, fastica_verbose = FALSE,CES_clustering = TRUE,
                                     consensus_clustering_of_data = FALSE,
                                     distance_function = 'pearson',no_mad_ge = 5000,
                                     consensus_maxK=20,consensus_reps=100,consensus_pItem=0.8,
                                     consensus_pFeature=1,consensus_clusterAlg="hc",
                                     permutation_testing_for_degr_indicator = TRUE,
                                     gene_level_mapping_indicator=data_with_chr_bp_mapping,
                                     intervals_for_gaus_kernel = seq(10000,2000000,by=10000),
                                     probe_no_for_gaus_kernel = 10,
                                     FDR = c(0.01,0.05,0.1,0.2),
                                     CL = 0.5,
                                     state_deciding_cutoff = 0.95,
                                     min_probesets = 5,
                                     set_seed = 12345678,
                                     gene_level_mapping_indicator = TRUE, TACNA_profiling = TRUE)



print(paste("time taken to do the whole analysis is", round((proc.time()[3]-total_time)/60,2),"mins"))

#Load the expression dataset for all samples from CCLE
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
ccle_expression_dataset <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt" , header = TRUE, colClasses = classes)


#Load consensus estimated sources matrix for fold2
ces_fold2 = read.table("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Consensus_Independent_Components_20191030_CCLE_data_fold2_.txt", sep = "\t", header = TRUE)
#Subset the expression dataset with genes from the consensus estimated sources matrix (in this case both are same)
ccle_expression_dataset = ccle_expression_dataset[rownames(ces_fold2),]
#standardize the expression dataset
ccle_expression_dataset_standardized = t(scale(t(ccle_expression_dataset)))

#Obtain mixing matrix weights for samples from CCLE using CES of only the samples from fold2 expression dataset
ces_fold2 = as.matrix(ces_fold2)
dot_prod = t(ces_fold2)%*%ces_fold2
V = solve(dot_prod)
mix_matrix = V%*%t(ces_fold2)%*%as.matrix(ccle_expression_dataset_standardized)
rownames(mix_matrix) = colnames(ces_fold2)

#load extreme valued genomic region information of consensus estimated sources
evr_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/SCNA_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_probe_no_for_gaus_kernel_10/Genelevel_20191030_CCLE_data_fold2_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_.txt", sep = "\t", header = TRUE))

#Create an indicator matrix with same dimension of the consensus estimated sources matrix having value 1 if a gene falls in an extreme valued genomic region in a consensus estimated source and 0 otherwise
ccle_ica_amp_del_indicator = matrix(0, dim(ces_fold2)[1], dim(ces_fold2)[2])
rownames(ccle_ica_amp_del_indicator) = rownames(ces_fold2)
colnames(ccle_ica_amp_del_indicator) = colnames(ces_fold2)
match_vec = match(evr_info$ENTREZID, rownames(ccle_ica_amp_del_indicator))
for(j in 1:dim(evr_info)[1])
{
  ccle_ica_amp_del_indicator[match_vec[j],evr_info$Component[j]] = abs(evr_info$Amp_or_Del[j])
}

#obtain CNA-CES (copy number alteration consensus estimated sources) using the indicator matrix created above
ccle_cna_ces_fold2 = ces_fold2*ccle_ica_amp_del_indicator

#obtain TACNA profiles of samples from CCLE using the CNA-CESs of samples ony from fold2
ccle_tacna_fold2 = ccle_cna_ces_fold2%*%mix_matrix

#load CN dataset (derived from SNP) for samples from CCLE
ccle_snp_data = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__median_centered_copy_number_data__genes_ranked_on_genomic_mapping.txt", sep = "\t", header = TRUE)

#convert CN dataset from probeset-level to gene-level
jetset_mapping = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Databases/Genomic\ mapping\ files/Genomic_Mapping_hgu133plus2_using_jetscore_30032018.txt", sep = "\t", header = TRUE))
rownames(jetset_mapping) = jetset_mapping$PROBESET
jetset_mapping = jetset_mapping[which(jetset_mapping$PROBESET%in%ccle_snp_data$PROBESET),]
jetset_mapping = jetset_mapping[which(jetset_mapping$top_probe_indicator==1),]
rownames(ccle_snp_data) = ccle_snp_data$PROBESET
ccle_snp_data = ccle_snp_data[rownames(jetset_mapping),]
rownames(ccle_snp_data) = jetset_mapping$ENTREZID


#list of genes common between SNP data and TACNA profiles of CCLE using CNA-CES of fold2
common_genes_snp_tacna = intersect(rownames(ccle_snp_data), rownames(ccle_tacna_fold2))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold2 with common genes
ccle_snp_data = ccle_snp_data[common_genes_snp_tacna,]
ccle_tacna_fold2 = ccle_tacna_fold2[common_genes_snp_tacna,]

common_samples = intersect(colnames(ccle_snp_data),colnames(ccle_tacna_fold2))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold2 with common samples
ccle_snp_data = ccle_snp_data[,common_samples]
ccle_tacna_fold2 = ccle_tacna_fold2[,common_samples]

#Obtain correlation between TACNA profiles using CNA-CES of fold2 and SNP profiles CCLE
ccle_snp_data = as.matrix(ccle_snp_data)
cor_between_snp_and_tacna_by_ccle = NULL

for(i in 1:dim(ccle_tacna_fold2)[2])
{
  cor_between_snp_and_tacna_by_ccle[i] = cor(ccle_snp_data[,i], ccle_tacna_fold2[,i], use = "pairwise.complete.obs")
}

#load list of samples used in fold2
samples_in_fold2 = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold2.txt", sep = "\t" , header = TRUE, nrows = 5)

#correlation coefficients for samples not belonging to fold2
cor_between_snp_and_tacna_fold2 = cor_between_snp_and_tacna_by_ccle[which(!(common_samples%in%colnames(samples_in_fold2)))]

samples_in_validation_fold2 = common_samples[which(!(common_samples%in%colnames(samples_in_fold2)))]

correlation_fold2 = as.data.frame(cbind(samples_in_validation_fold2, cor_between_snp_and_tacna_fold2))
correlation_fold2$cor_between_snp_and_tacna_fold2 = as.numeric(as.character(correlation_fold2$cor_between_snp_and_tacna_fold2))
#Export file containing correlation coefficients
write.table(correlation_fold2, file = "/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Correlation between TACNA and SNP fold2.txt", sep = "\t", row.names = FALSE)

#Plot histogram of the above correlation coefficients
pdf("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Plots/Histogram of correlation between SNP and regenerated CCLE TACNA fold2.pdf", height = 6, width = 12)
hist(cor_between_snp_and_tacna_fold2, breaks = 100,main = "Histogram of correlation between SNP and regenerated CCLE TACNA fold2"
     ,xlab = "Pearson correlation coefficient")

dev.off()






#TACNA profiling of fold3

#set the working directory
setwd("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Results/CCLE_cross_validation/")
library(sqldf)

total_time = proc.time()[3]
time1 = proc.time()[3]
#load expression data fold3
tab5rows  = read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold3.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
expression_data <- read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold3.txt" , header = TRUE, colClasses = classes)

print(head(rownames(expression_data)))
print(head(colnames(expression_data)))


#load genomic mapping file
data_with_chr_bp_mapping = read.csv("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/Final Non-missing 1-1 Genomic Mapping 20161130 hgu133plus2 Genelevel.csv")

print(paste("time taken to load the data is", round((proc.time()[3]-time1)/60,2),"mins"))

#load source TACNA profiling code
source("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Codes/Source code TACNA.R")

#Run the TACNA profiling source code
combat_consensus_ica_degr_clustering(expression_data = expression_data,probelevel_standardization = TRUE,
                                     combat_identifier = FALSE,
                                     Title="20191030_CCLE_data_fold3",send_email_indicator = FALSE,
                                     email_id="arkajyoti.bhattacharya@gmail.com",
                                     CES_clustering_algo1 = TRUE,
                                     batch_correction_check =FALSE,n_princomp_check = TRUE,
                                     prin_comp_cor = FALSE,choose_fastICA_ncomp = FALSE,
                                     var_cutoff_to_choose_ncomp = 0.85,ncomp_without_pca=1460,
                                     nsample_fastica = 25,
                                     seed = 12345678,fastica_alg_typ = "parallel",
                                     fastica_fun = "logcosh", fastica_alpha = 1,fastica_method = "R",
                                     fastica_row_norm = FALSE, fastica_maxit = 1000,
                                     fastica_tol = 0.0001, fastica_verbose = FALSE,CES_clustering = TRUE,
                                     consensus_clustering_of_data = FALSE,
                                     distance_function = 'pearson',no_mad_ge = 5000,
                                     consensus_maxK=20,consensus_reps=100,consensus_pItem=0.8,
                                     consensus_pFeature=1,consensus_clusterAlg="hc",
                                     permutation_testing_for_degr_indicator = TRUE,
                                     gene_level_mapping_indicator=data_with_chr_bp_mapping,
                                     intervals_for_gaus_kernel = seq(10000,2000000,by=10000),
                                     probe_no_for_gaus_kernel = 10,
                                     FDR = c(0.01,0.05,0.1,0.2),
                                     CL = 0.5,
                                     state_deciding_cutoff = 0.95,
                                     min_probesets = 5,
                                     set_seed = 12345678,
                                     gene_level_mapping_indicator = TRUE, TACNA_profiling = TRUE)



print(paste("time taken to do the whole analysis is", round((proc.time()[3]-total_time)/60,2),"mins"))

#Load the expression dataset for all samples from CCLE
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
ccle_expression_dataset <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt" , header = TRUE, colClasses = classes)


#Load consensus estimated sources matrix for fold3
ces_fold3 = read.table("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Consensus_Independent_Components_20191030_CCLE_data_fold3_.txt", sep = "\t", header = TRUE)
#Subset the expression dataset with genes from the consensus estimated sources matrix (in this case both are same)
ccle_expression_dataset = ccle_expression_dataset[rownames(ces_fold3),]
#standardize the expression dataset
ccle_expression_dataset_standardized = t(scale(t(ccle_expression_dataset)))

#Obtain mixing matrix weights for samples from CCLE using CES of only the samples from fold3 expression dataset
ces_fold3 = as.matrix(ces_fold3)
dot_prod = t(ces_fold3)%*%ces_fold3
V = solve(dot_prod)
mix_matrix = V%*%t(ces_fold3)%*%as.matrix(ccle_expression_dataset_standardized)
rownames(mix_matrix) = colnames(ces_fold3)

#load extreme valued genomic region information of consensus estimated sources
evr_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/SCNA_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_probe_no_for_gaus_kernel_10/Genelevel_20191030_CCLE_data_fold3_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_.txt", sep = "\t", header = TRUE))

#Create an indicator matrix with same dimension of the consensus estimated sources matrix having value 1 if a gene falls in an extreme valued genomic region in a consensus estimated source and 0 otherwise
ccle_ica_amp_del_indicator = matrix(0, dim(ces_fold3)[1], dim(ces_fold3)[2])
rownames(ccle_ica_amp_del_indicator) = rownames(ces_fold3)
colnames(ccle_ica_amp_del_indicator) = colnames(ces_fold3)
match_vec = match(evr_info$ENTREZID, rownames(ccle_ica_amp_del_indicator))
for(j in 1:dim(evr_info)[1])
{
  ccle_ica_amp_del_indicator[match_vec[j],evr_info$Component[j]] = abs(evr_info$Amp_or_Del[j])
}

#obtain CNA-CES (copy number alteration consensus estimated sources) using the indicator matrix created above
ccle_cna_ces_fold3 = ces_fold3*ccle_ica_amp_del_indicator

#obtain TACNA profiles of samples from CCLE using the CNA-CESs of samples ony from fold3
ccle_tacna_fold3 = ccle_cna_ces_fold3%*%mix_matrix

#load CN dataset (derived from SNP) for samples from CCLE
ccle_snp_data = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__median_centered_copy_number_data__genes_ranked_on_genomic_mapping.txt", sep = "\t", header = TRUE)

#convert CN dataset from probeset-level to gene-level
jetset_mapping = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Databases/Genomic\ mapping\ files/Genomic_Mapping_hgu133plus2_using_jetscore_30032018.txt", sep = "\t", header = TRUE))
rownames(jetset_mapping) = jetset_mapping$PROBESET
jetset_mapping = jetset_mapping[which(jetset_mapping$PROBESET%in%ccle_snp_data$PROBESET),]
jetset_mapping = jetset_mapping[which(jetset_mapping$top_probe_indicator==1),]
rownames(ccle_snp_data) = ccle_snp_data$PROBESET
ccle_snp_data = ccle_snp_data[rownames(jetset_mapping),]
rownames(ccle_snp_data) = jetset_mapping$ENTREZID


#list of genes common between SNP data and TACNA profiles of CCLE using CNA-CES of fold3
common_genes_snp_tacna = intersect(rownames(ccle_snp_data), rownames(ccle_tacna_fold3))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold3 with common genes
ccle_snp_data = ccle_snp_data[common_genes_snp_tacna,]
ccle_tacna_fold3 = ccle_tacna_fold3[common_genes_snp_tacna,]

common_samples = intersect(colnames(ccle_snp_data),colnames(ccle_tacna_fold3))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold3 with common samples
ccle_snp_data = ccle_snp_data[,common_samples]
ccle_tacna_fold3 = ccle_tacna_fold3[,common_samples]

#Obtain correlation between TACNA profiles using CNA-CES of fold3 and SNP profiles CCLE
ccle_snp_data = as.matrix(ccle_snp_data)
cor_between_snp_and_tacna_by_ccle = NULL

for(i in 1:dim(ccle_tacna_fold3)[2])
{
  cor_between_snp_and_tacna_by_ccle[i] = cor(ccle_snp_data[,i], ccle_tacna_fold3[,i], use = "pairwise.complete.obs")
}

#load list of samples used in fold3
samples_in_fold3 = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold3.txt", sep = "\t" , header = TRUE, nrows = 5)

#correlation coefficients for samples not belonging to fold3
cor_between_snp_and_tacna_fold3 = cor_between_snp_and_tacna_by_ccle[which(!(common_samples%in%colnames(samples_in_fold3)))]

samples_in_validation_fold3 = common_samples[which(!(common_samples%in%colnames(samples_in_fold3)))]

correlation_fold3 = as.data.frame(cbind(samples_in_validation_fold3, cor_between_snp_and_tacna_fold3))
correlation_fold3$cor_between_snp_and_tacna_fold3 = as.numeric(as.character(correlation_fold3$cor_between_snp_and_tacna_fold3))
#Export file containing correlation coefficients
write.table(correlation_fold3, file = "/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Correlation between TACNA and SNP fold3.txt", sep = "\t", row.names = FALSE)

#Plot histogram of the above correlation coefficients
pdf("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Plots/Histogram of correlation between SNP and regenerated CCLE TACNA fold3.pdf", height = 6, width = 12)
hist(cor_between_snp_and_tacna_fold3, breaks = 100,main = "Histogram of correlation between SNP and regenerated CCLE TACNA fold3"
     ,xlab = "Pearson correlation coefficient")

dev.off()





#TACNA profiling of fold4

#set the working directory
setwd("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Results/CCLE_cross_validation/")
library(sqldf)

total_time = proc.time()[3]
time1 = proc.time()[3]
#load expression data fold4
tab5rows  = read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold4.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
expression_data <- read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold4.txt" , header = TRUE, colClasses = classes)

print(head(rownames(expression_data)))
print(head(colnames(expression_data)))


#load genomic mapping file
data_with_chr_bp_mapping = read.csv("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/Final Non-missing 1-1 Genomic Mapping 20161130 hgu133plus2 Genelevel.csv")

print(paste("time taken to load the data is", round((proc.time()[3]-time1)/60,2),"mins"))

#load source TACNA profiling code
source("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Codes/Source code TACNA.R")

#Run the TACNA profiling source code
combat_consensus_ica_degr_clustering(expression_data = expression_data,probelevel_standardization = TRUE,
                                     combat_identifier = FALSE,
                                     Title="20191030_CCLE_data_fold4",send_email_indicator = FALSE,
                                     email_id="arkajyoti.bhattacharya@gmail.com",
                                     CES_clustering_algo1 = TRUE,
                                     batch_correction_check =FALSE,n_princomp_check = TRUE,
                                     prin_comp_cor = FALSE,choose_fastICA_ncomp = FALSE,
                                     var_cutoff_to_choose_ncomp = 0.85,ncomp_without_pca=1460,
                                     nsample_fastica = 25,
                                     seed = 12345678,fastica_alg_typ = "parallel",
                                     fastica_fun = "logcosh", fastica_alpha = 1,fastica_method = "R",
                                     fastica_row_norm = FALSE, fastica_maxit = 1000,
                                     fastica_tol = 0.0001, fastica_verbose = FALSE,CES_clustering = TRUE,
                                     consensus_clustering_of_data = FALSE,
                                     distance_function = 'pearson',no_mad_ge = 5000,
                                     consensus_maxK=20,consensus_reps=100,consensus_pItem=0.8,
                                     consensus_pFeature=1,consensus_clusterAlg="hc",
                                     permutation_testing_for_degr_indicator = TRUE,
                                     gene_level_mapping_indicator=data_with_chr_bp_mapping,
                                     intervals_for_gaus_kernel = seq(10000,2000000,by=10000),
                                     probe_no_for_gaus_kernel = 10,
                                     FDR = c(0.01,0.05,0.1,0.2),
                                     CL = 0.5,
                                     state_deciding_cutoff = 0.95,
                                     min_probesets = 5,
                                     set_seed = 12345678,
                                     gene_level_mapping_indicator = TRUE, TACNA_profiling = TRUE)



print(paste("time taken to do the whole analysis is", round((proc.time()[3]-total_time)/60,2),"mins"))

#Load the expression dataset for all samples from CCLE
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
ccle_expression_dataset <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt" , header = TRUE, colClasses = classes)


#Load consensus estimated sources matrix for fold4
ces_fold4 = read.table("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Consensus_Independent_Components_20191030_CCLE_data_fold4_.txt", sep = "\t", header = TRUE)
#Subset the expression dataset with genes from the consensus estimated sources matrix (in this case both are same)
ccle_expression_dataset = ccle_expression_dataset[rownames(ces_fold4),]
#standardize the expression dataset
ccle_expression_dataset_standardized = t(scale(t(ccle_expression_dataset)))

#Obtain mixing matrix weights for samples from CCLE using CES of only the samples from fold4 expression dataset
ces_fold4 = as.matrix(ces_fold4)
dot_prod = t(ces_fold4)%*%ces_fold4
V = solve(dot_prod)
mix_matrix = V%*%t(ces_fold4)%*%as.matrix(ccle_expression_dataset_standardized)
rownames(mix_matrix) = colnames(ces_fold4)

#load extreme valued genomic region information of consensus estimated sources
evr_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/SCNA_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_probe_no_for_gaus_kernel_10/Genelevel_20191030_CCLE_data_fold4_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_.txt", sep = "\t", header = TRUE))

#Create an indicator matrix with same dimension of the consensus estimated sources matrix having value 1 if a gene falls in an extreme valued genomic region in a consensus estimated source and 0 otherwise
ccle_ica_amp_del_indicator = matrix(0, dim(ces_fold4)[1], dim(ces_fold4)[2])
rownames(ccle_ica_amp_del_indicator) = rownames(ces_fold4)
colnames(ccle_ica_amp_del_indicator) = colnames(ces_fold4)
match_vec = match(evr_info$ENTREZID, rownames(ccle_ica_amp_del_indicator))
for(j in 1:dim(evr_info)[1])
{
  ccle_ica_amp_del_indicator[match_vec[j],evr_info$Component[j]] = abs(evr_info$Amp_or_Del[j])
}

#obtain CNA-CES (copy number alteration consensus estimated sources) using the indicator matrix created above
ccle_cna_ces_fold4 = ces_fold4*ccle_ica_amp_del_indicator

#obtain TACNA profiles of samples from CCLE using the CNA-CESs of samples ony from fold4
ccle_tacna_fold4 = ccle_cna_ces_fold4%*%mix_matrix

#load CN dataset (derived from SNP) for samples from CCLE
ccle_snp_data = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__median_centered_copy_number_data__genes_ranked_on_genomic_mapping.txt", sep = "\t", header = TRUE)

#convert CN dataset from probeset-level to gene-level
jetset_mapping = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Databases/Genomic\ mapping\ files/Genomic_Mapping_hgu133plus2_using_jetscore_30032018.txt", sep = "\t", header = TRUE))
rownames(jetset_mapping) = jetset_mapping$PROBESET
jetset_mapping = jetset_mapping[which(jetset_mapping$PROBESET%in%ccle_snp_data$PROBESET),]
jetset_mapping = jetset_mapping[which(jetset_mapping$top_probe_indicator==1),]
rownames(ccle_snp_data) = ccle_snp_data$PROBESET
ccle_snp_data = ccle_snp_data[rownames(jetset_mapping),]
rownames(ccle_snp_data) = jetset_mapping$ENTREZID


#list of genes common between SNP data and TACNA profiles of CCLE using CNA-CES of fold4
common_genes_snp_tacna = intersect(rownames(ccle_snp_data), rownames(ccle_tacna_fold4))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold4 with common genes
ccle_snp_data = ccle_snp_data[common_genes_snp_tacna,]
ccle_tacna_fold4 = ccle_tacna_fold4[common_genes_snp_tacna,]

common_samples = intersect(colnames(ccle_snp_data),colnames(ccle_tacna_fold4))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold4 with common samples
ccle_snp_data = ccle_snp_data[,common_samples]
ccle_tacna_fold4 = ccle_tacna_fold4[,common_samples]

#Obtain correlation between TACNA profiles using CNA-CES of fold4 and SNP profiles CCLE
ccle_snp_data = as.matrix(ccle_snp_data)
cor_between_snp_and_tacna_by_ccle = NULL

for(i in 1:dim(ccle_tacna_fold4)[2])
{
  cor_between_snp_and_tacna_by_ccle[i] = cor(ccle_snp_data[,i], ccle_tacna_fold4[,i], use = "pairwise.complete.obs")
}

#load list of samples used in fold4
samples_in_fold4 = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold4.txt", sep = "\t" , header = TRUE, nrows = 5)

#correlation coefficients for samples not belonging to fold4
cor_between_snp_and_tacna_fold4 = cor_between_snp_and_tacna_by_ccle[which(!(common_samples%in%colnames(samples_in_fold4)))]

samples_in_validation_fold4 = common_samples[which(!(common_samples%in%colnames(samples_in_fold4)))]

correlation_fold4 = as.data.frame(cbind(samples_in_validation_fold4, cor_between_snp_and_tacna_fold4))
correlation_fold4$cor_between_snp_and_tacna_fold4 = as.numeric(as.character(correlation_fold4$cor_between_snp_and_tacna_fold4))
#Export file containing correlation coefficients
write.table(correlation_fold4, file = "/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Correlation between TACNA and SNP fold4.txt", sep = "\t", row.names = FALSE)

#Plot histogram of the above correlation coefficients
pdf("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Plots/Histogram of correlation between SNP and regenerated CCLE TACNA fold4.pdf", height = 6, width = 12)
hist(cor_between_snp_and_tacna_fold4, breaks = 100,main = "Histogram of correlation between SNP and regenerated CCLE TACNA fold4"
     ,xlab = "Pearson correlation coefficient")

dev.off()





#TACNA profiling of fold5

#set the working directory
setwd("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Results/CCLE_cross_validation/")
library(sqldf)

total_time = proc.time()[3]
time1 = proc.time()[3]
#load expression data fold5
tab5rows  = read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold5.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
expression_data <- read.table("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold5.txt" , header = TRUE, colClasses = classes)

print(head(rownames(expression_data)))
print(head(colnames(expression_data)))


#load genomic mapping file
data_with_chr_bp_mapping = read.csv("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Data/Final Non-missing 1-1 Genomic Mapping 20161130 hgu133plus2 Genelevel.csv")

print(paste("time taken to load the data is", round((proc.time()[3]-time1)/60,2),"mins"))

#load source TACNA profiling code
source("/mnt/43cc8b63-bef6-4251-a432-50f54f68bd94/Arkajyoti/TACNA resubmission/Codes/Source code TACNA.R")

#Run the TACNA profiling source code
combat_consensus_ica_degr_clustering(expression_data = expression_data,probelevel_standardization = TRUE,
                                     combat_identifier = FALSE,
                                     Title="20191030_CCLE_data_fold5",send_email_indicator = FALSE,
                                     email_id="arkajyoti.bhattacharya@gmail.com",
                                     CES_clustering_algo1 = TRUE,
                                     batch_correction_check =FALSE,n_princomp_check = TRUE,
                                     prin_comp_cor = FALSE,choose_fastICA_ncomp = FALSE,
                                     var_cutoff_to_choose_ncomp = 0.85,ncomp_without_pca=1460,
                                     nsample_fastica = 25,
                                     seed = 12345678,fastica_alg_typ = "parallel",
                                     fastica_fun = "logcosh", fastica_alpha = 1,fastica_method = "R",
                                     fastica_row_norm = FALSE, fastica_maxit = 1000,
                                     fastica_tol = 0.0001, fastica_verbose = FALSE,CES_clustering = TRUE,
                                     consensus_clustering_of_data = FALSE,
                                     distance_function = 'pearson',no_mad_ge = 5000,
                                     consensus_maxK=20,consensus_reps=100,consensus_pItem=0.8,
                                     consensus_pFeature=1,consensus_clusterAlg="hc",
                                     permutation_testing_for_degr_indicator = TRUE,
                                     gene_level_mapping_indicator=data_with_chr_bp_mapping,
                                     intervals_for_gaus_kernel = seq(10000,2000000,by=10000),
                                     probe_no_for_gaus_kernel = 10,
                                     FDR = c(0.01,0.05,0.1,0.2),
                                     CL = 0.5,
                                     state_deciding_cutoff = 0.95,
                                     min_probesets = 5,
                                     set_seed = 12345678,
                                     gene_level_mapping_indicator = TRUE, TACNA_profiling = TRUE)



print(paste("time taken to do the whole analysis is", round((proc.time()[3]-total_time)/60,2),"mins"))

#Load the expression dataset for all samples from CCLE
tab5rows  = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
ccle_expression_dataset <- read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt" , header = TRUE, colClasses = classes)


#Load consensus estimated sources matrix for fold5
ces_fold5 = read.table("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Consensus_Independent_Components_20191030_CCLE_data_fold5_.txt", sep = "\t", header = TRUE)
#Subset the expression dataset with genes from the consensus estimated sources matrix (in this case both are same)
ccle_expression_dataset = ccle_expression_dataset[rownames(ces_fold5),]
#standardize the expression dataset
ccle_expression_dataset_standardized = t(scale(t(ccle_expression_dataset)))

#Obtain mixing matrix weights for samples from CCLE using CES of only the samples from fold5 expression dataset
ces_fold5 = as.matrix(ces_fold5)
dot_prod = t(ces_fold5)%*%ces_fold5
V = solve(dot_prod)
mix_matrix = V%*%t(ces_fold5)%*%as.matrix(ccle_expression_dataset_standardized)
rownames(mix_matrix) = colnames(ces_fold5)

#load extreme valued genomic region information of consensus estimated sources
evr_info = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/SCNA_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_probe_no_for_gaus_kernel_10/Genelevel_20191030_CCLE_data_fold5_All_Components_Amplified_or_Deleted_details_FDR_0.05_CL_0.5_state_deciding_cutoff_0.95_.txt", sep = "\t", header = TRUE))

#Create an indicator matrix with same dimension of the consensus estimated sources matrix having value 1 if a gene falls in an extreme valued genomic region in a consensus estimated source and 0 otherwise
ccle_ica_amp_del_indicator = matrix(0, dim(ces_fold5)[1], dim(ces_fold5)[2])
rownames(ccle_ica_amp_del_indicator) = rownames(ces_fold5)
colnames(ccle_ica_amp_del_indicator) = colnames(ces_fold5)
match_vec = match(evr_info$ENTREZID, rownames(ccle_ica_amp_del_indicator))
for(j in 1:dim(evr_info)[1])
{
  ccle_ica_amp_del_indicator[match_vec[j],evr_info$Component[j]] = abs(evr_info$Amp_or_Del[j])
}

#obtain CNA-CES (copy number alteration consensus estimated sources) using the indicator matrix created above
ccle_cna_ces_fold5 = ces_fold5*ccle_ica_amp_del_indicator

#obtain TACNA profiles of samples from CCLE using the CNA-CESs of samples ony from fold5
ccle_tacna_fold5 = ccle_cna_ces_fold5%*%mix_matrix

#load CN dataset (derived from SNP) for samples from CCLE
ccle_snp_data = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Data/CCLE__median_centered_copy_number_data__genes_ranked_on_genomic_mapping.txt", sep = "\t", header = TRUE)

#convert CN dataset from probeset-level to gene-level
jetset_mapping = as.data.frame(fread("/Users/arkajyotibhattacharya/Projects/Databases/Genomic\ mapping\ files/Genomic_Mapping_hgu133plus2_using_jetscore_30032018.txt", sep = "\t", header = TRUE))
rownames(jetset_mapping) = jetset_mapping$PROBESET
jetset_mapping = jetset_mapping[which(jetset_mapping$PROBESET%in%ccle_snp_data$PROBESET),]
jetset_mapping = jetset_mapping[which(jetset_mapping$top_probe_indicator==1),]
rownames(ccle_snp_data) = ccle_snp_data$PROBESET
ccle_snp_data = ccle_snp_data[rownames(jetset_mapping),]
rownames(ccle_snp_data) = jetset_mapping$ENTREZID


#list of genes common between SNP data and TACNA profiles of CCLE using CNA-CES of fold5
common_genes_snp_tacna = intersect(rownames(ccle_snp_data), rownames(ccle_tacna_fold5))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold5 with common genes
ccle_snp_data = ccle_snp_data[common_genes_snp_tacna,]
ccle_tacna_fold5 = ccle_tacna_fold5[common_genes_snp_tacna,]

common_samples = intersect(colnames(ccle_snp_data),colnames(ccle_tacna_fold5))

#subsetting both SNP data and TACNA profiles of CCLE using CNA-CES of fold5 with common samples
ccle_snp_data = ccle_snp_data[,common_samples]
ccle_tacna_fold5 = ccle_tacna_fold5[,common_samples]

#Obtain correlation between TACNA profiles using CNA-CES of fold5 and SNP profiles CCLE
ccle_snp_data = as.matrix(ccle_snp_data)
cor_between_snp_and_tacna_by_ccle = NULL

for(i in 1:dim(ccle_tacna_fold5)[2])
{
  cor_between_snp_and_tacna_by_ccle[i] = cor(ccle_snp_data[,i], ccle_tacna_fold5[,i], use = "pairwise.complete.obs")
}

#load list of samples used in fold5
samples_in_fold5 = read.table("/Users/arkajyotibhattacharya/Projects/Analysis\ on\ CCLE/Results/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore_fold5.txt", sep = "\t" , header = TRUE, nrows = 5)

#correlation coefficients for samples not belonging to fold5
cor_between_snp_and_tacna_fold5 = cor_between_snp_and_tacna_by_ccle[which(!(common_samples%in%colnames(samples_in_fold5)))]

samples_in_validation_fold5 = common_samples[which(!(common_samples%in%colnames(samples_in_fold5)))]

correlation_fold5 = as.data.frame(cbind(samples_in_validation_fold5, cor_between_snp_and_tacna_fold5))
correlation_fold5$cor_between_snp_and_tacna_fold5 = as.numeric(as.character(correlation_fold5$cor_between_snp_and_tacna_fold5))
#Export file containing correlation coefficients
write.table(correlation_fold5, file = "/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Results/CCLE\ cross\ validation/Correlation between TACNA and SNP fold5.txt", sep = "\t", row.names = FALSE)

#Plot histogram of the above correlation coefficients
pdf("/Users/arkajyotibhattacharya/Projects/TACNA\ resubmission/Plots/Histogram of correlation between SNP and regenerated CCLE TACNA fold5.pdf", height = 6, width = 12)
hist(cor_between_snp_and_tacna_fold5, breaks = 100,main = "Histogram of correlation between SNP and regenerated CCLE TACNA fold5"
     ,xlab = "Pearson correlation coefficient")

dev.off()


