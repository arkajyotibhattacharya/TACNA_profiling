library(data.table)
library(sqldf)
library(Rcpp)

expression_data = read.table("/data/bioinfo-fehrmann/Source codes/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", sep = " ", header = TRUE)

data_with_chr_bp_mapping = read.csv("/data/bioinfo-fehrmann/Source codes/Genomic\ mapping.csv")

rownames(data_with_chr_bp_mapping) = data_with_chr_bp_mapping$ENTREZID
common_genes = intersect(rownames(data_with_chr_bp_mapping), rownames(expression_data))
data_with_chr_bp_mapping = data_with_chr_bp_mapping[common_genes,]
expression_data = expression_data[common_genes,]

data_with_chr_bp_mapping = sqldf("select * from data_with_chr_bp_mapping order by 4,5")
rownames(data_with_chr_bp_mapping) = data_with_chr_bp_mapping$ENTREZID
expression_data = expression_data[rownames(data_with_chr_bp_mapping),1:300]
source("/data/bioinfo-fehrmann/Source codes/Source\ code\ TACNA\ genelevel.R")
setwd("/data/bioinfo-fehrmann/Source codes/")

combat_consensus_ica_degr_clustering(expression_data = expression_data, #input expression data probesets/genes in rows and samples in columns
                                                sample_batch_id = NULL, #Batch information of samples. Column names should be "SAMPLE_IDENTIFIER" and "BATCH_IDENTIFIER"
                                                probelevel_standardization = TRUE, #Logical parameter to indicate if the expression dataset is required to be standardized row-wise before application of the pipeline or not
                                                combat_identifier = FALSE, #Logical parameter to indicate if batch effect should be identified and removed or not.
                                                Title = "Test_dataset",
                                                send_email_indicator = FALSE,#Logical parameter to indicate if email should be sent after each step
                                                email_id,#Provide the email id if the above is true
                                                CES_clustering_algo1 = TRUE,#Logical parameter to indicate if the first algorithm to cluster all the estimated sources should be chosen or not (In this paper we choose it as TRUE)
                                                batch_correction_check =FALSE,#Logical parameter to indicate if it is required to check residual batch effect is significantly present or not
                                                n_princomp_check = TRUE,#Logical parameter to indicate if the input for number of estimated sources from ICA pipeline should be chosen using PCA or not.
                                                prin_comp_cor = FALSE,#Logical parameter to indicate if PCA should be done on correlation matrix or not.
                                                choose_fastICA_ncomp = FALSE,#Logical parameter to indicate if the input for number of estimated sources from ICA pipeline should be manually chosen or not.
                                                var_cutoff_to_choose_ncomp = 0.85,#If PCA is being used to choose the input for number of estimated sources from ICA pipeline, var_cutoff_to_choose_ncomp indicates proportion of variance explained by the targetted number of top principal components
                                                ncomp_without_pca,#If the input for number of estimated sources from ICA pipeline is chosen manually, then how many estimated sources are asked from the pipeline.
                                                nsample_fastica = 25, #Number of ICA runs
                                                seed = 123456,
                                                fastica_alg_typ = "parallel", #if alg.typ == "parallel" the components are extracted simultaneously (the default). if alg.typ == "deflation" the components are extracted one at a time.
                                                fastica_fun = "logcosh", #the functional form of the G function used in the approximation to neg-entropy
                                                fastica_alpha = 1,#constant in range [1, 2] used in approximation to neg-entropy when fun == "logcosh"
                                                fastica_method = "R", #if method == "R" then computations are done exclusively in R. The code allows the interested R user to see exactly what the algorithm does. if method == "C" then C code is used to perform most of the computations, which makes the algorithm run faster. During compilation the C code is linked to an optimized BLAS library if present, otherwise stand-alone BLAS routines are compiled.
                                                fastica_row_norm = FALSE, #a logical value indicating whether rows of the data matrix X should be standardized beforehand
                                                fastica_maxit = 2000, #maximum number of iterations to perform
                                                fastica_tol = 0.0001, #a positive scalar giving the tolerance at which the un-mixing matrix is considered to have converged.
                                                fastica_verbose = TRUE, #a logical value indicating the level of output as the algorithm runs
                                                CES_clustering = TRUE, #a logical value indicating the clustering of CESs from different runs should be done or not
                                                consensus_clustering_of_data = FALSE, #a logical value indicating if consensus clustering should be done on the input expression data or not
                                                distance_function = 'pearson', #distance function to be used. if correlation, then mention the method of correlation.
                                                no_mad_ge = 5000, #Number of top genes (from MAD analysis) to be clustered
                                                consensus_maxK=20, #integer value. maximum cluster number to evaluate.
                                                consensus_reps=100, #integer value. number of subsamples.
                                                consensus_pItem=0.8, #numerical value. proportion of items to sample.
                                                consensus_pFeature=1, #numerical value. proportion of features to sample.
                                                consensus_clusterAlg="hc", #character value. cluster algorithm. 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'kmdist' for k-means upon distance matrices (former km option), or a function that returns a clustering.
                                                permutation_testing_for_degr_indicator = TRUE, #a logical value indicating if extreme valued regions should be determined or not.
                                                gene_level_mapping_indicator = FALSE,#logical parameter to indicate if genelevel summary of DEGR is needed or not
                                                collapse_function = "max",#function to collapse probeset-level data to gene-level
                                                data_with_chr_bp_mapping, #Genomic mapping file with the first column as "PROBESET". Other columns are 'ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping'
                                                intervals_for_gaus_kernel = seq(10000,2000000,by=10000), #choose sequence of base pairs for determining standard deviation parameter of truncated normal distribution
                                                probe_no_for_gaus_kernel = 10, #minimum numbr of probesets/genes to be present in an extreme valued genomic region
                                                FDR = 0.05, #false discovery rate
                                                CL = 0.5, #confidence level
                                                state_deciding_cutoff = 0.95, #cutoff for Secondary indicator marks (sim)
                                                min_probesets = 5, #minimum number of probesets/genes's expression need to be present in the neighbourhood of a gene to be present in an extreme valued genomic region
                                                set_seed = 12345678,
                                                TACNA_profiling = TRUE,#logical parameter to indicate if TACNA-profiles are to be obtained or not
                                                large_data = FALSE,#logical parameter to indicate if the calculation of correlation requires a lot of memory or not.
                                                no_cores_v1 = 5 #number of cores to be used for parallelization
)



# 
# 
# expression_data = expression_data; #input expression data probesets/genes in rows and samples in columns
# sample_batch_id = NULL; #Batch information of samples. Column names should be "SAMPLE_IDENTIFIER" and "BATCH_IDENTIFIER"
# probelevel_standardization = TRUE; #Logical parameter to indicate if the expression dataset is required to be standardized row-wise before application of the pipeline or not
# combat_identifier = FALSE; #Logical parameter to indicate if batch effect should be identified and removed or not.
# Title = "Test_dataset";
# send_email_indicator = FALSE;#Logical parameter to indicate if email should be sent after each step
# email_id;#Provide the email id if the above is true
# CES_clustering_algo1 = TRUE;#Logical parameter to indicate if the first algorithm to cluster all the estimated sources should be chosen or not (In this paper we choose it as TRUE)
# batch_correction_check =FALSE;#Logical parameter to indicate if it is required to check residual batch effect is significantly present or not
# n_princomp_check = TRUE;#Logical parameter to indicate if the input for number of estimated sources from ICA pipeline should be chosen using PCA or not.
# prin_comp_cor = FALSE;#Logical parameter to indicate if PCA should be done on correlation matrix or not.
# choose_fastICA_ncomp = FALSE;#Logical parameter to indicate if the input for number of estimated sources from ICA pipeline should be manually chosen or not.
# var_cutoff_to_choose_ncomp = 0.85;#If PCA is being used to choose the input for number of estimated sources from ICA pipeline, var_cutoff_to_choose_ncomp indicates proportion of variance explained by the targetted number of top principal components
# ncomp_without_pca;#If the input for number of estimated sources from ICA pipeline is chosen manually, then how many estimated sources are asked from the pipeline.
# nsample_fastica = 25; #Number of ICA runs
# seed = 123456;
# fastica_alg_typ = "parallel"; #if alg.typ == "parallel" the components are extracted simultaneously (the default). if alg.typ == "deflation" the components are extracted one at a time.
# fastica_fun = "logcosh"; #the functional form of the G function used in the approximation to neg-entropy
# fastica_alpha = 1;#constant in range [1, 2] used in approximation to neg-entropy when fun == "logcosh"
# fastica_method = "R"; #if method == "R" then computations are done exclusively in R. The code allows the interested R user to see exactly what the algorithm does. if method == "C" then C code is used to perform most of the computations, which makes the algorithm run faster. During compilation the C code is linked to an optimized BLAS library if present, otherwise stand-alone BLAS routines are compiled.
# fastica_row_norm = FALSE; #a logical value indicating whether rows of the data matrix X should be standardized beforehand
# fastica_maxit = 2000; #maximum number of iterations to perform
# fastica_tol = 0.0001; #a positive scalar giving the tolerance at which the un-mixing matrix is considered to have converged.
# fastica_verbose = TRUE; #a logical value indicating the level of output as the algorithm runs
# CES_clustering = TRUE; #a logical value indicating the clustering of CESs from different runs should be done or not
# consensus_clustering_of_data = FALSE; #a logical value indicating if consensus clustering should be done on the input expression data or not
# distance_function = 'pearson'; #distance function to be used. if correlation, then mention the method of correlation.
# no_mad_ge = 5000; #Number of top genes (from MAD analysis) to be clustered 
# consensus_maxK=20; #integer value. maximum cluster number to evaluate.
# consensus_reps=100; #integer value. number of subsamples.
# consensus_pItem=0.8; #numerical value. proportion of items to sample.
# consensus_pFeature=1; #numerical value. proportion of features to sample.
# consensus_clusterAlg="hc"; #character value. cluster algorithm. 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'kmdist' for k-means upon distance matrices (former km option), or a function that returns a clustering.
# permutation_testing_for_degr_indicator = TRUE; #a logical value indicating if extreme valued regions should be determined or not.
# gene_level_mapping_indicator = FALSE;#logical parameter to indicate if genelevel summary of DEGR is needed or not
# collapse_function = "max";#function to collapse probeset-level data to gene-level
# data_with_chr_bp_mapping; #Genomic mapping file with the first column as "PROBESET". Other columns are 'ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping'
# intervals_for_gaus_kernel = seq(10000,2000000,by=10000); #choose sequence of base pairs for determining standard deviation parameter of truncated normal distribution
# probe_no_for_gaus_kernel = 10; #minimum numbr of probesets/genes to be present in an extreme valued genomic region
# FDR = 0.05; #false discovery rate
# CL = 0.5; #confidence level
# state_deciding_cutoff = 0.95; #cutoff for Secondary indicator marks (sim)
# min_probesets = 5; #minimum number of probesets/genes's expression need to be present in the neighbourhood of a gene to be present in an extreme valued genomic region
# set_seed = 12345678;
# TACNA_profiling = TRUE;#logical parameter to indicate if TACNA-profiles are to be obtained or not
# large_data = FALSE;#logical parameter to indicate if the calculation of correlation requires a lot of memory or not.
# no_cores_v1 = 5
# 
# final_avg_source_ic = read.table("/data/bioinfo-fehrmann/Source codes/ICA_0.85_Explained_variance_25_iterations_Probelevel_standardized__parallel/ICA_Consensus_results/Consensus_Estimated_Sources_Test_dataset_.txt",sep = "\t", header = TRUE)
# credibility_index = read.table("/data/bioinfo-fehrmann/Source codes/ICA_0.85_Explained_variance_25_iterations_Probelevel_standardized__parallel/ICA_Consensus_results/Consensus_credibility_indices_Test_dataset_.txt",sep = "\t", header = TRUE)
# matching_vec = match(data_with_chr_bp_mapping$ENTREZID,rownames(final_avg_source_ic))
# final_comps = which(credibility_index>=0.5)
# data_with_chr_bp_mapping = data_with_chr_bp_mapping[!(is.na(matching_vec)),]
# 
# matching_vec = matching_vec[!(is.na(matching_vec))]
# 
# final_avg_source_ic = cbind(as.data.frame(final_avg_source_ic[matching_vec,]),data_with_chr_bp_mapping[,c('ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping')])
# 
# final_avg_source_ic_v1 = final_avg_source_ic
# original_mapping_file = data_with_chr_bp_mapping
# 
# dataset=as.matrix(final_avg_source_ic_v1[,final_comps]);
# full_dataset = final_avg_source_ic
# mix_matrix = read.table("/data/bioinfo-fehrmann/Source codes/ICA_0.85_Explained_variance_25_iterations_Probelevel_standardized__parallel/ICA_Consensus_results/Consensus_Mixing_Matrix_Test_dataset_.txt",sep = "\t", header = TRUE)
# file_main = file.path(getwd(), paste("ICA",
#                                      ifelse(n_princomp_check,paste(var_cutoff_to_choose_ncomp,
#                                                                    "Explained_variance",sep = "_"), paste(ncomp_without_pca,
#                                                                                                           "ncomp_without_pca",sep = "_"))
#                                      ,nsample_fastica,
#                                      "iterations",
#                                      paste(ifelse(probelevel_standardization,"Probelevel_standardized","")),
#                                      paste(ifelse(combat_identifier,"Batch_effect_removed","")),
#                                      "parallel",sep = "_"))
# 
# 
# full_time = proc.time()[3]
# row_num = dim(dataset)[1]
# 
# # file_SCNA = file.path(file_main, paste("SCNA",
# #                                        "FDR",FDR,
# #                                        "CL", CL, "state_deciding_cutoff",state_deciding_cutoff,
# #                                        "probe_no_for_gaus_kernel",probe_no_for_gaus_kernel,
# #                                        "collapse_function", collapse_function,sep = "_"))
# # dir.create(file_SCNA, showWarnings = FALSE)
# # setwd(file_SCNA)
# # 
# # chromosome_seq = c(0,cumsum(table(as.numeric(as.character(data_with_chr_bp_mapping$CHR_Mapping)))))
# # 
# # data_with_chr_bp_mapping$BP_Mapping_v1 = data_with_chr_bp_mapping$BP_Mapping/10
# # 
# # for(i in 2:max(data_with_chr_bp_mapping$CHR_Mapping))
# # {
# #   data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]=data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]+data_with_chr_bp_mapping$BP_Mapping_v1[chromosome_seq[i]]
# # }
# # 
# # label_pos = sqldf("select distinct CHR_Mapping, avg(BP_Mapping_v1) as BP_Mapping from data_with_chr_bp_mapping group by 1 order by 1")
# # 
# # 
# # 
# # ########################
# # #Finding proper interval for sliding Gaussian Kernel
# # #Choosing that interval where 10 or more no. of probesets 
# # #corresponding to that chromosome are there in +/- 3*interval for 95% of the cases
# # ########################
# # quantile_5 = list()
# # for( k in 1:max(data_with_chr_bp_mapping$CHR_Mapping))
# # {
# #   quantile_5[[k]] = array(0,length(intervals_for_gaus_kernel))
# #   for(i in 1:length(intervals_for_gaus_kernel))
# #   {
# #     frequency_of_neigh_probes = NULL
# #     rows = list()
# #     for(j in (chromosome_seq[k]+1):(chromosome_seq[k+1]))
# #     {
# #       frequency_of_neigh_probes[j-chromosome_seq[k]] = length(which(abs(data_with_chr_bp_mapping$BP_Mapping[j] - data_with_chr_bp_mapping$BP_Mapping[c((chromosome_seq[k]+1):(chromosome_seq[k+1]))]) < 3*intervals_for_gaus_kernel[i]))
# #       
# #     }
# #     quantile_5[[k]][i] = quantile(frequency_of_neigh_probes,0.05)%/%1
# #     
# #     if(quantile_5[[k]][i]>=probe_no_for_gaus_kernel){
# #       print(paste("found",probe_no_for_gaus_kernel, "or more number of probe sets in 95% of the chromosome at interval", 
# #                   intervals_for_gaus_kernel[i],"for Chromosome",k))
# #       quantile_5[[k]][which(quantile_5[[k]]==0)] = NA
# #       break
# #     } 
# #   }
# #   
# # }
# # 
# # 
# # ###########################
# # #Creating the Density Matrix to get sliding gaussian Kernel
# # ###########################
# # 
# # density_matrix = matrix(0, row_num, row_num)
# # 
# # for(j in 1:row_num)
# # {
# #   rows[[j]] = which(data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[j])
# #   k = data_with_chr_bp_mapping$CHR_Mapping[j]
# #   s_x = as.numeric(data_with_chr_bp_mapping$BP_Mapping[rows[[j]]])
# #   s_mean = as.numeric(data_with_chr_bp_mapping$BP_Mapping[j]) 
# #   
# #   den_tnorm = dtruncnorm(s_x,
# #                          a = s_x[1],
# #                          b = s_x[length(s_x)],
# #                          mean = s_mean, 
# #                          sd = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[k]])))])
# #   
# #   density_matrix[j,rows[[j]]] = den_tnorm/sum(den_tnorm)
# # }
# # 
# # 
# # ###########################
# # #Assigning intervals to a vector for each chromosome
# # ###########################
# # 
# # final_intervals = NULL
# # 
# # for(chr_no in 1:max(data_with_chr_bp_mapping$CHR_Mapping)){
# #   final_intervals[chr_no] = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[chr_no]])))]
# # }
# # 
# # 
# # ################
# # #Creating different directories for different outputs
# # ################
# # file_probelevel_bp = file.path(getwd(), "DEGR_sorted_by_base_pair_number_probelevel")
# # file_probelevel_value = file.path(getwd(), "DEGR_sorted_by_CES_value_probelevel")
# # dir.create(file_probelevel_bp, showWarnings = FALSE)
# # dir.create(file_probelevel_value, showWarnings = FALSE)
# # if(gene_level_mapping_indicator)
# # {
# #   file_genelevel_bp = file.path(getwd(), "DEGR_sorted_by_base_pair_number_genelevel")
# #   file_genelevel_value = file.path(getwd(), "DEGR_sorted_by_CES_value_genelevel")
# #   dir.create(file_genelevel_bp, showWarnings = FALSE)
# #   dir.create(file_genelevel_value, showWarnings = FALSE)  
# # }
# # 
# # sourceCpp("/data/bioinfo-fehrmann/SCNA project/Codes/CCLE/Matrix_mult.cpp")
# # 
# # 
# # extreme_valued_probesets_summary = NULL
# # extreme_valued_section_summary = NULL
# # for(CES_no in 1:10)
# # {
# #   set.seed(set_seed+CES_no)
# #   time1 = proc.time()[3]
# #   # 1000 permutation of the non-smoothened weights of the gene-level CES are retained in a matrix p_CESMpx1000.
# #   N <- cbind(dataset[,CES_no],replicate(1000, sample(dataset[,CES_no])))
# #   # Obtain smoothened permuted CESM (sp_CESM) using the smoothing method explained in the step b.
# #   permute_1 = eigenMatMult(density_matrix, N)
# #   # permute_1 = density_matrix%*%N
# #   # Sort the absolute values of the weights of all the columns of sp_CESMpx1000 in decreasing order. Also 3.	Sort the absolute values of the weights of the smoothened CES in decreasing order (sorted_CES).
# #   permute_1_v1 = apply(abs(permute_1),2,function(x) sort(x,decreasing=TRUE))
# #   
# #   # For every column of sp_CESMpx1000 (sp_CESi),
# #   # i.	For every weight of sorted_CES (sorted_CESj), obtain the number of weights of sp_CESi greater than sorted_CESj (f>j)
# #   # ii.	Obtain the optimal cutoff for the sp_CESi (oci) as the maximum value of the weights of sorted_CES for which f>j/j > 5%.
# #   
# #   cutoffs = array(0,1000)
# #   for(j in 2:1001)
# #   {
# #     for(i in 1:dim(permute_1_v1)[1])
# #     {
# #       if(length(which(permute_1_v1[,j]>permute_1_v1[i,1]))/i>=FDR)
# #       {
# #         cutoffs[j-1] = permute_1_v1[i,1]
# #         break
# #       }
# #     }
# #   }
# #   
# #   
# #   # Initial indicator marks (iim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way
# #   # i.	If smoothened_CESs > median(oc) then iims = 1
# #   # ii.	If smoothened_CESs < -median(oc) then iims = -1
# #   # iii.	Otherwise zero.
# #   
# #   indicator_ica_2 = ifelse(permute_1[,1]>quantile(cutoffs,CL),1,ifelse(permute_1[,1]< -quantile(cutoffs,CL),-1,0))
# #   
# #   # Smoothen iim (smoothened_iim) using the smoothing method described in step b.
# #   indicator_ica_3 = density_matrix%*%indicator_ica_2
# #   
# #   # Secondary indicator marks (sim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way:
# #   # i.	If smoothened_iims > 0.85 then sims = 1
# #   # ii.	If smoothened_iims < -0.85 then sims = -1
# #   # iii.	Otherwise zero.
# #   
# #   indicator_ica_4 = ifelse(indicator_ica_3>state_deciding_cutoff,1, ifelse(indicator_ica_3< -state_deciding_cutoff,-1,0))
# #   
# #   # Final indicator marks (fim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way
# #   # Obtain the number of genes (ngs) mapped to the corresponding chromosome which have a distance from gene s in terms of base pair number < corresponding optimal interval length (oil) as described in step b.3.ii.
# #   # If ngs < 10 then fims = 0
# #   # Otherwise, fims = sims
# #   
# #   if(min_probesets>0)
# #   {
# #     state_of_the_weight_pos = which(abs(indicator_ica_4)>0)
# #     
# #     for(no in 1:length(state_of_the_weight_pos)){
# #       if(length(intersect(which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]] >= data_with_chr_bp_mapping$BP_Mapping[state_of_the_weight_pos[no]]-3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]]),
# #                           which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]] <= data_with_chr_bp_mapping$BP_Mapping[state_of_the_weight_pos[no]]+3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]])))<min_probesets){
# #         indicator_ica_4[state_of_the_weight_pos[no]] = 0
# #       }
# #     }
# #     
# #   }
# #   indicator_run = rle(as.vector(indicator_ica_4)) 
# #   
# #   limit = max(abs(min(dataset[,CES_no])),abs(max(dataset[,CES_no]))) 
# #   
# #   par(mar=c(3.1,3.1,3.1,3.1))
# #   plot(data_with_chr_bp_mapping$BP_Mapping_v1,
# #        dataset[,CES_no],
# #        pch=20,cex=0.4,ylim = c(-limit,limit),xaxt = "n",
# #        xlab = NA,
# #        ylab = NA, main = paste("CES", colnames(dataset)[CES_no], "FDR", FDR,  "CL", CL, "state_deciding_cutoff",state_deciding_cutoff), col=rgb(0, 0, 0, 0.2))
# #   mtext(side = 2, line = 2, 'CES Weight',cex = 0.8)
# #   mtext(side = 1, line = 2, paste("Chromosome No.s","of CES", colnames(dataset)[CES_no]),cex = 0.8)
# #   
# #   par(new = T)
# #   plot(data_with_chr_bp_mapping$BP_Mapping_v1,
# #        indicator_ica_4,xaxt = "n",
# #        pch=20,cex=0.4,
# #        ylim = c(-1,
# #                 1),
# #        xlab = NA,
# #        ylab = NA,axes = F, col=rgb(1, 0, 0, 0.2))
# #   axis(side = 4)
# #   mtext(side = 4, line = 2, 'State of the weight',cex = 0.8)
# #   axis(side = 1, at = label_pos[,2],labels = label_pos[,1])
# #   
# #   
# #   
# #   
# #   ######################
# #   #getting the extreme valued genomic positions
# #   ######################
# #   which_probe_extreme_valued = which(abs(indicator_ica_4)>0)
# #   if(length(which_probe_extreme_valued)>0){
# #     extreme_valued_probesets = as.data.frame(cbind(
# #       dataset[which_probe_extreme_valued,CES_no],
# #       row.names(dataset)[which_probe_extreme_valued],
# #       unlist(lapply(data_with_chr_bp_mapping$GENETITLE[which_probe_extreme_valued],as.character)),
# #       unlist(lapply(data_with_chr_bp_mapping$SYMBOL[which_probe_extreme_valued],as.character)),
# #       data_with_chr_bp_mapping$CHR_Mapping[which_probe_extreme_valued],
# #       data_with_chr_bp_mapping$BP_Mapping[which_probe_extreme_valued],
# #       data_with_chr_bp_mapping$BP_Mapping_v1[which_probe_extreme_valued],
# #       
# #       indicator_ica_4[which_probe_extreme_valued],
# #       
# #       rep(CES_no,length(which_probe_extreme_valued))))
# #     
# #     colnames(extreme_valued_probesets) = c( "Value", "ENTREZID", 'GENETITLE', 'SYMBOL',"CHR_Mapping", "BP_Mapping", "BP_Mapping_for_plot", "state_of_the_weight","CES")
# #     
# #     extreme_valued_probesets_summary = rbind(extreme_valued_probesets_summary,extreme_valued_probesets)
# #     
# #     run_vector = c(0, cumsum(indicator_run$lengths))
# #     
# #     extreme_valued_runs = which(abs(indicator_run$values)>0)
# #     
# #     for(num in 1:length(extreme_valued_runs)){
# #       
# #       if(gene_level_mapping_indicator){
# #         rows_vec = which(original_mapping_file$ENTREZID%in%data_with_chr_bp_mapping$ENTREZID[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1]])
# #         extreme_valued_section_gene_level = as.data.frame(cbind(dataset[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1],colnames(dataset)[CES_no]],data_with_chr_bp_mapping[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1],])) 
# #         
# #         colnames(extreme_valued_section_gene_level)[1] = "Value"
# #         
# #         if(indicator_run$values[extreme_valued_runs[num]]>0)
# #         {
# #           extreme_valued_section_gene_level_sort_by_value = sqldf("select * from extreme_valued_section_gene_level order by value DESC")
# #           
# #         }else{
# #           extreme_valued_section_gene_level_sort_by_value = sqldf("select * from extreme_valued_section_gene_level order by value")
# #           
# #         }
# #         
# #         write.table(extreme_valued_section_gene_level,
# #                     file = paste(file_genelevel_bp,paste(Title,"_CHR", extreme_valued_section_gene_level$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_gene_level$BP_Mapping[1]), "_BPEND",
# #                                                          round(extreme_valued_section_gene_level$BP_Mapping[dim(extreme_valued_section_gene_level)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
# #                     sep = "\t", row.names = FALSE)
# #         
# #         write.table(extreme_valued_section_gene_level_sort_by_value,
# #                     file = paste(file_genelevel_value,paste(Title,"_CHR", extreme_valued_section_gene_level_sort_by_value$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_gene_level_sort_by_value$BP_Mapping[1]), "_BPEND",
# #                                                             round(extreme_valued_section_gene_level_sort_by_value$BP_Mapping[dim(extreme_valued_section_gene_level_sort_by_value)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
# #                     sep = "\t", row.names = FALSE)
# #         
# #       }else{
# #         rows_vec = which(original_mapping_file$PROBESET%in%data_with_chr_bp_mapping$PROBESET[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1]])
# #         
# #       }
# #       
# #       
# #       
# #       
# #       extreme_valued_section = as.data.frame(full_dataset[rows_vec,c(colnames(dataset)[CES_no],'ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping')])
# #       
# #       extreme_valued_section$PROBESET = as.character(row.names(extreme_valued_section))
# #       
# #       
# #       colnames(extreme_valued_section)[1] = "Value"
# #       
# #       if(indicator_run$values[extreme_valued_runs[num]]>0)
# #       {
# #         extreme_valued_section_sort_by_value = sqldf("select * from extreme_valued_section order by value DESC")
# #       }else{
# #         extreme_valued_section_sort_by_value = sqldf("select * from extreme_valued_section order by value")
# #       }
# #       
# #       write.table(extreme_valued_section,
# #                   file = paste(file_probelevel_bp,paste(Title,"_CHR", extreme_valued_section$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section$BP_Mapping[1]), "_BPEND",
# #                                                         round(extreme_valued_section$BP_Mapping[dim(extreme_valued_section)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
# #                   sep = "\t", row.names = FALSE)
# #       
# #       write.table(extreme_valued_section_sort_by_value,
# #                   file = paste(file_probelevel_value,paste(Title,"_CHR", extreme_valued_section_sort_by_value$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_sort_by_value$BP_Mapping[1]), "_BPEND",
# #                                                            round(extreme_valued_section_sort_by_value$BP_Mapping[dim(extreme_valued_section_sort_by_value)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
# #                   sep = "\t", row.names = FALSE)
# #       CES = rep(CES_no,dim(extreme_valued_section)[1])
# #       
# #       extreme_valued_section_summary = rbind(extreme_valued_section_summary,cbind(extreme_valued_section,CES))
# #       
# #     }
# #     
# #     
# #     print(paste("Extreme valued identified in some genomic postiions of CES",colnames(dataset)[CES_no]))
# #   }else{
# #     print(paste("No Genomic position has got Extreme valued region for CES",colnames(dataset)[CES_no]))
# #   }
# #   print(paste("Processing of CES", colnames(dataset)[CES_no], "FDR", FDR,  "CL", CL, "state_deciding_cutoff",state_deciding_cutoff, "is done"))
# #   
# #   print(paste("Time taken for this iteration is",round((proc.time()[3]-time1)/60,3),"mins"))
# #   
# # }
# # 
# 
# 
# #############try
# 
# ##########################
# #Adjustment for Plot
# ###########################
# 
# file_SCNA = file.path(file_main, paste("SCNA",
#                                        "FDR",FDR,
#                                        "CL", CL, "state_deciding_cutoff",state_deciding_cutoff,
#                                        "probe_no_for_gaus_kernel",probe_no_for_gaus_kernel,
#                                        "collapse_function", collapse_function,sep = "_"))
# dir.create(file_SCNA, showWarnings = FALSE)
# setwd(file_SCNA)
# 
# chromosome_seq = c(0,cumsum(table(as.numeric(as.character(data_with_chr_bp_mapping$CHR_Mapping)))))
# 
# data_with_chr_bp_mapping$BP_Mapping_v1 = data_with_chr_bp_mapping$BP_Mapping/10
# 
# for(i in 2:max(data_with_chr_bp_mapping$CHR_Mapping))
# {
#   data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]=data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]+data_with_chr_bp_mapping$BP_Mapping_v1[chromosome_seq[i]]
# }
# 
# label_pos = sqldf("select distinct CHR_Mapping, avg(BP_Mapping_v1) as BP_Mapping from data_with_chr_bp_mapping group by 1 order by 1")
# 
# 
# 
# ########################
# #Finding proper interval for sliding Gaussian Kernel
# #Choosing that interval where 10 or more no. of probesets 
# #corresponding to that chromosome are there in +/- 3*interval for 95% of the cases
# ########################
# quantile_5 = list()
# for( k in 1:max(data_with_chr_bp_mapping$CHR_Mapping))
# {
#   quantile_5[[k]] = array(0,length(intervals_for_gaus_kernel))
#   for(i in 1:length(intervals_for_gaus_kernel))
#   {
#     frequency_of_neigh_probes = NULL
#     rows = list()
#     for(j in (chromosome_seq[k]+1):(chromosome_seq[k+1]))
#     {
#       frequency_of_neigh_probes[j-chromosome_seq[k]] = length(which(abs(data_with_chr_bp_mapping$BP_Mapping[j] - data_with_chr_bp_mapping$BP_Mapping[c((chromosome_seq[k]+1):(chromosome_seq[k+1]))]) < 3*intervals_for_gaus_kernel[i]))
#       
#     }
#     quantile_5[[k]][i] = quantile(frequency_of_neigh_probes,0.05)%/%1
#     
#     if(quantile_5[[k]][i]>=probe_no_for_gaus_kernel){
#       print(paste("found",probe_no_for_gaus_kernel, "or more number of probe sets in 95% of the chromosome at interval", 
#                   intervals_for_gaus_kernel[i],"for Chromosome",k))
#       quantile_5[[k]][which(quantile_5[[k]]==0)] = NA
#       break
#     } 
#   }
#   
# }
# 
# 
# ###########################
# #Creating the Density Matrix to get sliding gaussian Kernel
# ###########################
# 
# density_matrix = matrix(0, row_num, row_num)
# 
# for(j in 1:row_num)
# {
#   rows[[j]] = which(data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[j])
#   k = data_with_chr_bp_mapping$CHR_Mapping[j]
#   s_x = as.numeric(data_with_chr_bp_mapping$BP_Mapping[rows[[j]]])
#   s_mean = as.numeric(data_with_chr_bp_mapping$BP_Mapping[j]) 
#   
#   den_tnorm = dtruncnorm(s_x,
#                          a = s_x[1],
#                          b = s_x[length(s_x)],
#                          mean = s_mean, 
#                          sd = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[k]])))])
#   
#   density_matrix[j,rows[[j]]] = den_tnorm/sum(den_tnorm)
# }
# 
# 
# ###########################
# #Assigning intervals to a vector for each chromosome
# ###########################
# 
# final_intervals = NULL
# 
# for(chr_no in 1:max(data_with_chr_bp_mapping$CHR_Mapping)){
#   final_intervals[chr_no] = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[chr_no]])))]
# }
# 
# 
# ###########################
# #Permutation Test for obtaining extreme valued region
# ###########################
# 
# ################
# #Creating different directories for different outputs
# ################
# 
# file_genelevel_bp = file.path(getwd(), "DEGR_sorted_by_base_pair_number_genelevel")
# file_genelevel_value = file.path(getwd(), "DEGR_sorted_by_CES_value_genelevel")
# dir.create(file_genelevel_bp, showWarnings = FALSE)
# dir.create(file_genelevel_value, showWarnings = FALSE)  
# 
# 
# 
# 
# # pdf(paste(Title,"All_Independent_CESs_plot.pdf",sep = "_"), height = 20,width= 30)
# extreme_valued_probesets_summary = NULL
# extreme_valued_section_summary = NULL
# for(CES_no in 1:2)
# {
#   set.seed(set_seed+CES_no)
#   time1 = proc.time()[3]
#   # 1000 permutation of the non-smoothened weights of the gene-level CES are retained in a matrix p_CESMpx1000.
#   N <- cbind(dataset[,CES_no],replicate(1000, sample(dataset[,CES_no])))
#   # Obtain smoothened permuted CESM (sp_CESM) using the smoothing method explained in the step b.
#   permute_1 = density_matrix%*%N
#   # Sort the absolute values of the weights of all the columns of sp_CESMpx1000 in decreasing order. Also 3.	Sort the absolute values of the weights of the smoothened CES in decreasing order (sorted_CES).
#   permute_1_v1 = apply(abs(permute_1),2,function(x) sort(x,decreasing=TRUE))
#   
#   # For every column of sp_CESMpx1000 (sp_CESi),
#   # i.	For every weight of sorted_CES (sorted_CESj), obtain the number of weights of sp_CESi greater than sorted_CESj (f>j)
#   # ii.	Obtain the optimal cutoff for the sp_CESi (oci) as the maximum value of the weights of sorted_CES for which f>j/j > 5%.
#   
#   cutoffs = array(0,1000)
#   for(j in 2:1001)
#   {
#     for(i in 1:dim(permute_1_v1)[1])
#     {
#       if(length(which(permute_1_v1[,j]>permute_1_v1[i,1]))/i>=FDR)
#       {
#         cutoffs[j-1] = permute_1_v1[i,1]
#         break
#       }
#     }
#   }
#   
#   
#   # Initial indicator marks (iim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way
#   # i.	If smoothened_CESs > median(oc) then iims = 1
#   # ii.	If smoothened_CESs < -median(oc) then iims = -1
#   # iii.	Otherwise zero.
#   
#   indicator_ica_2 = ifelse(permute_1[,1]>quantile(cutoffs,CL),1,ifelse(permute_1[,1]< -quantile(cutoffs,CL),-1,0))
#   
#   # Smoothen iim (smoothened_iim) using the smoothing method described in step b.
#   indicator_ica_3 = density_matrix%*%indicator_ica_2
#   
#   # Secondary indicator marks (sim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way:
#   # i.	If smoothened_iims > 0.85 then sims = 1
#   # ii.	If smoothened_iims < -0.85 then sims = -1
#   # iii.	Otherwise zero.
#   
#   indicator_ica_4 = ifelse(indicator_ica_3>state_deciding_cutoff,1, ifelse(indicator_ica_3< -state_deciding_cutoff,-1,0))
#   
#   # Final indicator marks (fim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way
#   # Obtain the number of genes (ngs) mapped to the corresponding chromosome which have a distance from gene s in terms of base pair number < corresponding optimal interval length (oil) as described in step b.3.ii.
#   # If ngs < 10 then fims = 0
#   # Otherwise, fims = sims
#   
#   if(min_probesets>0)
#   {
#     state_of_the_weight_pos = which(abs(indicator_ica_4)>0)
#     
#     for(no in 1:length(state_of_the_weight_pos)){
#       if(length(intersect(which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]] >= data_with_chr_bp_mapping$BP_Mapping[state_of_the_weight_pos[no]]-3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]]),
#                           which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]] <= data_with_chr_bp_mapping$BP_Mapping[state_of_the_weight_pos[no]]+3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]])))<min_probesets){
#         indicator_ica_4[state_of_the_weight_pos[no]] = 0
#       }
#     }
#     
#   }
#   indicator_run = rle(as.vector(indicator_ica_4)) 
#   
#   limit = max(abs(min(dataset[,CES_no])),abs(max(dataset[,CES_no]))) 
#   
#   # par(mar=c(3.1,3.1,3.1,3.1))
#   # plot(data_with_chr_bp_mapping$BP_Mapping_v1,
#   #      dataset[,CES_no],
#   #      pch=20,cex=0.4,ylim = c(-limit,limit),xaxt = "n",
#   #      xlab = NA,
#   #      ylab = NA, main = paste("CES", colnames(dataset)[CES_no], "FDR", FDR,  "CL", CL, "state_deciding_cutoff",state_deciding_cutoff), col=rgb(0, 0, 0, 0.2))
#   # mtext(side = 2, line = 2, 'CES Weight',cex = 0.8)
#   # mtext(side = 1, line = 2, paste("Chromosome No.s","of CES", colnames(dataset)[CES_no]),cex = 0.8)
#   # 
#   # par(new = T)
#   # plot(data_with_chr_bp_mapping$BP_Mapping_v1,
#   #      indicator_ica_4,xaxt = "n",
#   #      pch=20,cex=0.4,
#   #      ylim = c(-1,
#   #               1),
#   #      xlab = NA,
#   #      ylab = NA,axes = F, col=rgb(1, 0, 0, 0.2))
#   # axis(side = 4)
#   # mtext(side = 4, line = 2, 'State of the weight',cex = 0.8)
#   # axis(side = 1, at = label_pos[,2],labels = label_pos[,1])
#   # 
#   
#   
#   
#   ######################
#   #getting the extreme valued genomic positions
#   ######################
#   which_probe_extreme_valued = which(abs(indicator_ica_4)>0)
#   if(length(which_probe_extreme_valued)>0){
#     extreme_valued_probesets = as.data.frame(cbind(
#       dataset[which_probe_extreme_valued,CES_no],
#       row.names(dataset)[which_probe_extreme_valued],
#       unlist(lapply(data_with_chr_bp_mapping$GENETITLE[which_probe_extreme_valued],as.character)),
#       unlist(lapply(data_with_chr_bp_mapping$SYMBOL[which_probe_extreme_valued],as.character)),
#       data_with_chr_bp_mapping$CHR_Mapping[which_probe_extreme_valued],
#       data_with_chr_bp_mapping$BP_Mapping[which_probe_extreme_valued],
#       data_with_chr_bp_mapping$BP_Mapping_v1[which_probe_extreme_valued],
#       
#       indicator_ica_4[which_probe_extreme_valued],
#       
#       rep(CES_no,length(which_probe_extreme_valued))))
#     
#     colnames(extreme_valued_probesets) = c( "Value", "ENTREZID", 'GENETITLE', 'SYMBOL',"CHR_Mapping", "BP_Mapping", "BP_Mapping_for_plot", "state_of_the_weight","CES")
#     
#     extreme_valued_probesets_summary = rbind(extreme_valued_probesets_summary,extreme_valued_probesets)
#     
#     run_vector = c(0, cumsum(indicator_run$lengths))
#     
#     extreme_valued_runs = which(abs(indicator_run$values)>0)
#     
#     for(num in 1:length(extreme_valued_runs)){
#       
#       
#       rows_vec = which(original_mapping_file$ENTREZID%in%data_with_chr_bp_mapping$ENTREZID[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1]])
#       extreme_valued_section_gene_level = as.data.frame(cbind(dataset[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1],colnames(dataset)[CES_no]],data_with_chr_bp_mapping[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1],])) 
#       
#       colnames(extreme_valued_section_gene_level)[1] = "Value"
#       
#       if(indicator_run$values[extreme_valued_runs[num]]>0)
#       {
#         extreme_valued_section_gene_level_sort_by_value = sqldf("select * from extreme_valued_section_gene_level order by value DESC")
#         
#       }else{
#         extreme_valued_section_gene_level_sort_by_value = sqldf("select * from extreme_valued_section_gene_level order by value")
#         
#       }
#       
#       # write.table(extreme_valued_section_gene_level,
#       #             file = paste(file_genelevel_bp,paste(Title,"_CHR", extreme_valued_section_gene_level$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_gene_level$BP_Mapping[1]), "_BPEND",
#       #                                                  round(extreme_valued_section_gene_level$BP_Mapping[dim(extreme_valued_section_gene_level)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
#       #             sep = "\t", row.names = FALSE)
#       # 
#       # write.table(extreme_valued_section_gene_level_sort_by_value,
#       #             file = paste(file_genelevel_value,paste(Title,"_CHR", extreme_valued_section_gene_level_sort_by_value$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_gene_level_sort_by_value$BP_Mapping[1]), "_BPEND",
#       #                                                     round(extreme_valued_section_gene_level_sort_by_value$BP_Mapping[dim(extreme_valued_section_gene_level_sort_by_value)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
#       #             sep = "\t", row.names = FALSE)
#       # 
#       
#       
#       
#       
#       extreme_valued_section = as.data.frame(full_dataset[rows_vec,c(colnames(dataset)[CES_no],'ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping')])
#       
#       colnames(extreme_valued_section)[1] = "Value"
#       
#       if(indicator_run$values[extreme_valued_runs[num]]>0)
#       {
#         extreme_valued_section_sort_by_value = sqldf("select * from extreme_valued_section order by value DESC")
#       }else{
#         extreme_valued_section_sort_by_value = sqldf("select * from extreme_valued_section order by value")
#       }
#       
#       # write.table(extreme_valued_section,
#       #             file = paste(file_probelevel_bp,paste(Title,"_CHR", extreme_valued_section$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section$BP_Mapping[1]), "_BPEND",
#       #                                                   round(extreme_valued_section$BP_Mapping[dim(extreme_valued_section)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
#       #             sep = "\t", row.names = FALSE)
#       # 
#       # write.table(extreme_valued_section_sort_by_value,
#       #             file = paste(file_probelevel_value,paste(Title,"_CHR", extreme_valued_section_sort_by_value$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_sort_by_value$BP_Mapping[1]), "_BPEND",
#       #                                                      round(extreme_valued_section_sort_by_value$BP_Mapping[dim(extreme_valued_section_sort_by_value)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
#       #             sep = "\t", row.names = FALSE)
#       CES = rep(CES_no,dim(extreme_valued_section)[1])
#       
#       extreme_valued_section_summary = rbind(extreme_valued_section_summary,cbind(extreme_valued_section,CES))
#       
#     }
#     
#     
#     print(paste("Extreme valued identified in some genomic postiions of CES",colnames(dataset)[CES_no]))
#   }else{
#     print(paste("No Genomic position has got Extreme valued region for CES",colnames(dataset)[CES_no]))
#   }
#   print(paste("Processing of CES", colnames(dataset)[CES_no], "FDR", FDR,  "CL", CL, "state_deciding_cutoff",state_deciding_cutoff, "is done"))
#   
#   print(paste("Time taken for this iteration is",round((proc.time()[3]-time1)/60,3),"mins"))
#   
# } 
