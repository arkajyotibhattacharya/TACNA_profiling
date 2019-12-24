library(data.table)
library(sqldf)
expression_data = read.table("/Applications/AMPPS/www/data/DownloadData/Source\ codes/CCLE__Affy_hgu133plus2_QCed_mRNA_NoDuplicates_CleanedIdentifiers_RMA-sketch_genelevel_using_jetscore.txt", sep = " ", header = TRUE)

data_with_chr_bp_mapping = read.csv("/Applications/AMPPS/www/data/DownloadData/Source\ codes/Genomic\ mapping.csv")

rownames(data_with_chr_bp_mapping) = data_with_chr_bp_mapping$ENTREZID
common_genes = intersect(rownames(data_with_chr_bp_mapping), rownames(expression_data))
data_with_chr_bp_mapping = data_with_chr_bp_mapping[common_genes,]
expression_data = expression_data[common_genes,]

data_with_chr_bp_mapping = sqldf("select * from data_with_chr_bp_mapping order by 4,5")
rownames(data_with_chr_bp_mapping) = data_with_chr_bp_mapping$ENTREZID
expression_data = expression_data[rownames(data_with_chr_bp_mapping),1:100]
source("/Applications/AMPPS/www/data/DownloadData/Source\ codes/Source\ code\ TACNA.R")
setwd("/Applications/AMPPS/www/data/DownloadData/Source\ codes/")

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
