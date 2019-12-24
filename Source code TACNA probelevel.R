###################
#Loading Libraries
###################
library(sva)
library(fastICA)
library(evd)
library(rARPACK)
library(DBI)
library(sqldf)
library(caret)
library(gPCA)
library(ConsensusClusterPlus)
library(foreach)
library(doMC)
library(biganalytics)
library(clValid)
library(amap)
library(sendmailR)
library(truncnorm)
library(parallel)


combat_consensus_ica_degr_clustering = function(expression_data, #input expression data probesets/genes in rows and samples in columns
                                           sample_batch_id, #Batch information of samples. Column names should be "SAMPLE_IDENTIFIER" and "BATCH_IDENTIFIER"
                                           probelevel_standardization = TRUE, #Logical parameter to indicate if the expression dataset is required to be standardized row-wise before application of the pipeline or not
                                           combat_identifier = TRUE, #Logical parameter to indicate if batch effect should be identified and removed or not.
                                           Title,
                                           send_email_indicator = TRUE,#Logical parameter to indicate if email should be sent after each step
                                           email_id,#Provide the email id if the above is true
                                           CES_clustering_algo1 = TRUE,#Logical parameter to indicate if the first algorithm to cluster all the estimated sources should be chosen or not (In this paper we choose it as TRUE)
                                           batch_correction_check =TRUE,#Logical parameter to indicate if it is required to check residual batch effect is significantly present or not
                                           n_princomp_check = TRUE,#Logical parameter to indicate if the input for number of estimated sources from ICA pipeline should be chosen using PCA or not.
                                           prin_comp_cor = FALSE,#Logical parameter to indicate if PCA should be done on correlation matrix or not.
                                           choose_fastICA_ncomp = FALSE,#Logical parameter to indicate if the input for number of estimated sources from ICA pipeline should be manually chosen or not.
                                           var_cutoff_to_choose_ncomp = 0.99,#If PCA is being used to choose the input for number of estimated sources from ICA pipeline, var_cutoff_to_choose_ncomp indicates proportion of variance explained by the targetted number of top principal components
                                           ncomp_without_pca,#If the input for number of estimated sources from ICA pipeline is chosen manually, then how many estimated sources are asked from the pipeline.
                                           nsample_fastica = 25, #Number of ICA runs
                                           seed = 123456,
                                           fastica_alg_typ = "parallel", #if alg.typ == "parallel" the components are extracted simultaneously (the default). if alg.typ == "deflation" the components are extracted one at a time.
                                           fastica_fun = "logcosh", #the functional form of the G function used in the approximation to neg-entropy
                                           fastica_alpha = 1,#constant in range [1, 2] used in approximation to neg-entropy when fun == "logcosh"
                                           fastica_method = "C", #if method == "R" then computations are done exclusively in R. The code allows the interested R user to see exactly what the algorithm does. if method == "C" then C code is used to perform most of the computations, which makes the algorithm run faster. During compilation the C code is linked to an optimized BLAS library if present, otherwise stand-alone BLAS routines are compiled.
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
                                           gene_level_mapping_indicator = TRUE,#logical parameter to indicate if genelevel summary of DEGR is needed or not
                                           collapse_function = "max",#function to collapse probeset-level data to gene-level
                                           data_with_chr_bp_mapping, #Genomic mapping file with the first column as "PROBESET". Other columns are 'ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping'
                                           intervals_for_gaus_kernel = seq(10000,1000000,by=10000), #choose sequence of base pairs for determining standard deviation parameter of truncated normal distribution
                                           probe_no_for_gaus_kernel = 10, #minimum numbr of probesets/genes to be present in an extreme valued genomic region
                                           FDR = 0.1, #false discovery rate
                                           CL = 0.5, #confidence level
                                           state_deciding_cutoff = 0.95, #cutoff for Secondary indicator marks (sim)
                                           min_probesets = 5, #minimum number of probesets/genes's expression need to be present in the neighbourhood of a gene to be present in an extreme valued genomic region
                                           set_seed = 12345678,
                                           TACNA_profiling = TRUE,#logical parameter to indicate if TACNA-profiles are to be obtained or not
                                           large_data = FALSE,#logical parameter to indicate if the calculation of correlation requires a lot of memory or not.
                                           no_cores_v1 = detectCores() #number of cores to be used for parallelization
)
{
  if((choose_fastICA_ncomp)&&(var_cutoff_to_choose_ncomp!=0.99))
  {
    print("User either chooses no. of CESs by his/her own or selects the variance cut-off. Please keep var_cutoff_to_choose_ncomp as 0.99 if you want to choose number of CESs by your own")
  }else{
    
    file_main = file.path(getwd(), paste("ICA",
                                         ifelse(n_princomp_check,paste(var_cutoff_to_choose_ncomp,
                                                                       "Explained_variance",sep = "_"), paste(ncomp_without_pca,
                                                                                                              "ncomp_without_pca",sep = "_"))
                                         ,nsample_fastica,
                                         "iterations",
                                         paste(ifelse(probelevel_standardization,"Probelevel_standardized","")),
                                         paste(ifelse(combat_identifier,"Batch_effect_removed","")),
                                         "parallel",sep = "_"))
    dir.create(file_main, showWarnings = FALSE)
    
    setwd(file_main)
    
    if(probelevel_standardization){
      expression_data = t(scale(t(expression_data)))
      print("Probelevel standardization is done")
    }
    if(combat_identifier)
    {
      file_batch_effect = file.path(getwd(), "Batch_Effect_Correction")
      dir.create(file_batch_effect, showWarnings = FALSE)
      
      sample_batch_id$SAMPLE_IDENTIFIER = as.character(sample_batch_id$SAMPLE_IDENTIFIER)
      sample_batch_id$BATCH_IDENTIFIER = as.character(sample_batch_id$BATCH_IDENTIFIER)
      if(sum(colnames(expression_data)==as.character(sample_batch_id[,1]))==dim(expression_data)[2])
      {
        if(min(table(sample_batch_id$BATCH_IDENTIFIER))==1)
        {
          single_sample_no = which(sample_batch_id$BATCH_IDENTIFIER%in%names(which(table(sample_batch_id$BATCH_IDENTIFIER)==1)))
          
          single_samples = sample_batch_id[single_sample_no,1]
          expression_data = expression_data[,-single_sample_no]
          sample_batch_id = sample_batch_id[-single_sample_no,]
          
          print(paste("The following samples are removed from the data as they belong to single sample batch"))
          print(single_samples)
        }
        time2 = proc.time()[3]
        batch = as.factor(as.character(sample_batch_id$BATCH_IDENTIFIER))
        modcombat = model.matrix(~1, data=sample_batch_id)
        combat_edata = ComBat(dat=expression_data, batch=batch, mod=modcombat, 
                              par.prior=TRUE, prior.plots=TRUE)
        write.table(combat_edata, file = paste(file_batch_effect,paste("Combat_",Title,".txt",sep=""),sep = "/"), sep="\t")
        print(paste("Batch effects removed using Combat in",round((proc.time()[3]- time2)/60,2),"mins"))
        if(send_email_indicator)
        {
          sendmail(from = sprintf("<Combat_batch_correction@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), subject = "Notification",
                   msg = paste("Batch effects removed using Combat in",round((proc.time()[3]- time2)/60,2),"mins"),
                   control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
        }
      }else{
        print("samples are not in order")
        stop()
      }
    }else{
      combat_edata = expression_data
    }
    
    ################
    #Batch correction check
    ################
    if(batch_correction_check)
    {
      time1 = proc.time()[3]
      gpca_result<-gPCA.batchdetect(x=t(combat_edata),batch=sample_batch_id$BATCH_IDENTIFIER,center=FALSE,
                                    scaleY=FALSE,filt=NULL,nperm=1000,seed=NULL)
      
      save(gpca_result, file = paste(file_batch_effect,paste("gpca_result",Title,".RDATA",sep=""),sep = "/"))
      print(paste("Total time spent in gPCA to check Batch correction is",round((proc.time()[3]-time1)/60,2),"mins"))
      print(paste("value of the test statistic delta is",gpca_result$delta))
      print(paste("p-value of the test is",gpca_result$p.val))
      
      if(send_email_indicator)
      {
        sendmail(from = sprintf("<gpca@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), 
                 subject = "Notification",
                 msg = paste(paste("Total time spent in gPCA to check Batch correction is"
                                   ,round((proc.time()[3]-time1)/60,2),"mins"),
                             paste("value of the test statistic delta is",gpca_result$delta),
                             paste("p-value of the test is",gpca_result$p.val),sep = " ||| "),
                 control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
      }
      rm(gpca_result)
      
    }
    ################
    #PCA To select no. of componenets in fastICA
    ################
    
    if(n_princomp_check)
    {
      file_pca = file.path(getwd(), "PCA_results")
      dir.create(file_pca, showWarnings = FALSE)
      
      time1 = proc.time()[3]
      prcomp = princomp(x=combat_edata, cor = prin_comp_cor)
      ev = (prcomp$sdev)^2
      ev_prop = ev/sum(ev)
      ev_prop_cum = cumsum(ev_prop)
      print(paste("first",min(which(ev_prop_cum>=0.95)),"components explain 95% of the total variance"))
      print(paste("first",min(which(ev_prop_cum>=0.96)),"components explain 96% of the total variance"))
      print(paste("first",min(which(ev_prop_cum>=0.97)),"components explain 97% of the total variance"))
      print(paste("first",min(which(ev_prop_cum>=0.98)),"components explain 98% of the total variance"))
      print(paste("first",min(which(ev_prop_cum>=0.99)),"components explain 99% of the total variance"))
      print(paste("first",min(which(ev_prop_cum>=0.995)),"components explain 99.5% of the total variance"))
      save(prcomp, file = paste(file_pca,paste("pca_result",Title,".RData",sep=""),sep = "/"))
      rm(prcomp)
      print(paste("Principal Components Analysis is done in",round((proc.time()[3]- time1)/60,2),"mins"))
      
      if(send_email_indicator)
      {
        sendmail(from = sprintf("<pca@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), 
                 subject = "Notification",
                 msg = paste(paste(paste("first",min(which(ev_prop_cum>=0.95)),"components explain 95% of the total variance"),
                                   paste("first",min(which(ev_prop_cum>=0.96)),"components explain 96% of the total variance"),
                                   paste("first",min(which(ev_prop_cum>=0.97)),"components explain 97% of the total variance"),
                                   paste("first",min(which(ev_prop_cum>=0.98)),"components explain 98% of the total variance"),
                                   paste("first",min(which(ev_prop_cum>=0.99)),"components explain 99% of the total variance"),
                                   paste("first",min(which(ev_prop_cum>=0.995)),"components explain 99.5% of the total variance"),
                                   paste("Principal Components Analysis is done in",round((proc.time()[3]- time1)/60,2),"mins"),
                                   sep = " ||| ")),
                 control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
      }
      if(choose_fastICA_ncomp)
      {
        fica_n_comp = readinteger("Choose number of components")
      }else{
        
        fica_n_comp = min(which(ev_prop_cum>=var_cutoff_to_choose_ncomp))
        print(paste("first",min(which(ev_prop_cum>=var_cutoff_to_choose_ncomp)),"components explain",var_cutoff_to_choose_ncomp, "of the total variance"))
      }
    }else{
      fica_n_comp = ncomp_without_pca
    }
    
    
    
    #######################
    #Implementation of fastICA algorithm
    #######################
    
    first_time = proc.time()[3]
    file_ica = file.path(getwd(), "ICA_all_iterations")
    dir.create(file_ica, showWarnings = FALSE)
    
    setwd(file_ica)
    
    fastICA_PARALLEL(X = combat_edata, n.comp = fica_n_comp, 
                     alg.typ = fastica_alg_typ, fun = fastica_fun, alpha = fastica_alpha,
                     method = fastica_method, row.norm = fastica_row_norm, maxit = fastica_maxit, 
                     tol = fastica_tol, verbose = fastica_verbose, w.init=NULL, seed = seed, Title = paste(Title,"ica", sep = "_"),
                     nsample_fastica = nsample_fastica, file_return_needed = FALSE, no_cores = no_cores_v1)
    
    total_elapsed_time = proc.time()[3] - first_time
    print(paste("total elapsed time for ica is", round(total_elapsed_time/60,2) , "mins"))
    if(send_email_indicator)
    {
      sendmail(from = sprintf("<fastica@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), 
               subject = "Notification",
               msg = paste("total elapsed time for ica is", round(total_elapsed_time/60,2) , "mins"),
               control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
    }
    
    if(CES_clustering)
    {
      
      ###################
      #Clustering of estimated sources by consensus sources estimation
      ###################
      time1 = proc.time()[3]
      
      # Combine the 〖ESM〗_(p×i)’s of all ICA runs together into a single matrix with p rows and i× number of ICA runs columns.
      
      CES_matrix = as.data.frame(matrix(0, dim(expression_data)[1], fica_n_comp*nsample_fastica))
      for(samp_no in 1:nsample_fastica)
      {
        load(paste(paste(Title,"ica", sep = "_"),samp_no,".RData", sep = "_"))
        
        CES_matrix[,((samp_no-1)*fica_n_comp+1):(samp_no*fica_n_comp)] = ica_result$S
        colnames(CES_matrix)[((samp_no-1)*fica_n_comp+1):(samp_no*fica_n_comp)] = paste("ICA",1:fica_n_comp,samp_no,sep = "_")
        rm(ica_result)
      }
      
      rownames(CES_matrix) = rownames(expression_data)
      CES_matrix = scale(CES_matrix)
      
      #calculating correlation coefficients between all ESM's from all runs.
      
      if(large_data)
      {
        correlation_matrix = fastCor(CES_matrix, nSplit = 25, upperTri = FALSE, verbose = TRUE)
        
      }else{
        correlation_matrix = cor(CES_matrix)
        
      }
      
      # Clustering of highly correlated ESs. In the present study, ESs were clustered together when the absolute values of the Pearson correlation coefficients between them were > 0.9
      
      similarity_matrix = ifelse(abs(correlation_matrix)>0.9,1,0)  
      
      diag(similarity_matrix) = 0
      
      print(paste("Total time spent in data preparation for CES clustering",round((proc.time()[3]-time1)/60,2),"mins"))
      if(send_email_indicator)
      {
        sendmail(from = sprintf("<CES_clustering@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), 
                 subject = "Notification",
                 msg = paste("Total time spent in data preparation for CES clustering",round((proc.time()[3]-time1)/60,2),"mins"),
                 control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
      }
      
      time1 = proc.time()[3]
      clusters = list()
      n_count = list()
      final_clusters = list()
      credibility_index = array(0,nsample_fastica*fica_n_comp)
      avg_source_ic = matrix(0, dim(CES_matrix)[1],fica_n_comp)
      if(CES_clustering_algo1)
      {
        for(iter in 1:(fica_n_comp))
        {
          time = proc.time()[3]
          clusters[[iter]] = list()
          n_count[[iter]] = array(0,dim(similarity_matrix)[2])
          for(samp_ica in 1:dim(similarity_matrix)[2])
          {
            x = union(samp_ica,which(similarity_matrix[,samp_ica]==1))
            
            if(is.null(x))
            {
              clusters[[iter]][[samp_ica]]= NA
              n_count[[iter]][samp_ica] = 0
            }else{
              clusters[[iter]][[samp_ica]]= x
              n_count[[iter]][samp_ica] = length(clusters[[iter]][[samp_ica]])
            }
            
          }
          
          
          final_clusters[[iter]] = clusters[[iter]][[which(n_count[[iter]]==max(n_count[[iter]]))[1]]]
          
          credibility_index[iter] = max(n_count[[iter]])/nsample_fastica
          
          avg_source_ic[,iter] = rowSums(CES_matrix[,final_clusters[[iter]]]%*%
                                           diag(sign(correlation_matrix[min(which(n_count[[iter]]==max(n_count[[iter]]))),
                                                                        final_clusters[[iter]]])))/(nsample_fastica*credibility_index[iter])
          
          print(paste("cluster found of length", max(n_count[[iter]]),"in",round((proc.time()[3]-time)/60,2),"mins"))
          for( i in final_clusters[[iter]])
          {
            similarity_matrix[i,] = 0
            similarity_matrix[,i] = 0
          }
          
          # The number of ESs in each cluster can be at most the total number of ICA runs (for the present study, that is 25). 
          # The higher the credibility index, the higher the chance of obtaining the ESs for which the negentropy converges to its global maxima. 
          # In the present study, the cut-off for the credibility index was fixed at 50%.
          # That is, clusters with a credibility index greater than 50% were only considered for obtaining the consensus estimated source matrix (CESM). 
          # The 〖CESM〗_(p×m) contains 〖CES〗_(p×1)’s from m clusters which had a credibility index greater than 50%. 
          # The characteristics of CESs are similar to those of ESs.
          
          if(credibility_index[iter]<0.5)
            break
          
          print(paste("iteration",iter,"completed in",round((proc.time()[3]-time)/60,2),"mins"))
          
        }
        final_avg_source_ic = avg_source_ic[,1:(max(which(credibility_index>0.5)))]
        rownames(final_avg_source_ic) = rownames(combat_edata)
        
        final_credibility_index = credibility_index[1:(max(which(credibility_index>0.5)))]
        file_consensus_ica = file.path(file_main, "ICA_Consensus_results")
        dir.create(file_consensus_ica, showWarnings = FALSE)
        setwd(file_consensus_ica)
        
        write.table(final_avg_source_ic, file = paste("Consensus_Estimated_Sources",Title, ".txt",sep = "_"), sep="\t")
        write.table(final_credibility_index, file = paste("Consensus_credibility_indices",Title, ".txt",sep = "_"), sep="\t")
        save(final_clusters,file = paste("Consensus_list_of_clusters",Title, ".RData",sep = "_"))
        
        ####################
        #updating Mixing matrix
        ####################
        
        dot_prod = t(final_avg_source_ic)%*%final_avg_source_ic
        V = solve(dot_prod)
        mix_matrix = V%*%t(final_avg_source_ic)%*%as.matrix(combat_edata)
        rownames(mix_matrix) = paste("V",1:dim(final_avg_source_ic)[2],sep = "")
        write.table(mix_matrix, file = paste("Consensus_Mixing_Matrix",Title, ".txt",sep = "_"), sep="\t")
        
      }else{
        for(iter in 1:(nsample_fastica*fica_n_comp))
        {
          #Another algorithm to cluster the ES's, but this was not used in the current paper
          
          time = proc.time()[3]
          clusters[[iter]] = list()
          count = 100000000
          n_count[[iter]] = array(0,dim(similarity_matrix)[2])
          
          for(samp_ica in 1:dim(similarity_matrix)[2])
          {
            if(!(samp_ica%in%count))
            {
              x = union(samp_ica,find_similar_components(data=similarity_matrix,
                                                         set_of_rows = which(similarity_matrix[,samp_ica]==1),set_of_cluster=NA))
              
              if(is.null(x))
              {
                clusters[[iter]][[samp_ica]]= NA
                n_count[[iter]][samp_ica] = 0
              }else{
                clusters[[iter]][[samp_ica]]= x
                n_count[[iter]][samp_ica] = length(clusters[[iter]][[samp_ica]])
              }
              count = union(count,x)
            }
            
          }
          
          
          final_clusters[[iter]] = clusters[[iter]][[which(n_count[[iter]]==max(n_count[[iter]]))[1]]]
          
          credibility_index[iter] = max(n_count[[iter]])/nsample_fastica
          
          avg_source_ic[,iter] = rowSums(CES_matrix[,final_clusters[[iter]]]%*%
                                           diag(sign(correlation_matrix[min(which(n_count[[iter]]==max(n_count[[iter]]))),
                                                                        final_clusters[[iter]]])))/(nsample_fastica*credibility_index[iter])
          
          print(paste("cluster found of length", max(n_count[[iter]]),"in",round((proc.time()[3]-time)/60,2),"mins"))
          for( i in final_clusters[[iter]])
          {
            similarity_matrix[i,] = 0
            similarity_matrix[,i] = 0
          }
          
          if(length(final_clusters[[iter]])==1)
            break
          
          print(paste("iteration",iter,"completed in",round((proc.time()[3]-time)/60,2),"mins"))
          
        }
        
        final_avg_source_ic = avg_source_ic[,1:(min(which(credibility_index==0))-2)]
        
        rownames(final_avg_source_ic) = rownames(combat_edata)
        final_credibility_index = credibility_index[1:(min(which(credibility_index==0))-2)]
        
        write.table(final_avg_source_ic, file = paste("Final Average Independent Components algo 2",Title, ".txt"), sep="\t")
        write.table(final_credibility_index, file = paste("Final Credibility Indices algo 2",Title, ".txt"), sep="\t")
        save(final_clusters,file = paste("Final list of clusters algo 2",Title, ".RData"))
        
        
        
        ####################
        #updating Mixing matrix
        ####################
        dot_prod = t(final_avg_source_ic)%*%final_avg_source_ic
        V = solve(dot_prod)
        mix_matrix = V%*%t(final_avg_source_ic)%*%as.matrix(combat_edata)
        rownames(mix_matrix) = paste("V",1:dim(final_avg_source_ic)[2],sep = "")
        write.table(mix_matrix, file = paste("Final Average Mix matrix algo 2",Title, ".txt"), sep="\t")
        
      }
      
      print(paste("Total time spent in CES clustering",round((proc.time()[3]-time1)/60,2),"mins"))
      if(send_email_indicator)
      {
        sendmail(from = sprintf("<CES_clustering@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), 
                 subject = "Notification",
                 msg = paste("Total time spent in CES clustering",round((proc.time()[3]-time1)/60,2),"mins"),
                 control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
      }
    }
    
    
    ####################
    #Consensus Clustering of the Data
    ####################
    
    if(consensus_clustering_of_data)
    {
      file_consensus_clustering = file.path(file_main, "Consensus_clustering_of_the_data")
      dir.create(file_consensus_clustering, showWarnings = FALSE)
      setwd(file_consensus_clustering)
      
      ###################
      #top mad gene exressions
      ###################
      
      time5 = proc.time()[3]
      mads = apply(combat_edata,1,mad)
      
      d = combat_edata[rev(order(mads))[1:no_mad_ge],]
      
      d = as.matrix(sweep(d,1, apply(d,1,median,na.rm=T)))
      print(paste("time taken to gene median center for gene expressions with 5000 mad the data is", round((proc.time()[3]-time5)/60,2),"mins"))
      if(send_email_indicator)
      {
        sendmail(from = sprintf("<consensus_clustering@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), 
                 subject = "Notification",
                 msg = paste("time taken to gene median center for gene expressions data with top 5000 mad is", round((proc.time()[3]-time5)/60,2),"mins"),
                 control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
      }
      
      time3 = proc.time()[3]
      results = ConsensusClusterPlus(d,maxK=consensus_maxK,reps=consensus_reps,pItem=consensus_pItem,pFeature=consensus_pFeature,
                                     title=paste(Title,"consensus_clustering", sep = "_"),clusterAlg=consensus_clusterAlg,distance=distance_function,seed=seed,plot="pdf",
                                     verbose = TRUE, writeTable = TRUE)
      save(results, file = paste(Title,"consensus_clustering_results.RData", sep = "_"))
      rm(results)
      print(paste("time taken to do consensus clustering on the data is", round((proc.time()[3]-time3)/60,2),"mins for distance function"))
      if(send_email_indicator)
      {
        sendmail(from = sprintf("<consensus_clustering@%s>", Sys.info()[4]), to = paste("<",email_id,">",sep = ""), 
                 subject = "Notification",
                 msg = paste("time taken to do consensus clustering on the data is", round((proc.time()[3]-time3)/60,2),"mins for distance function"),
                 control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
      }
      
    }
    ############################
    #Detection of extreme-valued genomic region (DEGR)
    ############################
    
    if(permutation_testing_for_degr_indicator){
      
      
      # Collapsing probeset-level weights of the CESM.
      
      if(gene_level_mapping_indicator){
        
        glmi_time = proc.time()[3]
        
        matching_vec = match(data_with_chr_bp_mapping$PROBESET,rownames(final_avg_source_ic))
        
        data_with_chr_bp_mapping = data_with_chr_bp_mapping[!(is.na(matching_vec)),]
        
        matching_vec = matching_vec[!(is.na(matching_vec))]
        
        final_avg_source_ic = cbind(as.data.frame(final_avg_source_ic[matching_vec,]),data_with_chr_bp_mapping[,c('ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping')])
        
        final_avg_source_ic_v1 = final_avg_source_ic
        original_mapping_file = data_with_chr_bp_mapping
        
        if(collapse_function == "max")
        {
          define_function = function(x) x[which(abs(x)==max(abs(x)))]
        }else if(collapse_function == "mean"){
          define_function = function(x) mean(x)
        }else{
          define_function = function_for_collapse
        }
        
        unique_genes = as.character(unique(final_avg_source_ic_v1$ENTREZID))
        gene_position_in_data = NULL
        gene_position_in_data[1] = 0
        final_comps = which(credibility_index>=0.5)
        collapsed_data = matrix(0,length(unique_genes), length(final_comps))
        for(uniq_gene_num in 1:length(unique_genes)){
          gene_position_in_data[uniq_gene_num+1] = max(which(final_avg_source_ic_v1$ENTREZID==unique_genes[uniq_gene_num]))
          
          collapsed_data[uniq_gene_num,] = apply(final_avg_source_ic_v1[(gene_position_in_data[uniq_gene_num]+1):(gene_position_in_data[uniq_gene_num+1]),final_comps],2,function(x) define_function(x))
        }
        
        row.names(collapsed_data) = unique_genes
        
        sub_data_gene_mapping = sqldf("select distinct ENTREZID,GENETITLE, SYMBOL, CHR_Mapping, BP_Mapping from data_with_chr_bp_mapping order by CHR_Mapping, BP_Mapping")
        data_with_chr_bp_mapping = sub_data_gene_mapping
        final_avg_source_ic_v1 = collapsed_data
        
        setwd(file_consensus_ica)
        
        write.table(collapsed_data, file = paste("Genelevel_Consensus_Estimated_Sources",Title, ".txt",sep = "_"), sep="\t")
        
        print(paste("gene level mapping is done in", round((proc.time()[3]-glmi_time)/60,3),"mins"))
        
      }else{
        final_comps = which(credibility_index>=0.5)
        
        matching_vec = match(data_with_chr_bp_mapping$ENTREZID,rownames(final_avg_source_ic))
        
        data_with_chr_bp_mapping = data_with_chr_bp_mapping[!(is.na(matching_vec)),]
        
        matching_vec = matching_vec[!(is.na(matching_vec))]
        
        final_avg_source_ic = cbind(as.data.frame(final_avg_source_ic[matching_vec,]),data_with_chr_bp_mapping[,c('ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping')])
        
        final_avg_source_ic_v1 = final_avg_source_ic
        original_mapping_file = data_with_chr_bp_mapping
        
      }
      
      for(fdr in FDR)
      {
        for(cl in CL)
        {
          for(state_deciding_cutoffs in state_deciding_cutoff)
          {
            for(probe_no_for_gaus_kernels in probe_no_for_gaus_kernel)
            {
              for(collapse_functions in collapse_function)
              {
                permutation_testing_to_find_degr(dataset=as.matrix(final_avg_source_ic_v1[,final_comps]), 
                                                 data_with_chr_bp_mapping = data_with_chr_bp_mapping,
                                                 intervals_for_gaus_kernel = intervals_for_gaus_kernel,
                                                 probe_no_for_gaus_kernel = probe_no_for_gaus_kernels, 
                                                 Title = Title,
                                                 FDR = fdr,
                                                 CL = cl,
                                                 state_deciding_cutoff = state_deciding_cutoffs,
                                                 min_probesets=min_probesets,
                                                 set_seed = set_seed,
                                                 full_dataset = final_avg_source_ic,
                                                 original_mapping_file = original_mapping_file,
                                                 mix_matrix = mix_matrix,
                                                 collapse_function = collapse_functions,
                                                 gene_level_mapping_indicator = gene_level_mapping_indicator,TACNA_profiling = TACNA_profiling, file_main = file_main)
                
              }
            }
          }
        }
      }
    }
  }
}


gene_ranking_in_consensus_sources = function(ica_extreme_valued_dataset, title)
{
  colnames(ica_extreme_valued_dataset) = unlist(lapply(ica_extreme_valued_dataset[1,], as.character))
  
  
  ica_extreme_valued_dataset = ica_extreme_valued_dataset[-1,]
  
  ica_extreme_valued_dataset = ica_extreme_valued_dataset[,c('SYMBOL'
                                               , 'ENTREZID'
                                               , 'GENETITLE'
                                               , 'CHR_Mapping'
                                               , 'BP_Mapping'
                                               , 'BP_Mapping_for_plot'
                                               , 'CES'
                                               , 'Value'
                                               , 'state_of_the_weight' )]
  
  for(cols in c( 'CHR_Mapping'
                 , 'BP_Mapping'
                 , 'BP_Mapping_for_plot'
                 , 'CES'
                 , 'Value'
                 , 'state_of_the_weight' ))
  {
    ica_extreme_valued_dataset[,cols]= as.numeric(as.character(ica_extreme_valued_dataset[,cols]))
  }
  
  ica_extreme_valued_dataset$abs_value = abs(ica_extreme_valued_dataset$Value)
  
  ica_extreme_valued_dataset_v1 = sqldf("select a.* from ica_extreme_valued_dataset a order by CES,abs_value desc")
  
  unique_comps = table(ica_extreme_valued_dataset_v1$CES)
  ica_extreme_valued_dataset_v1$rank = NA
  ica_extreme_valued_dataset_v1$rank_percentile = NA
  ica_extreme_valued_dataset_v1$dev_from_local_max = NA
  
  for( i in names(unique_comps))
  {
    
    ica_extreme_valued_dataset_v1$rank[which(ica_extreme_valued_dataset_v1$CES==i)] = rank(-ica_extreme_valued_dataset_v1$abs_value[which(ica_extreme_valued_dataset_v1$CES==i)], ties.method = "average")
    ica_extreme_valued_dataset_v1$rank_percentile[which(ica_extreme_valued_dataset_v1$CES==i)] = ica_extreme_valued_dataset_v1$rank[which(ica_extreme_valued_dataset_v1$CES==i)]/unique_comps[which(names(unique_comps)==i)]
    ica_extreme_valued_dataset_v1$dev_from_local_max[which(ica_extreme_valued_dataset_v1$CES==i)] = 1 - ica_extreme_valued_dataset_v1$abs_value[which(ica_extreme_valued_dataset_v1$CES==i)]/max(ica_extreme_valued_dataset_v1$abs_value[which(ica_extreme_valued_dataset_v1$CES==i)])
  }
  
  ica_extreme_valued_dataset_v1$dev_from_global_max =  1 - ica_extreme_valued_dataset_v1$abs_value/max(ica_extreme_valued_dataset_v1$abs_value)
  
  ica_extreme_valued_dataset_v2 = sqldf("select SYMBOL, ENTREZID, GENETITLE, CHR_Mapping, BP_Mapping, BP_Mapping_for_plot, min(rank_percentile) as min_rank_percentile 
                                 , count(*) as no_of_CESs from ica_extreme_valued_dataset_v1 group by 1,2,3,4,5,6")
  
  
  ica_extreme_valued_dataset_v3 = sqldf("select a.*, b.dev_from_local_max, b.dev_from_global_max from ica_extreme_valued_dataset_v2 a 
                                 left join ica_extreme_valued_dataset_v1 b
                                 on a.ENTREZID = b.ENTREZID
                                 and a.min_rank_percentile = b.rank_percentile")
  
  print(ica_extreme_valued_dataset_v3[which(ica_extreme_valued_dataset_v3$ENTREZID%in%names(table(ica_extreme_valued_dataset_v3$ENTREZID))[which(table(ica_extreme_valued_dataset_v3$ENTREZID)>1)]),])
  
  ica_extreme_valued_dataset_v4 = sqldf("select SYMBOL, ENTREZID, GENETITLE, CHR_Mapping, BP_Mapping, BP_Mapping_for_plot, min_rank_percentile 
                                 , no_of_CESs, min(dev_from_local_max) as dev_from_local_max, min(dev_from_global_max) as dev_from_global_max
                                 from ica_extreme_valued_dataset_v3 group by 1,2,3,4,5,6,7,8")
  
  write.table(ica_extreme_valued_dataset_v1, file = paste(title, "extreme_valued_all_gene_info.txt", sep = "_"), sep = "\t", row.names = FALSE)
  write.table(ica_extreme_valued_dataset_v3, file = paste(title, "extreme_valued_min_rank_gene_info.txt", sep = "_"), sep = "\t", row.names = FALSE)
  write.table(ica_extreme_valued_dataset_v4, file = paste(title, "extreme_valued_min_rank_min_dev_from_max_gene_info.txt", sep = "_"), sep = "\t", row.names = FALSE)
  
  
  
}


fastICA_PARALLEL <- function (X, n.comp, alg.typ = c("parallel","deflation"),
                     fun = c("logcosh", "exp"),
                     alpha = 1, method = c("R", "C"),
                     row.norm = FALSE, maxit = 200, tol = 1e-04,
                     verbose = FALSE, w.init=NULL, seed = 12345, Title = "",nsample_fastica = 25, file_return_needed = TRUE, no_cores = detectCores())
{
  #copied from FastICA source code in R and updated the parallelization part
  
  dd <- dim(X)
  d <- dd[dd != 1L]
  if (length(d) != 2L)
    stop("data must be matrix-conformal")
  X <- if (length(d) != length(dd)) matrix(X, d[1L], d[2L])
  else as.matrix(X)
  
  if (alpha < 1 || alpha > 2)
    stop("alpha must be in range [1,2]")
  method <- match.arg(method)
  alg.typ <- match.arg(alg.typ)
  fun <- match.arg(fun)
  n <- nrow(X)
  p <- ncol(X)
  
  if (n.comp > min(n, p)) {
    message("'n.comp' is too large: reset to ", min(n, p))
    n.comp <- min(n, p)
  }
  
  if (verbose) message("Centering")
  
  X <- scale(X, scale = FALSE)
  
  X <- if (row.norm) t(scale(X, scale=row.norm)) else t(X)
  
  if (verbose) message("Whitening")
  V <- X %*% t(X)/n
  
  s <- La.svd(V)
  D <- diag(c(1/sqrt(s$d)))
  
  K <- D %*% t(s$u)
  K <- matrix(K[1:n.comp, ], n.comp, p)
  X1 <- K %*% X
  
  if (method == "R") {
    ica_r = function(iter,X1_ =  X1,alg.typ_ = alg.typ, n.comp_ = n.comp, tol_ = tol, fun_ = fun,
                     alpha_ = alpha, maxit_ = maxit, verbose_ = verbose, w.init_ = w.init, 
                     seed_ = seed, Title_ = Title,
                     file_return_needed_ = file_return_needed, nsample_fastica_=nsample_fastica)
    {
      set.seed(seed=seed_+iter)
      current_seed = seed_+iter
      if(is.null(w.init_)){
        w.init_ <- matrix(rnorm(n.comp_^2),n.comp_,n.comp_)
        
      }else {
        if(!is.matrix(w.init_) || length(w.init_) != (n.comp_^2))
          stop("w.init is not a matrix or is the wrong size")
      }
      
      a <- if (alg.typ_ == "deflation"){
        ica.R.def(X = X1_, n.comp = n_comp_, tol = tol_, fun = fun_,
                  alpha = alpha_, maxit = maxit_, verbose = verbose_, w.init = w.init_)
      }else if (alg.typ_ == "parallel"){
        ica.R.par(X = X1_, n.comp = n_comp_, tol = tol_, fun = fun_,
                  alpha = alpha_, maxit = maxit_, verbose = verbose_, w.init = w.init_)
      }
      
      w <- a[['W']] %*% K
      S <- w %*% X
      A <- t(w) %*% solve(w %*% t(w))
      
      ica_result  = list(X = t(X), K = t(K), W = t(a[['W']]), A = t(A), S = t(S), comments = a[['comments']], seed = current_seed)
      
      if(length(a[['comments']])<(maxit_ - 1))
      {
        save(ica_result, file = paste(Title_,iter,".RData", sep = "_"))
      }else{
        rm(ica_result)
        iter_ = iter
        X1 =  X1_
        alg.typ = alg.typ_
        n.comp = n.comp_
        tol = tol_
        fun = fun_
        alpha = alpha_
        maxit = maxit_
        verbose = verbose_
        w.init = NULL
        Title = Title_
        file_return_needed = file_return_needed_
        nsample_fastica = nsample_fastica_
        
        ica_r(iter = iter_,X1_ =  X1,alg.typ_ = alg.typ, n.comp_ = n.comp, tol_ = tol, fun_ = fun,
              alpha_ = alpha, maxit_ = maxit, verbose_ = verbose, w.init_ = w.init, 
              seed_ = current_seed+nsample_fastica-iter, Title_ = Title,
              file_return_needed_ = file_return_needed, nsample_fastica_ = nsample_fastica)
      }
      
      if(file_return_needed_){
        return(ica_result)
      }else{
        return(1)
      }
      
    }
    # no_cores <- detectCores() 
    cl <- makeCluster(no_cores, type = "FORK")
    
    parLapply(cl, 1:nsample_fastica,ica_r)
    stopCluster(cl)
    
    
  } else if (method == "C") {
    a <- .C("icainc_JM",
            as.double(X),
            as.double(w.init),
            as.integer(p),
            as.integer(n),
            as.integer(n.comp),
            as.double(alpha),
            as.integer(1),
            as.integer(row.norm),
            as.integer(1L + (fun == "exp")),
            as.integer(maxit),
            as.double(tol),
            as.integer(alg.typ != "parallel"),
            as.integer(verbose),
            X = double(p * n),
            K = double(n.comp * p),
            W = double(n.comp * n.comp),
            A = double(p * n.comp),
            S = double(n.comp * n))
    X1 <- matrix(a$X, n, p)
    K <- matrix(a$K, p, n.comp)
    W <- matrix(a$W, n.comp, n.comp)
    A <- matrix(a$A, n.comp, p)
    S <- matrix(a$S, n, n.comp)
    list(X = X1, K = K, W = W, A = A, S = S)
  }
}

ica.R.def <-
  function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init)
  {
    if (verbose && fun == "logcosh")
      message("Deflation FastICA using logcosh approx. to neg-entropy function")
    if (verbose && fun =="exp")
      message("Deflation FastICA using exponential approx. to neg-entropy function")
    n <- nrow(X)
    p <- ncol(X)
    W <- matrix(0, n.comp, n.comp)
    for (i in 1:n.comp) {
      if (verbose) message("Component ", i)
      w <- matrix(w.init[i,], n.comp, 1)
      if (i > 1) {
        t <- w
        t[1:length(t)] <- 0
        for (u in 1:(i - 1)) {
          k <- sum(w * W[u, ])
          t <- t + k * W[u, ]
        }
        w <- w - t
      }
      w <- w/sqrt(sum(w^2))
      lim <- rep(1000, maxit)
      it <- 1
      if (fun == "logcosh") {
        comments = NULL
        while (lim[it] > tol && it < maxit) {
          wx <- t(w) %*% X
          gwx <- tanh(alpha * wx)
          gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
          xgwx <- X * gwx
          v1 <- apply(xgwx, 1, FUN = mean)
          g.wx <- alpha * (1 - (tanh(alpha * wx))^2)
          v2 <- mean(g.wx) * w
          w1 <- v1 - v2
          w1 <- matrix(w1, n.comp, 1)
          it <- it + 1
          if (i > 1) {
            t <- w1
            t[1:length(t)] <- 0
            for (u in 1:(i - 1)) {
              k <- sum(w1 * W[u, ])
              t <- t + k * W[u, ]
            }
            w1 <- w1 - t
          }
          w1 <- w1/sqrt(sum(w1^2))
          lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
          comments[it-1] = paste("Iteration ", it - 1, " tol = ", format(lim[it]))
          if (verbose)
            message("Iteration ", it - 1, " tol = ", format(lim[it]))
          w <- matrix(w1, n.comp, 1)
        }
      }
      if (fun == "exp") {
        comments = NULL
        while (lim[it] > tol && it < maxit) {
          wx <- t(w) %*% X
          gwx <- wx * exp(-(wx^2)/2)
          gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
          xgwx <- X * gwx
          v1 <- apply(xgwx, 1, FUN = mean)
          g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
          v2 <- mean(g.wx) * w
          w1 <- v1 - v2
          w1 <- matrix(w1, n.comp, 1)
          it <- it + 1
          if (i > 1) {
            t <- w1
            t[1:length(t)] <- 0
            for (u in 1:(i - 1)) {
              k <- sum(w1 * W[u, ])
              t <- t + k * W[u, ]
            }
            w1 <- w1 - t
          }
          w1 <- w1/sqrt(sum(w1^2))
          lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
          comments[it-1] = paste("Iteration ", it - 1, " tol = ", format(lim[it]))
          if (verbose)
            message("Iteration ", it - 1, " tol = ", format(lim[it]))
          w <- matrix(w1, n.comp, 1)
        }
      }
      W[i, ] <- w
    }
    return(list(W = W, comments = comments))
  }

ica.R.par <- function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init)
{
  Diag <- function(d) if(length(d) > 1L) diag(d) else as.matrix(d)
  n <- nrow(X)
  p <- ncol(X)
  W <- w.init
  sW <- La.svd(W)
  W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
  W1 <- W
  lim <- rep(1000, maxit)
  it <- 1
  comments = NULL
  if (fun == "logcosh") {
    if (verbose)
      message("Symmetric FastICA using logcosh approx. to neg-entropy function")
    
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- tanh(alpha * wx)
      v1 <- gwx %*% t(X)/p
      g.wx <- alpha * (1 - (gwx)^2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      W <- W1
      comments[it] = paste("Iteration ", it, " tol = ", format(lim[it]))
      if (verbose)
        message("Iteration ", it, " tol = ", format(lim[it + 1]))
      it <- it + 1
    }
  }
  if (fun == "exp") {
    if (verbose)
      message("Symmetric FastICA using exponential approx. to neg-entropy function")
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- wx * exp(-(wx^2)/2)
      v1 <- gwx %*% t(X)/p
      g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      W <- W1
      
      comments[it] = paste("Iteration ", it, " tol = ", format(lim[it]))
      if (verbose)
        message("Iteration ", it, " tol = ", format(lim[it + 1]))
      it <- it + 1
    }
  }
  return(list(W = W, comments = comments))
}


readinteger <- function(text_to_show)
{ 
  n <- readline(prompt=paste(text_to_show,": "))
  if(!grepl("^[0-9]+$",n))
  {
    return(readinteger())
  }
  
  return(as.integer(n))
}

#The below function identifies netwrok of rows and columns of a dataset which have connection between them using the value "1". But this function was not used in this paper.
find_similar_components = function(data, set_of_rows, set_of_cluster)
{set_of_rows = set_of_rows[which(!(set_of_rows%in%set_of_cluster))]
if(length(set_of_rows)==0)
{
  cluster_points = NULL
  return(cluster_points)
}else{
  cluster_set = NULL
  for(i in set_of_rows)
  {
    cluster_set = rbind(cluster_set, cbind(i,which(data[i,]==1)))
  }
  
  cluster_set_2 = NULL
  set_of_columns = unique(cluster_set[which(!(cluster_set[,2]%in%set_of_rows)),2])
  for(i in set_of_columns)
  {
    cluster_set_2 = rbind(cluster_set_2, cbind(which(data[,i]==1),i))
  }
  
  new_set_of_rows = unique(cluster_set_2[which(!(cluster_set_2[,1]%in%union(set_of_columns,set_of_rows))),1])
  
  set = union(unique(cluster_set_2[,1]),set_of_columns)
  if(length(set_of_cluster)==1)
  {
    set_of_cluster_v1 = set[which(!(set%in%new_set_of_rows))]
  }else{
    set_of_cluster_v1 = union(set[which(!(set%in%new_set_of_rows))],set_of_cluster)
    
  }
  cluster_points = union(new_set_of_rows,union(set_of_cluster_v1,find_similar_components(data,new_set_of_rows,set_of_cluster_v1)))
  return(cluster_points)
}

}

# Smoothing the CESM. 
# To minimize the effect of outliers, smoothing is applied on the gene-level weights of CESs. 
# Smoothing of the weights is performed per chromosome. 
# If there are g number of genes in a single chromosome, then the smoothing coefficients C_k1, C_k2, …, C_kg are generated for the kth gene from the truncated normal distribution (TND). 

sliding_window_smoothing = function(dataset, data_with_chr_bp_mapping, CN_dataset_genelevel,CN_dataset,
                                    intervals_for_gaus_kernel = seq(100000,2000000,by=10000),collapse_function = "max",function_for_collapse,
                                    probe_no_for_gaus_kernel = 10, Title = "",
                                    set_seed = 123456, min_probesets = 5, file_main = file.path(getwd())){
  
  full_time = proc.time()[3]
  
  file_SCNA = file.path(file_main, paste("TACNAP_CN_comparison",sep = "_"))
  dir.create(file_SCNA, showWarnings = FALSE)
  setwd(file_SCNA)
  
  #arrange expression dataset rows in accordance with genomic maping file
  
  matching_vec = match(data_with_chr_bp_mapping$PROBESET,rownames(dataset))
  data_with_chr_bp_mapping = data_with_chr_bp_mapping[!(is.na(matching_vec)),]
  matching_vec = matching_vec[!(is.na(matching_vec))]
  
  #arrange expression dataset columns in accordance with copy number data
  
  match_vec_2 = match(colnames(CN_dataset),colnames(dataset))
  dataset_v1 = dataset[matching_vec,match_vec_2[!is.na(match_vec_2)]]
  CN_dataset_v1 = CN_dataset[,which(!is.na(match_vec_2))]
  
  na_values = which(is.na(CN_dataset_v1[,1]))
  
  dataset_v1 = dataset_v1[-na_values,]
  CN_dataset_v1 = CN_dataset_v1[-na_values,]
  cor_vec_probe = array(0,dim(CN_dataset_v1)[1])
  cor_mat_probe = as.data.frame(cbind(unlist(lapply(rownames(dataset_v1),as.character)),cor_vec_probe))
  cor_mat_probe$cor_vec_probe = as.numeric(as.character(cor_mat_probe$cor_vec_probe))
  colnames(cor_mat_probe)[1] = "Probes"
  for(i in 1:dim(cor_mat_probe)[1])
  {
    cor_mat_probe$cor_vec_probe[i] = cor(t(dataset_v1[i,]), t(CN_dataset_v1[i,]), use = "complete.obs")
  }
  write.table(cor_mat_probe, file = paste(Title,"Probe_wise_correlation.txt", sep = "_"), sep = "\t", row.names = FALSE)
  write.table(dataset_v1, file = paste(Title,"Non_missing_TACNAP_probelevel_expression_data.txt", sep = "_"), sep = "\t")
  write.table(CN_dataset_v1, file = paste(Title,"Non_missing_probelevel_CN_data.txt", sep = "_"), sep = "\t", row.names = FALSE)
  
  dataset = cbind(as.data.frame(dataset[matching_vec,]),data_with_chr_bp_mapping[,c('ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping')])
  
  
  # Collapsing probeset-level weights of the CESM. 
  # The CESs identified in the GEO-dataset, CCLE-dataset and GDSC-dataset contain probeset-level weights. 
  # As multiple probesets can target a single gene, many genomic regions have a high chance of co-localizing extreme values corresponding to these probeset-level weights. 
  # Therefore, probeset-level weights need to be collapsed to gene-level weights. 
  # In the DEGR algorithm, out of multiple probeset-level weights corresponding to the same gene, the probeset-level weight with the highest absolute value is retained as its gene-level weight.
  
    if(collapse_function == "max")
    {
      define_function = function(x) x[which(abs(x)==max(abs(x)))[1]]
    }else if(collapse_function == "mean"){
      define_function = function(x) mean(x)
    }else{
      define_function = function_for_collapse
    }
    
    unique_genes = as.character(unique(dataset$ENTREZID))
    gene_position_in_data = NULL
    gene_position_in_data[1] = 0
    final_comps = 1:(dim(dataset)[2]-5)
    collapsed_data = matrix(0,length(unique_genes), length(final_comps))
    for(uniq_gene_num in 1:length(unique_genes)){
      gene_position_in_data[uniq_gene_num+1] = max(which(dataset$ENTREZID==unique_genes[uniq_gene_num]))
      
      collapsed_data[uniq_gene_num,] = apply(dataset[(gene_position_in_data[uniq_gene_num]+1):(gene_position_in_data[uniq_gene_num+1]),final_comps],2,function(x) define_function(x))
    }
    
    rownames(collapsed_data) = unique_genes
    colnames(collapsed_data) = colnames(dataset)[final_comps]
    sub_data_gene_mapping = sqldf("select distinct ENTREZID,GENETITLE, SYMBOL, CHR_Mapping, BP_Mapping from data_with_chr_bp_mapping order by CHR_Mapping, BP_Mapping")
    data_with_chr_bp_mapping = sub_data_gene_mapping
    dataset = collapsed_data
    
  
  
  row_num = dim(dataset)[1]
  if(row_num!=dim(data_with_chr_bp_mapping)[1])
  {
    print("two datasets are not of equal rows")
  }else{
    if(is.null(colnames(dataset)))
    {
      colnames(dataset) = paste("V",1:dim(dataset)[2],sep ="")
    }
    ##########################
    #Adjustment for Plot
    ###########################
    
    
    chromosome_seq = c(0,cumsum(table(as.numeric(as.character(data_with_chr_bp_mapping$CHR_Mapping)))))
    
    data_with_chr_bp_mapping$BP_Mapping_v1 = data_with_chr_bp_mapping$BP_Mapping/10
    
    for(i in 2:max(data_with_chr_bp_mapping$CHR_Mapping))
    {
      data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]=data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]+data_with_chr_bp_mapping$BP_Mapping_v1[chromosome_seq[i]]
    }
    
    label_pos = sqldf("select distinct CHR_Mapping, avg(BP_Mapping_v1) as BP_Mapping from data_with_chr_bp_mapping group by 1 order by 1")
    
    
    
    ########################
    #Finding proper interval for sliding Gaussian Kernel
    #Choosing that interval where 10 or more no. of probesets 
    #corresponding to that chromosome are there in +/- 3*interval for 95% of the cases
    ########################
    quantile_5 = list()
    for( k in 1:max(data_with_chr_bp_mapping$CHR_Mapping))
    {
      quantile_5[[k]] = array(0,length(intervals_for_gaus_kernel))
      for(i in 1:length(intervals_for_gaus_kernel))
      {
        frequency_of_neigh_probes = NULL
        rows = list()
        for(j in (chromosome_seq[k]+1):(chromosome_seq[k+1]))
        {
          frequency_of_neigh_probes[j-chromosome_seq[k]] = length(which(abs(data_with_chr_bp_mapping$BP_Mapping[j] - data_with_chr_bp_mapping$BP_Mapping[c((chromosome_seq[k]+1):(chromosome_seq[k+1]))]) < 3*intervals_for_gaus_kernel[i]))
          
        }
        quantile_5[[k]][i] = quantile(frequency_of_neigh_probes,0.05)%/%1
        
        if(quantile_5[[k]][i]>=probe_no_for_gaus_kernel){
          print(paste("found",probe_no_for_gaus_kernel, "or more number of probe sets in 95% of the chromosome at interval", 
                      intervals_for_gaus_kernel[i],"for Chromosome",k))
          quantile_5[[k]][which(quantile_5[[k]]==0)] = NA
          break
        } 
      }
      
    }
    
    
    ###########################
    #Creating the Density Matrix to get sliding gaussian Kernel
    ###########################
    
    density_matrix = matrix(0, row_num, row_num)
    
    for(j in 1:row_num)
    {
      rows[[j]] = which(data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[j])
      k = data_with_chr_bp_mapping$CHR_Mapping[j]
      s_x = as.numeric(data_with_chr_bp_mapping$BP_Mapping[rows[[j]]])
      s_mean = as.numeric(data_with_chr_bp_mapping$BP_Mapping[j]) 
      
      den_tnorm = dtruncnorm(s_x,
                             a = s_x[1],
                             b = s_x[length(s_x)],
                             mean = s_mean, 
                             sd = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[k]])))])
      
      density_matrix[j,rows[[j]]] = den_tnorm/sum(den_tnorm)
    }
    
    
    smoothed_dataset = density_matrix%*%as.matrix(dataset[,final_comps])
    rownames(smoothed_dataset) = rownames(dataset)
    write.table(smoothed_dataset, file = paste(Title,"smoothed_expression_data.txt", sep = "_"), sep = "\t", row.names = FALSE)
    
    match_vec_2 = match(colnames(CN_dataset_genelevel),colnames(smoothed_dataset))
    smoothed_dataset_v1 = smoothed_dataset[,match_vec_2[!is.na(match_vec_2)]]
    CN_dataset_genelevel_v1 = CN_dataset_genelevel[,which(!is.na(match_vec_2))]
    
    cor_vec_sample = array(0,dim(CN_dataset_genelevel_v1)[2])
    cor_mat_sample = as.data.frame(cbind(unlist(lapply(colnames(CN_dataset_genelevel_v1),as.character)),cor_vec_sample))
    cor_mat_sample$cor_vec_sample = as.numeric(as.character(cor_mat_sample$cor_vec_sample))
    colnames(cor_mat_sample)[1] = "Sample"
    for(i in 1:dim(cor_mat_sample)[1])
    {
      cor_mat_sample[i,2] = cor(smoothed_dataset_v1[,i], CN_dataset_genelevel_v1[,i], use = "complete.obs")
    }
    write.table(cor_mat_sample, file = paste(Title,"Sample_wise_correlation_pearson.txt", sep = "_"), sep = "\t", row.names = FALSE)
    
    cor_vec_sample = array(0,dim(CN_dataset_genelevel_v1)[2])
    cor_mat_sample = as.data.frame(cbind(unlist(lapply(colnames(CN_dataset_genelevel_v1),as.character)),cor_vec_sample))
    cor_mat_sample$cor_vec_sample = as.numeric(as.character(cor_mat_sample$cor_vec_sample))
    colnames(cor_mat_sample)[1] = "Sample"
    for(i in 1:dim(cor_mat_sample)[1])
    {
      cor_mat_sample[i,2] = cor(smoothed_dataset_v1[,i], CN_dataset_genelevel_v1[,i], use = "complete.obs", method = "kendall")
    }
    write.table(cor_mat_sample, file = paste(Title,"Sample_wise_correlation_kendall.txt", sep = "_"), sep = "\t", row.names = FALSE)
    
    cor_vec_sample = array(0,dim(CN_dataset_genelevel_v1)[2])
    cor_mat_sample = as.data.frame(cbind(unlist(lapply(colnames(CN_dataset_genelevel_v1),as.character)),cor_vec_sample))
    cor_mat_sample$cor_vec_sample = as.numeric(as.character(cor_mat_sample$cor_vec_sample))
    colnames(cor_mat_sample)[1] = "Sample"
    for(i in 1:dim(cor_mat_sample)[1])
    {
      cor_mat_sample[i,2] = cor(smoothed_dataset_v1[,i], CN_dataset_genelevel_v1[,i], use = "complete.obs", method = "spearman")
    }
    write.table(cor_mat_sample, file = paste(Title,"Sample_wise_correlation_spearman.txt", sep = "_"), sep = "\t", row.names = FALSE)
    
    cor_vec_GENE_pearson = array(0,dim(CN_dataset_genelevel_v1)[1])
    cor_vec_GENE_kendall = array(0,dim(CN_dataset_genelevel_v1)[1])
    cor_vec_GENE_spearman = array(0,dim(CN_dataset_genelevel_v1)[1])
    var_vec_GENE = array(0,dim(CN_dataset_genelevel_v1)[1])
    
    cor_mat_GENE = as.data.frame(cbind(unlist(lapply(rownames(smoothed_dataset_v1),as.character)),cor_vec_GENE_pearson,cor_vec_GENE_kendall,cor_vec_GENE_spearman, var_vec_GENE))
    cor_mat_GENE$cor_vec_GENE_pearson = as.numeric(as.character(cor_mat_GENE$cor_vec_GENE_pearson))
    cor_mat_GENE$cor_vec_GENE_kendall = as.numeric(as.character(cor_mat_GENE$cor_vec_GENE_kendall))
    cor_mat_GENE$cor_vec_GENE_spearman = as.numeric(as.character(cor_mat_GENE$cor_vec_GENE_spearman))
    cor_mat_GENE$var_vec_GENE = as.numeric(as.character(cor_mat_GENE$var_vec_GENE))
    colnames(cor_mat_GENE)[1] = "GENEs"
    CN_dataset_genelevel_v1 = as.matrix(CN_dataset_genelevel_v1)
    for(i in 1:dim(cor_mat_GENE)[1])
    {
      cor_mat_GENE$cor_vec_GENE_pearson[i] = cor(smoothed_dataset_v1[i,], CN_dataset_genelevel_v1[i,], use = "na.or.complete")
      cor_mat_GENE$cor_vec_GENE_kendall[i] = cor(smoothed_dataset_v1[i,], CN_dataset_genelevel_v1[i,], use = "na.or.complete", method = "kendall")
      cor_mat_GENE$cor_vec_GENE_spearman[i] = cor(smoothed_dataset_v1[i,], CN_dataset_genelevel_v1[i,], use = "na.or.complete", method = "spearman")
      cor_mat_GENE$var_vec_GENE[i] = var(CN_dataset_genelevel_v1[i,], na.rm = TRUE)
    }
    write.table(cor_mat_GENE, file = paste(Title,"Gene_wise_correlation_vector.txt", sep = "_"), sep = "\t", row.names = FALSE)
    rm(cor_mat_GENE)
    
    probe_gene_cormat= probe_gene_wise_cn_TACNAP_cor_mat(exp_data = smoothed_dataset_v1,cn_data=CN_dataset_genelevel_v1, cor_method = "pearson", no_cores = no_cores_v1)
    
    
    write.table(probe_gene_cormat, file = paste(Title,"Gene_wise_correlation_matrix_pearson.txt", sep = "_"), sep = "\t")
    rm(probe_gene_cormat)
    probe_gene_cormat= probe_gene_wise_cn_TACNAP_cor_mat(exp_data = smoothed_dataset_v1,cn_data=CN_dataset_genelevel_v1, cor_method = "kendall", no_cores = no_cores_v1)
    
    
    write.table(probe_gene_cormat, file = paste(Title,"Gene_wise_correlation_matrix_kendall.txt", sep = "_"), sep = "\t")
    rm(probe_gene_cormat)
    probe_gene_cormat= probe_gene_wise_cn_TACNAP_cor_mat(exp_data = smoothed_dataset_v1,cn_data=CN_dataset_genelevel_v1, cor_method = "spearman", no_cores = no_cores_v1)
    
    
    write.table(probe_gene_cormat, file = paste(Title,"Gene_wise_correlation_matrix_spearman.txt", sep = "_"), sep = "\t")
    
  }
  print(paste("Total time taken is", (proc.time()[3]-full_time)/60,"mins"))
}

probe_gene_wise_cn_TACNAP_cor_mat = function(exp_data, cn_data, cor_method, no_cores = detectCores())
{
  
  no_probe_gene = dim(cn_data)[1]
  
  factors_no_probe_gene = find_factors(no_probe_gene)
  
  no_of_cases_to_divide = factors_no_probe_gene[min(which(factors_no_probe_gene>=no_cores))]
  
  cor_mat_calculation = function(iter, no_of_cases_to_divide,exp_data, cn_data, cor_method1 )
  {
    nrows = (dim(cn_data)[1]/no_of_cases_to_divide)
    row_num = c(((iter-1)*nrows+1):(iter*nrows))
    cn_vec = as.matrix(cn_data[row_num,])
    
    if(nrows==1)
    {
      cor_vec = cor(cn_vec, t(exp_data), use = "pairwise.complete.obs", method = cor_method1)
    }else{
      cor_vec = cor(t(cn_vec), t(exp_data), use = "pairwise.complete.obs", method = cor_method1)
    }
    
    return(cor_vec)
  }
  
  if(no_of_cases_to_divide==no_probe_gene)
  {
    no_probe_gene = dim(cn_data)[1] - 1
    
    factors_no_probe_gene = find_factors(no_probe_gene)
    
    no_of_cases_to_divide = factors_no_probe_gene[min(which(factors_no_probe_gene>=no_cores))]
    
    cl <- makeCluster(no_cores, type = "FORK")
    
    x = parLapply(cl, 1:no_of_cases_to_divide,no_of_cases_to_divide = no_of_cases_to_divide, exp_data= exp_data,cn_data = cn_data,cor_method1 = cor_method,  cor_mat_calculation)
    stopCluster(cl)
    
    cn_vec = as.matrix(cn_data[dim(cn_data)[1],])
    
    z = cor(t(cn_vec), t(exp_data), use = "pairwise.complete.obs", method = cor_method)
    y = do.call("rbind", x)
    y = rbind(y,z)
  }else{
    cl <- makeCluster(no_cores, type = "FORK")
    
    x = parLapply(cl, 1:no_of_cases_to_divide,no_of_cases_to_divide = no_of_cases_to_divide, exp_data= exp_data,cn_data = cn_data,cor_method1 = cor_method,  cor_mat_calculation)
    stopCluster(cl)
    y = do.call("rbind", x)
    
  }
  
  
  
  rownames(y) = rownames(exp_data)
  colnames(y) = rownames(exp_data)
  return(y)
}

find_factors <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  return(factors)
}

permutation_testing_to_find_degr = function(dataset #Consensus estimated sources matrix
                                                     , data_with_chr_bp_mapping #genomic mapping file
                                                     ,mix_matrix #consensus mixing matrix 
                                                     ,intervals_for_gaus_kernel = seq(100000,2000000,by=10000) #Fix a set of possible standard deviation (sd) values as input. In the present study, input values were all integers from 10,000 to 2,000,000 with a gap of 10,000
                                                     ,probe_no_for_gaus_kernel = 10, #Obtain the optimal interval length for chromosome ch (〖oil〗_ch), where 〖oil〗_ch is the minimum of the il’s for which 〖quantile_nd〗_(ch,g,il)> 10
                                                     Title = "",
                                                     FDR = 0.1, #false discovery rate
                                                     CL = 0.5, #confidence level
                                                     state_deciding_cutoff = 0.95, #cutoff for Secondary indicator marks (sim)
                                                     set_seed = 123456, 
                                                     min_probesets = 5, #minimum number of genes/probesets to be present in the neighbourhood to consider a gene in EVR
                                                     full_dataset, #CES with genomic mapping variables
                                                     original_mapping_file, #probeset level genomic mapping file
                                                     gene_level_mapping_indicator = TRUE, #logical parameter to indicate if genelevel summary of DEGR is needed or not
                                                     collapse_function = "max", #function to collapse probeset-level data to gene-level
                                                     function_for_collapse, #other function to collapse probeset-level data to gene-level
                                                     TACNA_profiling = TRUE, #logical parameter to indicate if TACNA-profiles are to be obtained or not
                                                     file_main = file.path(getwd()) #directory to save the files
                                            ){
  full_time = proc.time()[3]
  row_num = dim(dataset)[1]
  if(row_num!=dim(data_with_chr_bp_mapping)[1])
  {
    print("two datasets are not of equal rows")
  }else{
    if(is.null(colnames(dataset)))
    {
      colnames(dataset) = paste("V",1:dim(dataset)[2],sep ="")
    }
    ##########################
    #Adjustment for Plot
    ###########################
    
    file_SCNA = file.path(file_main, paste("SCNA",
                                           "FDR",FDR,
                                           "CL", CL, "state_deciding_cutoff",state_deciding_cutoff,
                                           "probe_no_for_gaus_kernel",probe_no_for_gaus_kernel,
                                           "collapse_function", collapse_function,sep = "_"))
    dir.create(file_SCNA, showWarnings = FALSE)
    setwd(file_SCNA)
    
    chromosome_seq = c(0,cumsum(table(as.numeric(as.character(data_with_chr_bp_mapping$CHR_Mapping)))))
    
    data_with_chr_bp_mapping$BP_Mapping_v1 = data_with_chr_bp_mapping$BP_Mapping/10
    
    for(i in 2:max(data_with_chr_bp_mapping$CHR_Mapping))
    {
      data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]=data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]+data_with_chr_bp_mapping$BP_Mapping_v1[chromosome_seq[i]]
    }
    
    label_pos = sqldf("select distinct CHR_Mapping, avg(BP_Mapping_v1) as BP_Mapping from data_with_chr_bp_mapping group by 1 order by 1")
    
    
    
    ########################
    #Finding proper interval for sliding Gaussian Kernel
    #Choosing that interval where 10 or more no. of probesets 
    #corresponding to that chromosome are there in +/- 3*interval for 95% of the cases
    ########################
    quantile_5 = list()
    for( k in 1:max(data_with_chr_bp_mapping$CHR_Mapping))
    {
      quantile_5[[k]] = array(0,length(intervals_for_gaus_kernel))
      for(i in 1:length(intervals_for_gaus_kernel))
      {
        frequency_of_neigh_probes = NULL
        rows = list()
        for(j in (chromosome_seq[k]+1):(chromosome_seq[k+1]))
        {
          frequency_of_neigh_probes[j-chromosome_seq[k]] = length(which(abs(data_with_chr_bp_mapping$BP_Mapping[j] - data_with_chr_bp_mapping$BP_Mapping[c((chromosome_seq[k]+1):(chromosome_seq[k+1]))]) < 3*intervals_for_gaus_kernel[i]))
          
        }
        quantile_5[[k]][i] = quantile(frequency_of_neigh_probes,0.05)%/%1
        
        if(quantile_5[[k]][i]>=probe_no_for_gaus_kernel){
          print(paste("found",probe_no_for_gaus_kernel, "or more number of probe sets in 95% of the chromosome at interval", 
                      intervals_for_gaus_kernel[i],"for Chromosome",k))
          quantile_5[[k]][which(quantile_5[[k]]==0)] = NA
          break
        } 
      }
      
    }
    
    pdf(paste(Title,"graphs_to_select_interval_in_Gaussian_Kernel.pdf",sep = "_"))
    
    for(k in 1:max(data_with_chr_bp_mapping$CHR_Mapping))
    {
      plot(intervals_for_gaus_kernel, quantile_5[[k]],
           main = paste("Chromosome",k), 
           ylab = "5% quantile of No. of Probesets", 
           xlab = "intervals_for_gaus_kernel",pch=3,cex=0.1)
    }
    
    dev.off()
    
    ###########################
    #Creating the Density Matrix to get sliding gaussian Kernel
    ###########################
    
    density_matrix = matrix(0, row_num, row_num)
    
    for(j in 1:row_num)
    {
      rows[[j]] = which(data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[j])
      k = data_with_chr_bp_mapping$CHR_Mapping[j]
      s_x = as.numeric(data_with_chr_bp_mapping$BP_Mapping[rows[[j]]])
      s_mean = as.numeric(data_with_chr_bp_mapping$BP_Mapping[j]) 
      
      den_tnorm = dtruncnorm(s_x,
                             a = s_x[1],
                             b = s_x[length(s_x)],
                             mean = s_mean, 
                             sd = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[k]])))])
      
      density_matrix[j,rows[[j]]] = den_tnorm/sum(den_tnorm)
    }
    
    
    ###########################
    #Assigning intervals to a vector for each chromosome
    ###########################
    
    final_intervals = NULL
    
    for(chr_no in 1:max(data_with_chr_bp_mapping$CHR_Mapping)){
      final_intervals[chr_no] = intervals_for_gaus_kernel[max(which(!is.na(quantile_5[[chr_no]])))]
    }
    
    
    ###########################
    #Permutation Test for obtaining extreme valued region
    ###########################
    
    ################
    #Creating different directories for different outputs
    ################
    file_probelevel_bp = file.path(getwd(), "DEGR_sorted_by_base_pair_number_probelevel")
    file_probelevel_value = file.path(getwd(), "DEGR_sorted_by_CES_value_probelevel")
    dir.create(file_probelevel_bp, showWarnings = FALSE)
    dir.create(file_probelevel_value, showWarnings = FALSE)
    if(gene_level_mapping_indicator)
    {
      file_genelevel_bp = file.path(getwd(), "DEGR_sorted_by_base_pair_number_genelevel")
      file_genelevel_value = file.path(getwd(), "DEGR_sorted_by_CES_value_genelevel")
      dir.create(file_genelevel_bp, showWarnings = FALSE)
      dir.create(file_genelevel_value, showWarnings = FALSE)  
    }
    
    
    
    pdf(paste(Title,"All_Independent_CESs_plot.pdf",sep = "_"), height = 20,width= 30)
    extreme_valued_probesets_summary = NULL
    extreme_valued_section_summary = NULL
    for(CES_no in 1:dim(dataset)[2])
    {
      set.seed(set_seed+CES_no)
      time1 = proc.time()[3]
      # 1000 permutation of the non-smoothened weights of the gene-level CES are retained in a matrix p_CESMpx1000.
      N <- cbind(dataset[,CES_no],replicate(1000, sample(dataset[,CES_no])))
      # Obtain smoothened permuted CESM (sp_CESM) using the smoothing method explained in the step b.
      permute_1 = density_matrix%*%N
      # Sort the absolute values of the weights of all the columns of sp_CESMpx1000 in decreasing order. Also 3.	Sort the absolute values of the weights of the smoothened CES in decreasing order (sorted_CES).
      permute_1_v1 = apply(abs(permute_1),2,function(x) sort(x,decreasing=TRUE))
      
      # For every column of sp_CESMpx1000 (sp_CESi),
      # i.	For every weight of sorted_CES (sorted_CESj), obtain the number of weights of sp_CESi greater than sorted_CESj (f>j)
      # ii.	Obtain the optimal cutoff for the sp_CESi (oci) as the maximum value of the weights of sorted_CES for which f>j/j > 5%.
      
      cutoffs = array(0,1000)
      for(j in 2:1001)
      {
        for(i in 1:dim(permute_1_v1)[1])
        {
          if(length(which(permute_1_v1[,j]>permute_1_v1[i,1]))/i>=FDR)
          {
            cutoffs[j-1] = permute_1_v1[i,1]
            break
          }
        }
      }
      
      
      # Initial indicator marks (iim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way
      # i.	If smoothened_CESs > median(oc) then iims = 1
      # ii.	If smoothened_CESs < -median(oc) then iims = -1
      # iii.	Otherwise zero.
      
      indicator_ica_2 = ifelse(permute_1[,1]>quantile(cutoffs,CL),1,ifelse(permute_1[,1]< -quantile(cutoffs,CL),-1,0))
      
      # Smoothen iim (smoothened_iim) using the smoothing method described in step b.
      indicator_ica_3 = density_matrix%*%indicator_ica_2
      
      # Secondary indicator marks (sim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way:
      # i.	If smoothened_iims > 0.85 then sims = 1
      # ii.	If smoothened_iims < -0.85 then sims = -1
      # iii.	Otherwise zero.
      
      indicator_ica_4 = ifelse(indicator_ica_3>state_deciding_cutoff,1, ifelse(indicator_ica_3< -state_deciding_cutoff,-1,0))
      
      # Final indicator marks (fim) for every weight of the smoothened CES (smoothened_CESs) are obtained in the following way
      # Obtain the number of genes (ngs) mapped to the corresponding chromosome which have a distance from gene s in terms of base pair number < corresponding optimal interval length (oil) as described in step b.3.ii.
      # If ngs < 10 then fims = 0
      # Otherwise, fims = sims
      
      if(min_probesets>0)
      {
        state_of_the_weight_pos = which(abs(indicator_ica_4)>0)
        
        for(no in 1:length(state_of_the_weight_pos)){
          if(length(intersect(which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]] >= data_with_chr_bp_mapping$BP_Mapping[state_of_the_weight_pos[no]]-3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]]),
                              which(data_with_chr_bp_mapping$BP_Mapping[data_with_chr_bp_mapping$CHR_Mapping==data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]] <= data_with_chr_bp_mapping$BP_Mapping[state_of_the_weight_pos[no]]+3*final_intervals[data_with_chr_bp_mapping$CHR_Mapping[state_of_the_weight_pos[no]]])))<min_probesets){
            indicator_ica_4[state_of_the_weight_pos[no]] = 0
          }
        }
        
      }
      indicator_run = rle(as.vector(indicator_ica_4)) 
      
      limit = max(abs(min(dataset[,CES_no])),abs(max(dataset[,CES_no]))) 
      
      par(mar=c(3.1,3.1,3.1,3.1))
      plot(data_with_chr_bp_mapping$BP_Mapping_v1,
           dataset[,CES_no],
           pch=20,cex=0.4,ylim = c(-limit,limit),xaxt = "n",
           xlab = NA,
           ylab = NA, main = paste("CES", colnames(dataset)[CES_no], "FDR", FDR,  "CL", CL, "state_deciding_cutoff",state_deciding_cutoff), col=rgb(0, 0, 0, 0.2))
      mtext(side = 2, line = 2, 'CES Weight',cex = 0.8)
      mtext(side = 1, line = 2, paste("Chromosome No.s","of CES", colnames(dataset)[CES_no]),cex = 0.8)
      
      par(new = T)
      plot(data_with_chr_bp_mapping$BP_Mapping_v1,
           indicator_ica_4,xaxt = "n",
           pch=20,cex=0.4,
           ylim = c(-1,
                    1),
           xlab = NA,
           ylab = NA,axes = F, col=rgb(1, 0, 0, 0.2))
      axis(side = 4)
      mtext(side = 4, line = 2, 'State of the weight',cex = 0.8)
      axis(side = 1, at = label_pos[,2],labels = label_pos[,1])
      
      
      
      
      ######################
      #getting the extreme valued genomic positions
      ######################
      which_probe_extreme_valued = which(abs(indicator_ica_4)>0)
      if(length(which_probe_extreme_valued)>0){
        extreme_valued_probesets = as.data.frame(cbind(
          dataset[which_probe_extreme_valued,CES_no],
          row.names(dataset)[which_probe_extreme_valued],
          unlist(lapply(data_with_chr_bp_mapping$GENETITLE[which_probe_extreme_valued],as.character)),
          unlist(lapply(data_with_chr_bp_mapping$SYMBOL[which_probe_extreme_valued],as.character)),
          data_with_chr_bp_mapping$CHR_Mapping[which_probe_extreme_valued],
          data_with_chr_bp_mapping$BP_Mapping[which_probe_extreme_valued],
          data_with_chr_bp_mapping$BP_Mapping_v1[which_probe_extreme_valued],
          
          indicator_ica_4[which_probe_extreme_valued],
          
          rep(CES_no,length(which_probe_extreme_valued))))
        
        colnames(extreme_valued_probesets) = c( "Value", "ENTREZID", 'GENETITLE', 'SYMBOL',"CHR_Mapping", "BP_Mapping", "BP_Mapping_for_plot", "state_of_the_weight","CES")
        
        extreme_valued_probesets_summary = rbind(extreme_valued_probesets_summary,extreme_valued_probesets)
        
        run_vector = c(0, cumsum(indicator_run$lengths))
        
        extreme_valued_runs = which(abs(indicator_run$values)>0)
        
        for(num in 1:length(extreme_valued_runs)){
          
          if(gene_level_mapping_indicator){
            rows_vec = which(original_mapping_file$ENTREZID%in%data_with_chr_bp_mapping$ENTREZID[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1]])
            extreme_valued_section_gene_level = as.data.frame(cbind(dataset[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1],colnames(dataset)[CES_no]],data_with_chr_bp_mapping[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1],])) 
            
            colnames(extreme_valued_section_gene_level)[1] = "Value"
            
            if(indicator_run$values[extreme_valued_runs[num]]>0)
            {
              extreme_valued_section_gene_level_sort_by_value = sqldf("select * from extreme_valued_section_gene_level order by value DESC")
              
            }else{
              extreme_valued_section_gene_level_sort_by_value = sqldf("select * from extreme_valued_section_gene_level order by value")
              
            }
            
            write.table(extreme_valued_section_gene_level,
                        file = paste(file_genelevel_bp,paste(Title,"_CHR", extreme_valued_section_gene_level$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_gene_level$BP_Mapping[1]), "_BPEND",
                                                              round(extreme_valued_section_gene_level$BP_Mapping[dim(extreme_valued_section_gene_level)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
                        sep = "\t", row.names = FALSE)
            
            write.table(extreme_valued_section_gene_level_sort_by_value,
                        file = paste(file_genelevel_value,paste(Title,"_CHR", extreme_valued_section_gene_level_sort_by_value$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_gene_level_sort_by_value$BP_Mapping[1]), "_BPEND",
                                                                 round(extreme_valued_section_gene_level_sort_by_value$BP_Mapping[dim(extreme_valued_section_gene_level_sort_by_value)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
                        sep = "\t", row.names = FALSE)
            
          }else{
            rows_vec = which(original_mapping_file$PROBESET%in%data_with_chr_bp_mapping$PROBESET[(run_vector[extreme_valued_runs[num]]+1):run_vector[extreme_valued_runs[num]+1]])
            
          }
          
          
          
          
          extreme_valued_section = as.data.frame(full_dataset[rows_vec,c(colnames(dataset)[CES_no],'ENTREZID','GENETITLE', 'SYMBOL','CHR_Mapping','BP_Mapping')])
          
          extreme_valued_section$PROBESET = as.character(row.names(extreme_valued_section))
          
          
          colnames(extreme_valued_section)[1] = "Value"
          
          if(indicator_run$values[extreme_valued_runs[num]]>0)
          {
            extreme_valued_section_sort_by_value = sqldf("select * from extreme_valued_section order by value DESC")
          }else{
            extreme_valued_section_sort_by_value = sqldf("select * from extreme_valued_section order by value")
          }
          
          write.table(extreme_valued_section,
                      file = paste(file_probelevel_bp,paste(Title,"_CHR", extreme_valued_section$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section$BP_Mapping[1]), "_BPEND",
                                   round(extreme_valued_section$BP_Mapping[dim(extreme_valued_section)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
                      sep = "\t", row.names = FALSE)
          
          write.table(extreme_valued_section_sort_by_value,
                      file = paste(file_probelevel_value,paste(Title,"_CHR", extreme_valued_section_sort_by_value$CHR_Mapping[1],"_BPSTART",round(extreme_valued_section_sort_by_value$BP_Mapping[1]), "_BPEND",
                                                            round(extreme_valued_section_sort_by_value$BP_Mapping[dim(extreme_valued_section_sort_by_value)[1]]),"_CES",colnames(dataset)[CES_no], ".txt", sep = ""),sep = "/"),
                      sep = "\t", row.names = FALSE)
          CES = rep(CES_no,dim(extreme_valued_section)[1])
          
          extreme_valued_section_summary = rbind(extreme_valued_section_summary,cbind(extreme_valued_section,CES))
          
            }
        
        
        print(paste("Extreme valued identified in some genomic postiions of CES",colnames(dataset)[CES_no]))
      }else{
        print(paste("No Genomic position has got Extreme valued region for CES",colnames(dataset)[CES_no]))
      }
      print(paste("Processing of CES", colnames(dataset)[CES_no], "FDR", FDR,  "CL", CL, "state_deciding_cutoff",state_deciding_cutoff, "is done"))
      
      print(paste("Time taken for this iteration is",round((proc.time()[3]-time1)/60,3),"mins"))
      
    }
    dev.off()
    
    write.table(extreme_valued_probesets_summary, file = paste("Genelevel",Title,
                                                         paste("All_CESs_extreme_valued_details",  "FDR", FDR, 
                                                               "CL",CL, "state_deciding_cutoff",
                                                               state_deciding_cutoff, sep = "_"),".txt", sep = "_"), sep = "\t", row.names = FALSE)
    
    write.table(extreme_valued_section_summary, file = paste("Probelevel",Title,
                                                                   paste("All_CESs_extreme_valued_details",  "FDR", FDR, 
                                                                         "CL",CL, "state_deciding_cutoff",
                                                                         state_deciding_cutoff, sep = "_"),".txt", sep = "_"), sep = "\t", row.names = FALSE)
    
    extreme_valued_probesets_summary = as.matrix(extreme_valued_probesets_summary)
    pdf(paste(Title,
              paste("Extreme_valued_Region_all_CESs_combined",  "FDR", FDR, 
                    "CL",CL, "state_deciding_cutoff",
                    state_deciding_cutoff, sep = "_"),".pdf", sep = "_"), height = 20,width= 30)
    
    par(mar=c(3.1,3.1,3.1,3.1))
    plot(extreme_valued_probesets_summary[,'BP_Mapping_for_plot'],
         extreme_valued_probesets_summary[,'Value'],
         pch=20,cex=0.4,
         xlim = c(min(data_with_chr_bp_mapping$BP_Mapping_v1), max(data_with_chr_bp_mapping$BP_Mapping_v1)),
         xaxt = "n",
         xlab = NA,
         ylab = NA, main = paste("Extreme valued Region by all CES combined weight"), col=rgb(0, 0, 0, 0.2))
    mtext(side = 2, line = 2, 'CES Weight',cex = 0.8)
    mtext(side = 1, line = 2, paste("Chromosome No.s of CES"),cex = 0.8)
    axis(side = 1, at = label_pos[,2],labels = label_pos[,1])
    
    
    plot(extreme_valued_probesets_summary[,'BP_Mapping_for_plot'],
         extreme_valued_probesets_summary[,'state_of_the_weight'],
         pch=20,cex=0.4,
         xlim = c(min(data_with_chr_bp_mapping$BP_Mapping_v1), max(data_with_chr_bp_mapping$BP_Mapping_v1)),
         xaxt = "n",
         xlab = NA,
         ylab = NA, main = paste("Extreme valued Region by all CES combined State"), col=rgb(0, 0, 0, 0.2))
    mtext(side = 2, line = 2, 'State of the weight',cex = 0.8)
    mtext(side = 1, line = 2, paste("Chromosome No.s of CES"),cex = 0.8)
    axis(side = 1, at = label_pos[,2],labels = label_pos[,1])
    
    dev.off()
    
    
    if(TACNA_profiling){
      
      file_TACNA_profile = file.path(file_SCNA, "TACNA_Profile_results")
      dir.create(file_TACNA_profile, showWarnings = FALSE)
      setwd(file_TACNA_profile)
      
      extreme_valued_comps = unique(as.numeric(as.character(extreme_valued_section_summary$CES)))
      
      updated_dataset = matrix(0, dim(full_dataset)[1],dim(full_dataset)[2])
      rownames(updated_dataset) = rownames(full_dataset)
      colnames(updated_dataset) = colnames(full_dataset)
      
      row_num_match = NULL
      for(alteration_probes in 1:dim(extreme_valued_section_summary)[1])
      {
        row_num_match[alteration_probes] = which(rownames(full_dataset)==extreme_valued_section_summary$PROBESET[alteration_probes])
        
      }
      
      for(ics in extreme_valued_comps)
      {
        rows = which(extreme_valued_section_summary$CES==ics)
        updated_dataset[row_num_match[rows],ics] = as.numeric(as.character(extreme_valued_section_summary$Value[rows]))
      }
      
      updated_dataset = updated_dataset[,extreme_valued_comps]
      
      updated_mix_matrix = mix_matrix[extreme_valued_comps,]
      
      updated_expression_data = as.matrix(updated_dataset)%*%as.matrix(updated_mix_matrix)
      
      rownames(updated_expression_data) = rownames(updated_dataset)
      colnames(updated_expression_data) = colnames(updated_mix_matrix)
      
      data_with_chr_bp_mapping = original_mapping_file
      
      matching_vec = match(data_with_chr_bp_mapping$PROBESET,rownames(full_dataset))
      
      data_with_chr_bp_mapping = data_with_chr_bp_mapping[!(is.na(matching_vec)),]
      
      matching_vec = matching_vec[!(is.na(matching_vec))]
      
      updated_expression_data_v1 = cbind(as.data.frame(updated_expression_data[matching_vec,]),data_with_chr_bp_mapping)
      write.table(updated_expression_data, file = paste(Title,
                                                        paste("TACNA_profiled_expression_data",  "FDR", FDR, 
                                                              "CL",CL, "state_deciding_cutoff",
                                                              state_deciding_cutoff,sep = "_"),".txt",sep = "_"), sep = "\t")
      
      updated_expression_data = as.data.frame(updated_expression_data[matching_vec,])
      chromosome_seq = c(0,cumsum(table(as.numeric(as.character(data_with_chr_bp_mapping$CHR_Mapping)))))
      
      data_with_chr_bp_mapping$BP_Mapping_v1 = data_with_chr_bp_mapping$BP_Mapping/10
      
      
      for(i in 2:max(data_with_chr_bp_mapping$CHR_Mapping))
      {
        data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]=data_with_chr_bp_mapping$BP_Mapping_v1[(chromosome_seq[i]+1):chromosome_seq[i+1]]+data_with_chr_bp_mapping$BP_Mapping_v1[chromosome_seq[i]]
      }
      
      label_pos = sqldf("select distinct CHR_Mapping, avg(BP_Mapping_v1) as BP_Mapping from data_with_chr_bp_mapping group by 1 order by 1")
      
      
      limit = 7
      pdf(paste(Title,
                paste("TACNA_profiled_expression_data",  "FDR", FDR, 
                      "CL",CL, "state_deciding_cutoff",
                      state_deciding_cutoff,sep = "_"),".pdf",sep = "_"))
      for(samp_no in 1:dim(updated_expression_data)[2]){
        
        par(mar=c(3.1,3.1,3.1,3.1))
        plot(data_with_chr_bp_mapping$BP_Mapping_v1,
             updated_expression_data[,samp_no],
             pch=20,cex=0.1,ylim = c(-limit,limit),xaxt = "n",
             xlab = NA,
             ylab = NA, main = paste("Sample", colnames(updated_expression_data)[samp_no]), col=rgb(0, 0, 0, 0.2))
        mtext(side = 2, line = 2, 'TACNA profiled expression',cex = 0.8)
        mtext(side = 1, line = 2, paste("Chromosome No.s","of Sample", colnames(updated_expression_data)[samp_no]),cex = 0.8)
        axis(side = 1, at = label_pos[,2],labels = label_pos[,1],pch=1,cex=0.12)
        
      }
      dev.off()
      
      
    }
    
    
    
  }
  print(paste("Time taken for the whole state detection analysis is",round((proc.time()[3]-full_time)/60,3),"mins"))
}
