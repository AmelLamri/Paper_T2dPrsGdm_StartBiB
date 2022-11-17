#!/usr/bin/env Rscript
# Written by Amel Lamri
# Last modified on 16/11/2022


###############
# DESCRIPTION #
###############  

  # #This script derives a T2D PRS in START and BiB using LDpred2 .
  # Large portions of this scripts follow the instructions & code described in
  # https://privefl.github.io/bigsnpr/articles/LDpred2.html or from by F. Prive
  #For more information contact Amel Lamri (lamria@mcmaster.ca)

# Load Libraries & set options parameters

  library(PGPmisc)   # my functions
  library(bigsnpr) # Load packages bigsnpr and bigstatsr
  library(runonce)
  library(rsample)   # for resampling procedures
  library(ggplot2) 
  library(sas7bdat) #to read SAS files
  library(dplyr) #for piping
  library(DescTools)#for cstat
  library(fmsb) #to calculate the pseudo R2

# Set option, Seed and define paths

  options(bigstatsr.check.parallel.blas = FALSE)
  options(default.nproc.blas = NULL)
  set.seed(131)# Set seed

  NCORES <- nb_cores() 
  lmd <- "220609"
  trait<- "Binary"

  # Paths 
  avr_p="/home/lamria/avr"
  srv_p="/home/lamria/srvstore"
  avrd_p<- paste0(avr_p,"/Data")
  stp_p 	<- paste0(avrd_p, "/START/Phenotypes")
  ms_p<- paste0(avr_p, "/Projects/byStudy/MixedStudies")
  r0_p <- paste0(ms_p, "/START_BiB_1KG_Mahajan2014/202107_ListCommonSNPs/Results") #results folder
  m_p <- paste0(srv_p, "/Merged_Studies/BiB_START_1KG_Mahajan2014/202108_Derive_T2Dprs")
  g_p <- paste0(m_p, "/Genotypes")
  r_p <- paste0(m_p, "/Results")
  p_p <- paste0(m_p, "/Phenotypes")
  o_p <- paste0(ms_p, "/START_BiB/202108_Derive_LDpred2_T2D_PRS/Results")
  s_p <- paste0(m_p, "/Scripts")
  s_fn <- "202107_START_BiB_LDpred2_T2D_PRS_S2_Derive_PRS"
  bo_fnpe <- paste0(o_p, "/BiB_T2D_PRS_Mah14_LDpred2BestGrid_", lmd, ".csv")
  so_fnpe <- paste0(o_p, "/START_T2D_PRS_Mah14_LDpred2BestGrid_", lmd, ".csv")
  sbo_fnpe <- paste0(o_p, "/START_BiB_T2D_PRS_Mah14_LDpred2BestGrid_", lmd, ".csv")
  # tmp <- tempfile(tmpdir = paste0(r_p, "/tmp-data"))
  # on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)


# define P values to be used in the grid search

  grid2p <- c(5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 0.05, seq(0.1,1,0.1)) 
  grid2p_c <- c("5e-8", "5e-7", "5e-6", "5e-5", "5e-4", "5e-3", as.character(c(0.05, seq(0.1,1,0.1)))) 

# Paths Part II 
  for (i in length(grid2p)){
    gi_p <- paste0(g_p, "/", grid2p[i])
    genoifnp <- paste0(gi_p, "/START_BiB_1KG_SNPs_in_START_BiB_Mah14_1KG_P_se_", grid2p[i])
    genoifnp2 <- paste0(gi_p, "/START_BiB_SNPs_in_START_BiB_Mah14_1KG_P_se_", grid2p[i], "_noMiss")
  }



# Load and prep genotypes of START BiB 1KG merged 

  # Read from bed/bim/fam, it generates .bk and .rds files.
  if(!file.exists(paste0(genoifnp2, ".rds"))) snp_readBed(paste0(genoifnp2, ".bed"))
  
  # Attach the "bigSNP" object
  obj.bigSNP <- snp_attach(paste0(genoifnp2,".rds"))
  str(obj.bigSNP, max.level = 2, strict.width = "cut")
  
  # Get aliases for useful slots
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  y   <- obj.bigSNP$fam$affection

# Read and prep extrenal summary statistics

  common <- read.table(paste0(r0_p, "/START_BiB_Mahajan2014_1KG_CommonSNPsOnly_infos.tsv"), head=T, string=F, sep="\t")
  common <- common[!is.na(common$Direction_SA),]
  
  common$Beta <- log(common$OR)
  #common$OR_seL=-(common$OR_95L-common$OR)/1.96
  #common$OR_seU=(common$OR_95U-common$OR)/1.96
  #common$OR_se<-mean(common$OR_seL,common$OR_seU)
  
  common$beta_L95 <- log(common$"OR_95L")
  common$beta_U95 <- log(common$"OR_95U")
  common$beta_se1 <- -(common$beta_L95-common$Beta)/1.959964
  common$beta_se2 <- (common$beta_U95-common$Beta)/1.959964
  common$beta_se <- rowMeans(common[,c("beta_se1", "beta_se2")])
  
  sumstats <- common[c("RSID_s_impd", "CHR","POS", "OtherAllele", "RiskAllele", "Beta", "beta_se", "N", "P")]
  colnames(sumstats) <-  c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "N", "p") # In the previous version of the script, I had named OtherAllele as a0 and RiskAllele as a0, which is the opposite of what I should have done since https://privefl.github.io/bigsnpr/articles/LDpred2.html clearly states that "a0" (reference allele) and "a1" (derived allele), as a results, the direction of my PRS was reversed ( association was in opposite direction), here, I fixed this.
   
  str(sumstats) #Visualize the data
  



# Open SNPs posoitions and infos 

  info <- readRDS(runonce::download_file("https://ndownloader.figshare.com/files/25503788", dir = paste0(r_p, "/tmp-data"), fname = "map_hm3_ldpred2.rds"))
  str(info)# Here I used the file for build GRCh37 / hg19)



# Get samples FIDs and IIDs 

  fam <- read.table(paste0(genoifnp2, ".fam"), head=F, string=F)[,1:2]
  fam$order<-1:nrow(fam)
  colnames(fam)[1:2] <- c("FID", "IID")
  fams <- read.table(paste0(p_p, "/START_BiB_1KG_FID_IIDs_Study_infos.txt"), head=T, string=F)
  fam <- merge(fam, fams, all.x=T, all.y=F)


# Open and prep phenotypes 
  
  
  # BiB 
  
  
    bibgdm <- read.csv(paste0(avrd_p, "/Born_in_Bradford/Phenotypes/GDM_RelatedVars/Bib_Mothers_GDMsaBiBwBfnT2DMaxVal_GDMsaBiBwBfnT2DPregNs_GDMsaBiBwBfnT2DTimeMaxVal_MotherLevel_20210824.csv"), head=T, string=F)
    
    bibids <-  read.table(paste0(avrd_p, "/Born_in_Bradford/Phenotypes/BiB_FID_IID_PregnancyID_MotherID_PersonID_correspondance_PregnancyLevel20201228.txt"), head=T, string=F)
    
    bibgdm2 <- merge(bibgdm[,c( "MotherID" ,"GDMsaBiBwBfnT2DMaxVal")], bibids[,c("PersonID", "IID")], by.x="MotherID", by.y="PersonID",all.x=F, all.y=F)
    
    bibgdm3 <- merge(bibgdm2,fam, all=F)
    nrow(bibgdm3) ==nrow(bibgdm2)
  
  
  
  
  
  # START 
    st1 <- read.sas7bdat(paste0(avrd_p, "/START/Phenotypes/Raw_or_Received/gdm_all_2019.sas7bdat")) #gdm and AUC data
    st1$Study <- "START"
    
    # Recode phentoype variables of interest 
    st1$IID <- st1$pmkitn 
    st1[which(st1$GDM_BiBp == 2), "GdmSaBiBp01"] <- 1
    st1[which(st1$GDM_BiBp == 1), "GdmSaBiBp01"] <- 0
    st1[which(is.nan(st1$GDM_BiBp)), "GdmSaBiBp01"] <- NA
    
    table(st1$IID %in% fam$IID) 
    
    st2 <- merge(fam, st1[,c("id", "IID", "GdmSaBiBp01")], all=F)
    table(st2$GdmSaBiBp01, useNA="ifany")
    
    #st2$GDM_BiBp <- NULL

# Merge FAM files infos 

  fam2 <- merge(bibgdm3, st2, all=T)
  colnames(fam)[colnames(fam) %in% colnames(fam2)] 
  fam3 <- merge(fam2, fam[which(!fam$IID %in% fam2$IID), ], all=T)
  fam4 <-fam3[order(fam3$order), ]
  table(fam4$order == 1:nrow(fam4))


# Create a homogenous GDM variable between START and BiB 
 
  fam4[which(fam4$Study == "START"), "GDM"] <- fam4[which(fam4$Study == "START"), "GdmSaBiBp01"]
  fam4[which(fam4$Study == "BiB"), "GDM"] <- fam4[which(fam4$Study == "BiB"), "GDMsaBiBwBfnT2DMaxVal"]
  table(fam4$GDM, useNA="ifany")
  table(fam4$GDM,fam4$Study, useNA="ifany")
  y<- fam4$GDM

# Splitting data to training and validation sets

  # I'll do this using Stratified sampling
  # outputs from this section are : 
  #ind.val which is the index of individauls to be included in the trainaing set
  #ind.test which is the index of individauls to be included in the testing set
  
  
  #only keep samples with GDM CC status 
  
  famwgdm<-fam4[!is.na(fam4$GDM),]
  split_strat  <- initial_split(famwgdm, prop = 0.7,strata = "Study")
  
  
  train_strat  <- training(split_strat)
  test_strat   <- testing(split_strat)
  
  ind.val <- which(fam4$"IID" %in% train_strat$IID) #sample( nrow(G), 350)
  ind.test <- which(fam4$"IID" %in% test_strat$IID) 
  
  #prop.table(table(fam$IID %in% train_strat$IID))
  #prop.table(table(fam$IID %in% test_strat$IID))



# Calculate effective sample size 

  # This wil be useed later when deriving the PRS 
  # From the summary statistics, you need to get "beta", "beta_se" (standard errors), and "n_eff" (effective sample size per variant for GWAS with logistic regression, and just total sample size for continuous traits).

  # Here, I don't have the number of cases and controls per SNP for Mahajan 2014, so I used the overal % of cases/controls mentionned in the aBstract of the manuscript , which is 21491 CASES /( 21491 CASES + 83964 CTRL) = 0.2037931,  ~ 20% 
 
  trait <- "binary"
  if(trait %in% c("Binary", "binary")){
    sumstats$n_case <- round((sumstats$N* 20)/100,0)
    sumstats$n_control <- sumstats$N - sumstats$n_case
    sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
    sumstats$n_case <- sumstats$n_control <- NULL
  }else if (trait %in% c("Continuous", "continuous")) {
    sumstats$n_eff <- sumstats$N
    }

# Matching variants between genotype data and summary statistics
    
  #To match variants contained in genotype data and summary statistics, the variables "chr" (chromosome number), "pos" (genetic position), "a0" (reference allele) and "a1" (derived allele) should be available in the summary statistics and in the genotype data. These 4 variables are used to match variants between the two data frames. 
  #head(sumstats[,c("N", "n_control","n_case")])
  
  # extract the SNP information from the genotype
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  df_beta <- snp_match(sumstats, map)


# Computing LDpred2 scores 
  # Compute Correlations

  # First, we need to compute correlations between variants. The authors recommend to use a window size of 3 cM (see paper).

  # POS2a <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

  POS2 <- obj.bigSNP$map$genetic.dist
  tmp <- tempfile(tmpdir = paste0(r_p, "/tmp-data_v3"))
  for (chr in 1:22) {
    
     print(chr)
    
    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'G'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    
    corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                     infos.pos = POS2[ind.chr2], ncores = NCORES)
    
    if (chr == 1) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp, compact = TRUE)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
      }
    }
  
  file.size(corr$sbk) / 1024^3  # file size in GB


  # Run LDpred2-inf: infinitesimal model

  ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))
  
  h2_est <- ldsc[["h2"]]
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
  pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
  cor(pred_inf, y[ind.test], use="complete.obs")

  # Calculate the null R2

  #this comes from  https://choishingwan.github.io/PRS-Tutorial/ldpred/
  # For binary traits : 
    # IMPORTANT Scripts for binary trait analysis only serve as a reference as we have not simulate any binary traits. In addition, Nagelkerke R2  is biased when there are ascertainment of samples. For more information, please refer to this paper

    # Reformat the phenotype file such that y is of the same order as the 
    # sample ordering in the genotype file
    #y <- pheno[fam.order, on = c("FID", "IID")]
    # Calculate the null R2
    # use glm for binary trait 
    # (will also need the fmsb package to calculate the pseudo R2)
    #null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    #    paste0("Height~Sex+", .) %>%
    #    as.formula %>%
    #    glm(., data = y, family=binomial) %>%
    #    summary
    #null.r2 <- fmsb::NagelkerkeR2(null.model)


  # For quantitiative traits 
  
  
    # Reformat the phenotype file such that y is of the same order as the 
    # sample ordering in the genotype file
    #y <- pheno[fam.order, on = c("FID", "IID")]
    # Calculate the null R2
    # use glm for binary trait 
    # (will also need the fmsb package to calculate the pseudo R2)
    #null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    #    paste0("Height~Sex+", .) %>%
    #    as.formula %>%
    #    lm(., data = y) %>%
    #    summary
    #null.r2 <- null.model$r.squared










  # Run LDpred2 grid

  
  h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
  p_seq <- unique(signif(c(seq_log(5e-8, 0.001, length.out = 8), seq_log(0.001, 0.1,15), seq(0.1, 0.7, 0.1)),2))
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
  
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
  
  
  pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
  
  params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
    if (all(is.na(x))) return(NA)
      #if(trait %in% c("Continuous", "continuous")){ 
      # summary(lm(y[ind.val] ~ x))$coef["x", 3] #to be used for continuous trait 
      #} else if (trait %in% c("Binary", "binary")){
        summary(glm(y[ind.val] ~ x, family = binomial(link='logit')))$coef["x", 3]
  	 # null.r2 <- fmsb::NagelkerkeR2(null.model)
      #}
      }
    )
  
  
  
  params$cstat <- apply(pred_grid[ind.val, ], 2, function(x) {
    if (all(is.na(x))) return(NA)
      #if(trait %in% c("Continuous", "continuous")){ 
      # summary(lm(y[ind.val] ~ x))$coef["x", 3] #to be used for continuous trait 
      #} else if (trait %in% c("Binary", "binary")){
        Cstat(glm(y[ind.val] ~ x, family = binomial(link='logit')))
  	 # null.r2 <- fmsb::NagelkerkeR2(null.model)
      #}
      }
    )

   # Note that missing values represent models that diverged substantially.
  
  # Plot results 
    
    
    jpeg(paste0(r_p, "/home/lamria/ldpred2_v2_grid_plot.jpeg"),  width = 960, height = 480, )
    myplot <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
      theme_bigstatsr() +
      geom_point() +
      geom_line() +
      scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
      facet_wrap(~ sparse, labeller = label_both) +
      labs(y = "GLM Z-Score ", color = "h2") +
      theme(legend.position = "top", panel.spacing = unit(1, "lines"))
    print(myplot)
    dev.off()
    
    
    jpeg(paste0(o_p, "/home/lamria/ldpred2_v2_grid_Cstat.jpeg"),  width = 960, height = 480, )
    myplot <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
      theme_bigstatsr() +
      geom_point() +
      geom_line() +
      scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
      facet_wrap(~ sparse, labeller = label_both) +
      labs(y = "GLM Cstat", color = "h2") +
      theme(legend.position = "top", panel.spacing = unit(1, "lines"))
    print(myplot)
    dev.off()
  

  # visualiwe results
  
  params %>%
    mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
    arrange(score) %>%
    mutate_at(c("score", "sparsity"), round, digits = 3) %>%
    slice(1:10)

  
  # Extract best PRS from grid search and run predictions 
  
    best_beta_grid <- params %>%
      mutate(id = row_number()) %>%
      # filter(sparse) %>% 
      arrange(desc(score)) %>%
      slice(1) %>%
      pull(id) %>% 
      beta_grid[, .]
    
    pred <- big_prodVec(G, best_beta_grid, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
    cor(pred, y[ind.test], use="complete.obs")



# Merge best score with FIDs and IIDs and clean 

  bestGridScore_d <- cbind(fam[,c("FID", "IID")], pred_grid[,65])
  # rename PRS 
  colnames(bestGridScore_d)[3] <- "T2dMah14Prs" 
  # save raw PRS
  write.csv(bestGridScore_d, paste0(o_p, "/START_BiB_LDpred2_BestGrid_PRS_oldScoreFilppedDirection20220610.csv"), row.names=F, quote=F)

# Standardize PRS and derive PRS categories
  # START 
    
    st1[is.nan(st1$IID) , "IID"] <- NA
    mstscore <- merge(st1, bestGridScoren_d, all=F)
    
    #PRS Z score 
    mstscore$"T2dMah14PrsZ" <- scale(mstscore$"T2dMah14Prs") 
    #PRS tertiles
    mstscore$T2dMah14PrsTer <- mstscore %>% mutate(quantile = ntile(mstscore$T2dMah14PrsZ, 3))  %>% .$quantile %>% factor() %>% as.character 
    table(mstscore$T2dMah14PrsTer)
    table(mstscore$"T2dMah14PrsTer", useNA="ifany")
    round(prop.table(table(mstscore$"T2dMah14PrsTer", useNA="ifany"))*100,1)
    #PRS T1,2,3 dummy coded
    mstscore[,"T2dMah14PrsT1vT2"] <- recode(mstscore$T2dMah14PrsTer, "1"=0, "2"=1) %>% as.character(na_if(mstscore$T2dMah14PrsTer, 3))
    table(mstscore$T2dMah14PrsT1vT2, useNA="ifany")
    mstscore[,"T2dMah14PrsT1vT3"] <- recode(mstscore$T2dMah14PrsTer, "1"=0, "3"=1) %>% as.character(na_if(mstscore$T2dMah14PrsTer, 2))
    table(mstscore$T2dMah14PrsT1vT3, useNA="ifany")
    
    #PRS T1 v T2+T3 and  T1 + T2 v. T3
    mstscore[,"T2dMah14PrsT1vT23"] <- recode(mstscore$T2dMah14PrsTer, "1"=0, "2"=1, "3"=1) 
    table(mstscore$T2dMah14PrsT1vT23, useNA="ifany")
    mstscore[,"T2dMah14PrsT12vT3"] <- recode(mstscore$T2dMah14PrsTer, "1"=0, "2"=0, "3"=1) 
    table(mstscore$T2dMah14PrsT12vT3, useNA="ifany")
    
    #PRS top decile vs bottom 90% 
    mstscore$T2dMah14PrsZtop10 <- ifelse(mstscore$T2dMah14PrsZ < quantile(mstscore$T2dMah14PrsZ, prob=9/10, na.rm=T), 0,1)
   
  
  
  # BiB 
  
    mbscore<- merge(bibgdm3, bestGridScore_d, all=F)
    
    mbscore$"T2dMah14PrsZ" <- scale(mbscore$"T2dMah14Prs") 
    mbscore$T2dMah14PrsTer <- mbscore %>% mutate(quantile = ntile(mbscore$T2dMah14PrsZ, 3))  %>% .$quantile %>% factor() %>% as.character 
    table(mbscore$T2dMah14PrsTer)
    table(mbscore$"T2dMah14PrsTer", useNA="ifany")
    round(prop.table(table(mbscore$"T2dMah14PrsTer", useNA="ifany"))*100,1)
    #PRS T1,2,3 dummy coded
    mbscore[,"T2dMah14PrsT1vT2"] <- recode(mbscore$T2dMah14PrsTer, "1"=0, "2"=1) %>% as.character(na_if(mbscore$T2dMah14PrsTer, 3))
    mbscore[,"T2dMah14PrsT1vT3"] <- recode(mbscore$T2dMah14PrsTer, "1"=0, "3"=1) %>% as.character(na_if(mbscore$T2dMah14PrsTer, 2))
    #PRS T1 v T2+T3 and  T1 + T2 v. T3
    mbscore[,"T2dMah14PrsT1vT23"] <- recode(mbscore$T2dMah14PrsTer, "1"=0, "2"=1, "3"=1) 
    mbscore[,"T2dMah14PrsT12vT3"] <- recode(mbscore$T2dMah14PrsTer, "1"=0, "2"=0, "3"=1) 
    
    #PRS top decile vs bottom 90% 
    mbscore$T2dMah14PrsZtop10 <- ifelse(mbscore$T2dMah14PrsZ < quantile(mbscore$T2dMah14PrsZ, prob=9/10, na.rm=T), 0,1)
 


# Save outputs  

  # START and BIB PRSs 
    write.csv(bestGridScoren_d, sbo_fnpe, row.names=F, quote=F)
  
  # START 
    write.csv(mstscore[,c("id","pmkitn", "T2dMah14Prs", "T2dMah14PrsZ","T2dMah14PrsTer", "T2dMah14PrsT1vT2", "T2dMah14PrsT1vT3", "T2dMah14PrsT1vT23", "T2dMah14PrsT12vT3", "T2dMah14PrsZtop10")], so_fnpe, row.names=F, quote=F)
  #BiB
    write.csv(mbscore[,c("MotherID", "T2dMah14Prs", "T2dMah14PrsZ","T2dMah14PrsTer", "T2dMah14PrsT1vT2", "T2dMah14PrsT1vT3", "T2dMah14PrsT1vT23", "T2dMah14PrsT12vT3", "T2dMah14PrsZtop10")], bo_fnpe, row.names=F, quote=F)

# derive .R script from .Rmd file
# knitr::purl(paste0(s_p, "/", s_fn, ".Rmd"))
