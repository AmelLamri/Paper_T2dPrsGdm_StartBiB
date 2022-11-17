#!/usr/bin/env Rscript
# Written by Amel Lamri (lamria@mcmaster.ca)
# Last modified on 16/11/2022

#################
## DESCRIPTION ##
#################

  # This script prepares all the files necessary to run LDpred2 using a list of common
  # SNPs between START and BiB genotypes and Mahajan 2014.
  # This is a subscript of: 2211_START_BiB_T2D_PRS_GDM_Main.sh
  # For more information, contact Amel Lamri (lamria@mcmaster.ca)


# Load required libraries 

  library(PGPmisc)

# Define paths and file names 

  h_p <- "/home/lamria"
  pdrv_p <- paste0(h_p, "/pdrive") 
  avr_p <- paste0(h_p, "/avr")
  avrd_p <- paste0(avr_p, "/Data")
  srv_p <- paste0(h_p, "/srvstore")
  
  m_p <- paste0(h_p, "/avr/Projects/MixedStudies/START_BiB_1KG_Mahajan2014")
  r0_p <- paste0(m_p, "/202107_ListCommonSNPs//Results") # results folder 
  r_p <- paste0(r_p, "/202108_Derive_LDpred2_T2D_PRS/Results") 
  r5_p <- paste0(srv_p,"/Merged_Studies/BiB_START_1KG_Mahajan2014/202108_Derive_T2Dprs")



# Create P-val grid values

  grid2p <- c(5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 0.05, seq(0.1,1,0.1)) 
  grid2p_c <- c("5e-8", "5e-7", "5e-6", "5e-5", "5e-4", "5e-3", as.character(
    c(0.05, seq(0.1,1,0.1)))) 
 

# Open list of common SNPS between START and BiB and Mahajan 2014
  common <- read.table(paste0(r0_p, "/START_BiB_Mahajan2014_1KG_CommonSNPsOnly_infos.tsv"), head=T, string=F, sep="\t")


# Prep FamilyIDs 

  fs <- read.table(paste0( srv_p,"/Nutrigen/Genotypes/START/GWAS_HCE_ICE/Filtered/V1/InfoScore0.7/START_GW_HCE_ICE_hg19_updated_phased_imputed_chr22_filtered_info0.7_Miss0.05_CallThres0_UpdatedAllelesAndIDsAndParents.fam"), head=F, string=F)[,1:2]
  fb<- read.table(paste0(srv_p, "/Born_in_Bradford/Genotypes/202101_1KG_Imputation/BiB_SouthAsian_Mothers_Imputed/CoreExome_GSA/Chr22/SNP_Subsets/Info0.7/BiB_SouthAsian_Mothers_CoreExomeGSA_Imputed_1KG_Chr22_Info0.7_geno0.95_hct0.5.fam"), head=F, string=F)[,1:2]
  f1<-read.table(paste0(avrd_p, "/1000Genomes/20130502_Release/Genotypes/PLINK/Sample_Subsets/South_Asians/1000G_SA_chr22.fam"), head=F, string=F)[,1:2]
  
  fs$Study<-"START"
  fb$Study<-"BiB"
  f1$Study<-"1KG"
  fams<-rbind(fs, fb, f1)
  colnames(fams)[1:2] <- c("FID", "IID")
  write.table(fams, paste0(srv_p, "/Merged_Studies/BiB_START_1KG_Mahajan2014/202108_Derive_T2Dprs/Phenotypes/START_BiB_1KG_FID_IIDs_Study_infos.txt") , col.names=T, row.names=F, quote=F)


# Compare alleles between studies #


  # BIB GSA vs BiB HCE 
    CheckSnpAlleles(common,"A1_b_gtypd_gsa", "A2_b_gtypd_gsa", "A1_b_gtypd_hce", "A2_b_gtypd_hce", "chrpos") # 0 SNPs with discordant alleles
  
  # BIB HCE vs START HCE 
    CheckSnpAlleles(common,"A1_b_gtypd_gsa", "A2_b_gtypd_gsa", "A1_s_impd", "A2_s_impd", "chrpos")# 0 SNPs with discordant alleles
  
  
  # BIB HCE GSA vs Mahajan 
    discordAll1 <- CheckSnpAlleles(common,"A1_b_impd_hceGsa", "A2_b_impd_hceGsa", "RiskAllele", "OtherAllele", "chrpos") #  7446 SNPs with discordant alleles

  
  # BIB HCE GSA vs START  
    CheckSnpAlleles(common2,"A1_b_impd_hceGsa", "A2_b_impd_hceGsa", "A1_s_impd", "A2_s_impd", "chrpos") #  0 SNPs with discordant alleles
  
  # BIB HCE GSA vs 1KG  
    CheckSnpAlleles(common2,"A1_b_impd_hceGsa", "A2_b_impd_hceGsa", "A1_kg", "A2_kg", "chrpos") #  0 SNPs with discordant alleles
  
  
  # START vs Mahajan 
    discordAll2<- CheckSnpAlleles(common,"A1_s_impd", "A2_s_impd", "RiskAllele", "OtherAllele", "chrpos") #  7446 SNPs with discordant alleles

# Get common SNPs (without discordant alleles )
    common2 <- common[-which(common$chrpos %in% c(discordAll1[[2]]$chrpos, discordAll2[[2]]$chrpos, discordAll1[[4]]$chrpos)),] #remove discordant alleles + 206 indels Indels 



# Keep SNPs with MAF >= 0.01 in one of the two studies 
  common3<- common2[which(common2$MAF_s_impd >= 0.01 | common2$MAF_bib_impd_hceGsa >= 0.01),]
  #Visual checks
  nrow(common3[which(common3$MAF_s_impd >= 0.01 & common3$MAF_bib_impd_hceGsa < 0.01),]) 
  nrow(common3[which(common3$MAF_s_impd < 0.01 & common3$MAF_bib_impd_hceGsa >= 0.01),]) 
  summary(common3[which(common3$MAF_s_impd >= 0.01 & common3$MAF_bib_impd_hceGsa < 0.01),"MAF_bib_impd_hceGsa"]) 
  summary(common3[which(common3$MAF_s_impd < 0.01 & common3$MAF_bib_impd_hceGsa >= 0.01),"MAF_s_impd"]) 

# Keep only Mahajan SNPs tested in South Asians

  common4 <- common3[!is.na(common3$Direction_SA),]

## Keep SNPs tested in >= 85% of Mahajan's N 

  maxn<-max (mah$N)
  common4$prc_sample_in_snp <- ((common4$N)/maxn)*100
  common5 <- common4[common4$N >= 85,] 
  nrow(common5)

# Save Output files 
  for (ps in length(grid2p):1){
    cat (ps, "\n")
    p <-grid2p[ps]
    if (p!=1) common5 <- common5[common5$P<= p, ]
    if (!file.exists(paste0(r5_p, "/Genotypes/", p)))  dir.create(paste0(r5_p, "/Genotypes/", p), recursive=T)
    if (!file.exists(paste0(r5_p, "/Phenotypes/START"))) dir.create(paste0(r5_p, "/Phenotypes/START"), recursive=T)
    if (!file.exists(paste0(r5_p, "/Phenotypes/BiB"))) dir.create(paste0(r5_p, "/Phenotypes/BiB"), recursive=T)
  # save the name of variants to be renamed
    write.table(common5[, c( "CHR", "POS", "POS",       "RSID_s_impd")], paste0(r5_p, "/Genotypes/", p, "/SNPs_in_START_BiB_Mah14_1KG_P_se_", p, "_ChrPos.txt"), col.names=F, row.names=F, quote=F) #START
    if(p==1){
      write.table(common5[, c("RSID_b_impd_hceGsa", "RSID_s_impd")],  paste0(r5_p, "/Genotypes/",p,"/BiB_SNPs_in_START_BiB_Mah14_1KG_P_se_", p, "_RSID_toRename.txt"), col.names=F, row.names=F, quote=F)#BiB
      write.table(common5[, c("RSID_kg", "RSID_s_impd")],  paste0(r5_p, "/Genotypes/",p,"/1KG_SNPs_in_START_BiB_Mah14_1KG_P_se_", p, "_RSID_toRename.txt"), col.names=F, row.names=F, quote=F)#BiB   
    }
  
  }
