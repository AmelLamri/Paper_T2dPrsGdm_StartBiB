#!/usr/bin/env bash
# Written by Amel Lamri
# Last modified on 16/11/2022

###############
# DESCRIPTION #
###############  

  # This script is associated to the manuscript titled "The genetic risk of 
  # gestational diabetes in South Asian women"  by lamri et al bublished at eLife
  # in 2022, which includes data from the START and BiB studies. More specifically,
  # this script run all the bash PLINK commands as well as the subscripts 
  # associated.
  # For more information, contact Amel Lamri (lamria@mcmaster.ca) and / or Sonia S Anand (anands@mcmaster.ca)

# Define input paths and variables 
  
  h_p=/home/lamria
  m_p=$h_p/srvstore/Merged_Studies/BiB_START_1KG_Mahajan2014/202108_Derive_T2Dprs
  g_p=$m_p/Genotypes
  extrctfnp=$g_p/1/SNPs_in_START_BiB_Mah14_1KG_P_se_1_ChrPos.txt
  omlfnp2=$g_p/1/START_BiB_1KG_SNPs_in_START_BiB_Mah14_1KG_P_se_1_mergeList.txt
  omlfnp3=$g_p/1/START_BiB_SNPs_in_START_BiB_Mah14_1KG_P_se_1_mergeList.txt
  ogfnp3=$g_p/1/START_BiB_1KG_SNPs_in_START_BiB_Mah14_1KG_P_se_1 
  ogfnp4=$g_p/1/START_BiB_SNPs_in_START_BiB_Mah14_1KG_P_se_1 

  s_p=$h_p/avr/Projects/byStudy/MixedStudies/START_BiB/202109_T2dPRSxRiskFactors_v3/Scripts/GDM_IADPSG/Public_Script_elife
  s_fne=202107_START_BiB_LDpred2_T2D_PRS_Main.bash
  Subs1_fne=202107_START_BiB_LDpred2_T2D_PRS_S1_PrepCommonSNPs.R
  Subs2_fne=202107_START_BiB_LDpred2_T2D_PRS_S2_Derive_PRS.R

   extrctfnp=$g_p/1/SNPs_in_START_BiB_Mah14_1KG_P_se_1_ChrPos.txt
   omlfnp2=$g_p/1/START_BiB_1KG_SNPs_in_START_BiB_Mah14_1KG_P_se_1_mergeList.txt
   omlfnp3=$g_p/1/START_BiB_SNPs_in_START_BiB_Mah14_1KG_P_se_1_mergeList.txt
   ogfnp3=$g_p/1/START_BiB_1KG_SNPs_in_START_BiB_Mah14_1KG_P_se_1 
   ogfnp4=$g_p/1/START_BiB_SNPs_in_START_BiB_Mah14_1KG_P_se_1 
 
  # Related Scripts file names 
   Subs1_fne=START_BiB_LDpred2_T2D_PRS_S1_PrepLDpredFiles.R
   Subs2_fne=START_BiB_LDpred2_T2D_PRS_S2_Derive_PRS.R
   Subs3_fne=START_BiB_LDpred2_T2D_PRS_S3_Merge_PRS_Phenos.R
   Subs4_fne=START_BiB_LDpred2_T2D_PRS_S4_RunAssociation.R

  
# Get the list of common SNPs between START BIB 1KG, Mahajan2014
  $s_p/$Subs1_fne

# Subset SNPs of interest and merge files  
 
  if [ ! -f $ogfnp3 ]; then 
    
    for std in START BiB 1KG; do 
      
      ogfnp2=$g_p/1/$std/${std}_SNPs_in_START_BiB_Mah14_1KG_P_se_1
    
      omlfnp1=$g_p/1/$std/${std}_SNPs_in_START_BiB_Mah14_1KG_P_se_1_MergeList.txt
        
      for chr in {1..22};do
    
        if [ $std = "START" ] ;then
          igfnp1=/START/GWAS_HCE_ICE/Filtered/V1/InfoScore0.7/START_GW_HCE_ICE_hg19_updated_phased_imputed_chr${chr}    _filtered_info0.7_Miss0.05_CallThres0_UpdatedAllelesAndIDsAndParents
        elif [ $std = "BiB" ] ;then
          igfnp1=/Born_in_Bradford/Genotypes/202101_1KG_Imputation/BiB_SouthAsian_Mothers_Imputed/CoreExome_GSA/Chr${chr}/    SNP_Subsets/Info0.7/BiB_SouthAsian_Mothers_CoreExomeGSA_Imputed_1KG_Chr${chr}_Info0.7
        elif [ $std = "1KG" ] ;then
          igfnp1=/1000Genomes/20130502_Release/Genotypes/PLINK/Sample_Subsets/South_Asians/1000G_SA_chr${chr}
        fi
        
        ogfnpref1=$g_p/1/$std/${std}_SNPs_in_START_BiB_Mah14_1KG_P_se_1_chr
        ogfnp1=$ogfnpref1${chr}
        
        # Subset genotypes of interest
        plink \
        --bfile $igfnp1 \
        --extract range $extrctfnp \
        --make-bed \
        --out $ogfnp1
        
        # Merge files by chromosome
        if [ $chr -eq 1 ];then
          echo $ogfnp1 > $omlfnp1    
        elif [ $chr -eq 22 ];then
          echo $ogfnp1 >> $omlfnp1
          plink \
          --merge-list $omlfnp1 \
          --make-bed \
          --out $ogfnp2    
          
          # Delete files by chromosome
          if [ -f $ogfnp2.bed ] ; then
            rm ${ogfnpref1}{1..22}.{bed,bim,fam,nosex}
          fi      
        else
          echo $ogfnp1 > $omlfnp1    
        fi 
    
      done
    
    
      # Merge START BiB and 1KG 
      if [ $std = START ];then
        echo $ogfnp2 > $omlfnp2    
      elif [ $std = BiB ];then
        echo $ogfnp2 >> $omlfnp2
      elif [ $std = 1KG ];then
        echo $ogfnp2 >> $omlfnp2
        plink \
        --merge-list $omlfnp2 \
        --make-bed \
        --out $ogfnp3
       
    # Merge START & BiB  
      if [ $std = START ];then
        echo $ogfnp2 > $omlfnp3    
      elif [ $std = BiB ];then
        echo $ogfnp2 >> $omlfnp3
        plink \
        --merge-list $omlfnp3 \
        --make-bed \
        --out $ogfnp4
        
        # Delete files by study
        if [ -f $ogfnp3.bed ] ; then
             rm $g_p/1/{START,BiB,1KG}/{START,BiB,1KG}_SNPs_in_START_BiB_Mah14_1KG_P_se_1.{bed,bim,fam,nosex}
        fi
      
      fi
      
    done
  fi

# Impute missing genotypes  

  plink \
    --bfile ${ogfnp3} \
    --fill-missing-a2 \
    --make-bed \
    --out ${ogfnp3}_noMiss
  
  plink \
    --bfile ${ogfnp4} \
    --fill-missing-a2 \
    --make-bed \
    --out ${ogfnp4}_noMiss

# Check number of variants is correct 

  nl1=$( wc -l < $ogfnp3.bim )
  nl2=$( wc -l <  )
  if [ $nl1 -eq $nl2 ];then echo NUMBER OF VARIANTS OK \n ;fi

# Derive PRS 
  $s_p/$Subs2_fne

# Merge PRS with Phenos
  $s_p/$Subs3_fne

# Run Association tests 

  $s_p/$Subs4_fne
