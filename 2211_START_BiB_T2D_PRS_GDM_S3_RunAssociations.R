#!/usr/env/bin/Rscript
# Written by Amel Lamri
# Last modified on 16/11/2022

###############
# DESCRIPTION #
###############  

  # This script is associated to the manuscript titled "The genetic risk of 
  # gestational diabetes in South Asian women"  by lamri et al bublished at eLife
  # in 2022, which includes data from the START and BiB studies. More specifically,
  # this script run all the association tests described in the manuscript. The 
  # outputs include Tables 1 to 4, figures 1 and 2, and supplementary tables a-f
  # For more information, contact Amel Lamri (lamria@mcmaster.ca) and / or Sonia S Anand (anands@mcmaster.ca)

# Load relevant libraries 
  library("tidyverse")
  library("jtools")
  library("haven")
  library("ggplot2")
  library("readxl")
  library("sas7bdat")
  library("sjPlot") # interaction plot
  library("metafor") #meta-analaysis 
  library("AF") # Attributable fraction 
  library("PGPmisc") #library("PGPmisc") # my own functions
  library("interactions") #for interact_plot
  library("InteractionPoweR")

# Define paths and files names 

  # General
    lmd <- "220922"
    avr_p <- "/home/lamria/avr"
    prjcts_p <- paste0(avr_p, "/Projects/byStudy")
    pmsp_p <- paste0(prjcts_p, "/MixedStudies")
    sbib_p <- paste0(pmsp_p, "/START_BiB/")
    m_p    <- paste0(sbib_p, "/202109_T2dPRSxRiskFactors_v3")
    s_p    <- paste0(m_p, "/Scripts")
    o_p    <- paste0(m_p, "/Outputs")
    o2_p    <- paste0(o_p, "/v2206")
      
  # Other input files & paths  
    prs_p <- paste0(pmsp_p, "/START_BiB/202108_Derive_LDpred2_T2D_PRS/Results")
    bib_p <- paste0(prjcts_p, "/Born_in_Bradford")
    bibfo_p <- paste0(bib_p, "/202109_Mothers_FoodIntake_v2/Outputs")
    comsnpinfp_p <- paste0(pmsp_p, "/START_BiB_1KG_Mahajan2014/202107_ListCommonSNPs/Results") #results folder
  
  #Input files 
  
  
    sifnpe1 <- paste0(o2_p, "/START_BIB_T2D_GDMiadpsg_PRS_START_Data_20220619.txt")# START  phenotypes for analysis
    bifnpe1 <-  paste0(o2_p, "/START_BIB_T2D_GDMiadpsg_PRS_BiB_Data_20220619.txt") # BiB phenotypes for analysis 
    sbifnpe1 <- paste0(sbib_p, "/VariableCorrespondence/Outputs/START_BiB_Variablescorrespondence_20220619.tbl") # START - BiB Variable correspondance table 
    
    suplt1fne<- paste0("TableS1_GDMiadpsg_SNPsInPRSInfo_", lmd, ".csv")
    suplt2fne<- paste0("TableS2_GDMiadpsg_UnivarAssocPRS_", lmd, ".csv")
    suplt3fne<- paste0("TableS3_GDMiadpsg_MultivarAssocPRSter_", lmd, ".csv")
    
    bprs_fnpe <- paste0(prs_p, "/BiB_T2D_PRS_Mah14_LDpred2BestGrid_20220609.csv")
    sprs_fnpe <- paste0(prs_p, "/START_T2D_PRS_Mah14_LDpred2BestGrid_20220609.csv")


# Open Datasets 

  # Phenotypes 
    # START 
      st0 <- read.delim(sifnpe1, string=F)
      st0$Study <- "START"
  
    # BiB
      bib0 <- read.delim(bifnpe1, string=F)
      bib0$Study <- "BiB"
  
  # PRS 
    # START 
      stprs <- read.csv( sprs_fnpe, string=F)
    #BiB
      bibprs <- read.csv( bprs_fnpe, string=F)


  # Variable correspondance table 
    
    cvt<- read.delim(sbifnpe1,  head=T, string=F)
    #The cvt table includes the correspondance of variables names between the 
    # two studies (START and BiB ) for example : 
    #________________________________________
    #|VarGleShrt | VarNameSTART | VarNameBiB |
    #|___________|______________|____________|
    #| sex       |momSex        | sexmb      |
    #| bmi       |momBmi        | bmimb      |
    #|___________|______________|____________|
    # in this illustrative example, the name of the sex vairable in START is 
    # momSex vs. sexmb in BiB  



# Merge PRS and phenotype data 
  # Visually check common columns
  colnames(stprs)[colnames(stprs)%in%colnames(st10)]
  colnames(bibprs)[colnames(bibprs)%in%colnames(bib0)]
  # Merge
  st1 <- merge(stprs, st0, all.y=T, all.y=F) 
  bib1 <- merge(bibprs, bib0, all.y=T, all.y=F) 

# Set levels of factor variables 
  # This is done so the different levels will be compared/showed in the right order
  st1$"BornCaSaOth" <- factor(st1$"BornCaSaOth", levels=c("South Asia", "Canada", "Other"))
  bib1$"BornUkSaOth_AL" <- factor(bib1$"BornUkSaOth_AL", levels=c("South Asia", "United Kingdom", "Other"))

# Store data in a list
  # This is done so I can apply the same functions to both START and BiB later on 
  
  sbl <- list()
  sbl[[1]] <- st1
  sbl[[2]] <- bib1
  names(sbl) <- c("START", "BiB")

# Subset variables of interest to the analysis
  vis <- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg","GDM", "AUCg", "AUClwz", "Ogtt2hG", "FPG", "T2dMah14Prs","BMI_v2", paste0("PC", 1:5), "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "lowDQS", "bmi_s23vle23_v2", "parity_2cat","Education_2c", "age_s29vle2931", "age_s29vle32","T2dMah14PrsT1vT2", "T2dMah14PrsT1vT3", "Ogtt1hG", "T2dMah14PrsT1vT23", "T2dMah14PrsTer", "Vegetarian", "YrInCtry", "CountryOrigin", "BirthCountry3c","Height", "Weight","Age", "parity_3cat", "ihunhDQS", "PregnancyNumber","CountryOrigin3c", "T2dMah14PrsZ", "T2dMah14PrsT12vT3", "T2dMah14PrsZtop10"), ]

# Table 1: Characteristics by GDM status 
  # Caracteristics table with the N (%), and Mean (SD) by GDM case/control status in START and BiB
  # Generate data for continuous variables 
    
    
    # Variables as called in cvt column "VarGleShrt"
     cvit1<- c("Age", "Height", "Weight", "BMI_v2", "T2dMah14PrsZ", "YrInCtry", "FPG", "Ogtt2hG", "Ogtt1hG", "AUCg")
    
      cvt[cvt$VarGleShrt %in% cvit1, c("VarGleShrt", "VarClass_Gle", "TableLegend", "Round")]
      # the charTblByCCtrl function produces mean(SD) for continuous vars
    
      t1s<-list()
      for(std in 1:2){
        t1s[[std]]<- PGPmisc:::charTblByCCtrl(OutVars=cvt[cvt$VarGleShrt %in% cvit1,std+1], CCvar=cvt[cvt$VarGleShrt=="GDM_iadpsg",std+1], data=sbl[[std]], legend=cvt      [cvt$VarGleShrt%in% cvit1, "TableLegend"] , VarTypes=cvt[cvt$VarGleShrt%in% cvit1, "VarClass_Gle"], rounding=5)
        }
    
    # round numbers 
      for (i in 1:length(cvit1)){
        vari <- cvt[cvt$VarGleShrt ==cvit1[i], "TableLegend"]
        roundi <- cvt[which(cvt$TableLegend == vari), "Round"]
        for (j in 1:2){
          cat(i, j, "rounding", vari, roundi, "\n")
          t1s[[j]][which(t1s[[j]]$Variable==vari),"Ctrl_MeanOrN"] <- round(as.numeric(t1s[[j]][which(t1s[[j]]$Variable==vari),"Ctrl_MeanOrN"]), roundi)
          t1s[[j]][which(t1s[[j]]$Variable==vari),"Ctrl_SDor%"  ] <- round(as.numeric(t1s[[j]][which(t1s[[j]]$Variable==vari),"Ctrl_SDor%"  ]), roundi)
          t1s[[j]][which(t1s[[j]]$Variable==vari),"Case_MeanOrN"] <- round(as.numeric(t1s[[j]][which(t1s[[j]]$Variable==vari),"Case_MeanOrN"]), roundi)
          t1s[[j]][which(t1s[[j]]$Variable==vari),"Case_SDor%"  ] <- round(as.numeric(t1s[[j]][which(t1s[[j]]$Variable==vari),"Case_SDor%"  ]), roundi)
        
          }
        }
    
    colnames(t1s[[1]]) <-c( "Variable",  c(paste0("START_", colnames(t1s[[1]])[2:6])))
    colnames(t1s[[2]]) <-c( "Variable",  c(paste0("BiB_", colnames(t1s[[2]])[2:6])))
    table(t1s[[1]]$Variable %in% t1s[[1]]$Variable)
    t1 <-merge( t1s[[1]],  t1s[[2]], all=T)
    
 
  # Add the first row ( N and % of GDM)   
    #This one is done separately because it does not rely on a variable, it's based on the number of rows
    t1r1 <- t1[0,]
    t1r1[1,1] <- "N (%)"
    t1r1[1,2] <- as.numeric(table(sbl[[1]]$anygdm, useNA="ifany"))[1]
    t1r1[1,3] <- as.numeric(round(prop.table(table(sbl[[1]]$anygdm, useNA="ifany"))*100,1))[1]
    t1r1[1,4] <- as.numeric(table(sbl[[1]]$anygdm, useNA="ifany"))[2]
    t1r1[1,5] <- as.numeric(round(prop.table(table(sbl[[1]]$anygdm, useNA="ifany"))*100,1))[2]
    
    t1r1[1,7] <- as.numeric(table(sbl[[2]]$drvgesdiab01, useNA="ifany"))[1]
    t1r1[1,8] <- as.numeric(round(prop.table(table(sbl[[2]]$drvgesdiab01, useNA="ifany"))*100,1))[1]
    t1r1[1,9] <- as.numeric(table(sbl[[2]]$drvgesdiab01, useNA="ifany"))[2]
    t1r1[1,10] <- as.numeric(round(prop.table(table(sbl[[2]]$drvgesdiab01, useNA="ifany"))*100,1))[2]
  
  

  
  # Generate data for categorical variables 
    ctvit1<- c("T2dMah14PrsTer", "ParHistDiabYN", "parity_3cat", "BornInSaTF", "Education_2c", "Vegetarian",   "CountryOrigin3c","Vegetarian")
    
    t1cat <- t1[0,] 
    #1:length(ctvit1)
    for( vari in 1:length(ctvit1)){
    vit <- ctvit1[vari] 
      for ( stdi in 1:2){
        tmpd<- sbl[[stdi]]
        gdmvni <- cvt[cvt$VarGleShrt=="GDM_iadpsg",stdi+1]
        vit2   <- cvt[cvt$VarGleShrt==vit, stdi+1]  
        tb <- table(tmpd[tmpd[,gdmvni]==0, vit2] )
        tb2 <- table(tmpd[tmpd[,gdmvni]==1, vit2] )
        
        if(stdi==1){
          colz<- 2:6
        }else if(stdi==2){
          colz<- 7:11
        }
        
        if(stdi==1) {
          tbl3orM <- as.data.frame(cbind(tb,round(prop.table(tb)*100,1), tb2,round(prop.table(tb2)*100,1)))
          colnames(tbl3orM) <- colnames(t1cat)[colz][1:4]
          tbl3orM[,colnames(t1cat)[colz][5]] <- PGPmisc:::format1Pval(as.numeric(chisq.test(cbind(tb,tb2))$"p.value"))
          tbl3orM[,colnames(t1)[1]] <- cvt[cvt$VarGleShrt==vit,"TableLegend"]
    	  tbl3orM[,"Var"] <- cvt[cvt$VarGleShrt==vit,1]
          tbl3orM[,"Subgroup"] <- names(tb)
        }else if(stdi==2) {
          tbl3orM2 <- as.data.frame(cbind(tb,round(prop.table(tb)*100,1), tb2,round(prop.table(tb2)*100,1)))
          colnames(tbl3orM2) <- colnames(t1cat)[colz][1:4]
          tbl3orM2[,colnames(t1cat)[colz][5]] <- PGPmisc:::format1Pval(as.numeric(chisq.test(cbind(tb,tb2))$"p.value"))
          tbl3orM2[,"Var"] <- cvt[cvt$VarGleShrt==vit,1]
    	  tbl3orM2[,colnames(t1)[1]] <- cvt[cvt$VarGleShrt==vit,"TableLegend"]
          tbl3orM2[,"Subgroup"] <- names(tb)
          colnames(tbl3orM2) [ colnames(tbl3orM2) %in% colnames(tbl3orM2) ]
          tbl3orM3 <- merge(tbl3orM2, tbl3orM, all=T)
        }
      }
      t1cat <- merge(t1cat,tbl3orM3, all=T)
    }
  

  
  # Run n(%) of lowDQS (low diet quality variable)
    # lowDQS is ran separately because the variable is available in sTART only (not BiB)
    
    dqlt <- t1[0,1:6] 
    vit<- "lowDQS"
    tmpd<- sbl[[1]]
    stdi <- 1
    gdmvni <- cvt[cvt$VarGleShrt=="GDM_iadpsg",stdi+1]
    vit2   <- cvt[cvt$VarGleShrt==vit, stdi+1]
    dqlt[1,1] <- cvt[cvt$VarGleShrt==vit, "TableLegend"]
    dqlt[1,"Var"]  <- vit2
    
    ## GDM Controls
    tmpxd1 <- table(tmpd[tmpd[,gdmvni] ==0, vit2])
    tmpxd2 <- table(tmpd[tmpd[,gdmvni] ==1, vit2])
    dqlt[1,2] <- tmpxd1[2]
    dqlt[1,3] <- round(prop.table(tmpxd1)*100,1)[2]
    dqlt[1,4] <- tmpxd2[2]
    dqlt[1,5] <- round(prop.table(tmpxd2)*100,1)[2]
    dqlt[1,6] <- PGPmisc:::format1Pval(as.numeric(chisq.test(cbind(tmpxd1,tmpxd2))$"p.value"))
  
  

  
  # sort rows in the desired order  

    #This is done so I can sort the rows of the table the way I want to  
    t1r1[, "order"] <- 1
    t1[t1$Variable ==  "Age, years"								, "order"] <- 2
    t1[t1$Variable ==  "Height, cm"       								, "order"] <- 3
    t1[t1$Variable ==  "Weight, kg"       								, "order"] <- 4
    t1[t1$Variable ==  "BMI, kg/m2"       								, "order"] <- 5
    t1cat[t1cat$"Var" == "Parity 3 categories"		& t1cat$Subgroup==0         , "order"] <- 6
    t1cat[t1cat$Var == "Parity 3 categories"		& t1cat$Subgroup==1         , "order"] <- 7
    t1cat[t1cat$Var == "Parity 3 categories"		& t1cat$Subgroup==2         , "order"] <- 8
    t1cat[t1cat$Var == "Education level 2 categories (Pre vs. Post-secondary)"		& t1cat$Subgroup==1         ,   "order"] <- 9
    t1cat[t1cat$Var == "Country of origin (3 categories)"	& t1cat$Subgroup=="India"   , "order"] <- 10
    t1cat[t1cat$Var == "Country of origin (3 categories)"	& t1cat$Subgroup=="Pakistan", "order"] <- 11
    t1cat[t1cat$Var == "Country of origin (3 categories)"	& t1cat$Subgroup=="Other"   , "order"] <- 12
    t1cat[t1cat$Var == "Born in South Asia [T/F]"			& t1cat$Subgroup==TRUE      , "order"] <- 13
    t1[t1$Variable ==  "Years in recruitment country (Canada/UK)"     , "order"] <- 14
    t1cat[t1cat$Var == "Parental History of Diabetes"	& t1cat$Subgroup==1             , "order"] <- 15
    t1cat[t1cat$Var == "Vegetarian status [T/F]" 		& t1cat$Subgroup==TRUE      , "order"] <- 16
    t1cat[t1cat$Var == "Interheart Unhealthy foods score"			& t1cat$Subgroup==0         , "order"] <- 17
    t1cat[t1cat$Var == "Interheart Unhealthy foods score"			& t1cat$Subgroup==1         , "order"] <- 18
    t1cat[t1cat$Var == "Interheart Unhealthy foods score"			& t1cat$Subgroup==2         , "order"] <- 19
    dqlt[, "order"] <- 20
    t1[t1$Variable == "Polygenic risk score (z-scores)"     , "order"] <- 21
    t1cat[t1cat$Var == "T2D PRS Tertiles (LDpred2)"	& t1cat$Subgroup==1         , "order"] <- 22
    t1cat[t1cat$Var == "T2D PRS Tertiles (LDpred2)"	& t1cat$Subgroup==2         , "order"] <- 23
    t1cat[t1cat$Var == "T2D PRS Tertiles (LDpred2)"	& t1cat$Subgroup==3         , "order"] <- 24
    t1[t1$Variable == "Fasting plasma glucose, mmol/L"		, "order"] <- 25
    t1[t1$Variable == "1h post-load glucose, mmol/L"		, "order"] <- 26
    t1[t1$Variable == "2h post-load glucose, mmol/L"		, "order"] <- 27
    t1[t1$Variable == "Area under glucose curve (mmolh)"		, "order"] <- 28
    t1_b <- merge(merge(merge(t1cat[!is.na(t1cat$order),], t1, all=T) ,dqlt, all=T),t1r1, all=T)


  # Format and clean Output table 
   
    # correct the wording in the first column of the table 
    t1_b[t1_b$Variable=="Parity (0 vs 1 vs 2+)", "Variable"] <- "Parity, n (%)"
    t1_b[t1_b$Variable=="Born in South Asia", "Variable"] <- "Born in South Asia, n (%)"
    t1_b[t1_b$Variable=="Polygenic risk Score (continuous, z-scores)", "Variable"] <- "Polygenic risk score (z-scores)"
    t1_b[t1_b$Variable=="Polygenic risk score (categorical), n (%)", "Variable"] <- "Polygenic risk score (Tertiles), n(%)"
    t1_b[t1_b$Variable=="Education level (Post-secondary)", "Variable"] <- "Post-secondary education, n(%)"
    t1_b[t1_b$Variable=="Country of origin/ancestry", "Variable"] <- "Country of origin/ancestry, n(%)"
    
    catwsub <- c("Parity, n (%)", "Country of origin/ancestry, n(%)", "Polygenic risk score (Tertiles), n(%)")
    
    vwprc <-c(t1_b[!is.na(t1_b$Subgroup),"Variable"], "n (%)")
    
    # Order rows
    t1_b <- t1_b[order(t1_b$order), ]
  
    # Formatting : paste the values with () 
    t1_b[t1_b$Variable %in% vwprc, "No GDM_st" ] <- paste0(t1_b[t1_b$Variable %in% vwprc, "START_Ctrl_MeanOrN"], " (",t1_b  [t1_b$Variable %in% vwprc, "START_Ctrl_SDor%"], "%)")
    t1_b[t1_b$Variable %in% vwprc, "GDM_st" 	] <- paste0(t1_b[t1_b$Variable %in% vwprc, "START_Case_MeanOrN"], " (",t1_b  [t1_b$Variable %in% vwprc, "START_Case_SDor%"], "%)")
    t1_b[t1_b$Variable %in% vwprc, "No GDM_bib"] <- paste0(t1_b[t1_b$Variable %in% vwprc, "BiB_Ctrl_MeanOrN"  ], " (",t1_b  [t1_b$Variable %in% vwprc, "BiB_Ctrl_SDor%"  ], "%)")
    t1_b[t1_b$Variable %in% vwprc, "GDM_bib" 	] <- paste0(t1_b[t1_b$Variable %in% vwprc, "BiB_Case_MeanOrN"  ], " (",t1_b  [t1_b$Variable %in% vwprc, "BiB_Case_SDor%"  ], "%)")
    
    t1_b[ !t1_b$Variable %in% vwprc, "No GDM_st" ] <- paste0(t1_b[ !t1_b$Variable %in% vwprc, "START_Ctrl_MeanOrN"], " (",  t1_b[ !t1_b$Variable %in% vwprc, "START_Ctrl_SDor%"], ")")
    t1_b[ !t1_b$Variable %in% vwprc, "GDM_st" 	] <- paste0(t1_b[ !t1_b$Variable %in% vwprc, "START_Case_MeanOrN"], " (",  t1_b[ !t1_b$Variable %in% vwprc, "START_Case_SDor%"], ")")
    t1_b[ !t1_b$Variable %in% vwprc, "No GDM_bib"] <- paste0(t1_b[ !t1_b$Variable %in% vwprc, "BiB_Ctrl_MeanOrN"  ], " (",  t1_b[ !t1_b$Variable %in% vwprc, "BiB_Ctrl_SDor%"  ], ")")
    t1_b[ !t1_b$Variable %in% vwprc, "GDM_bib" 	] <- paste0(t1_b[ !t1_b$Variable %in% vwprc, "BiB_Case_MeanOrN"  ], " (",  t1_b[ !t1_b$Variable %in% vwprc, "BiB_Case_SDor%"  ], ")")
    
    t1_b[t1_b$Variable %in% catwsub, "Variable" ] <- paste0(t1_b[t1_b$Variable %in% catwsub, "Variable" ], " - ", t1_b[t1_b$Variable %in% catwsub, "Subgroup" ])
  
    t1_b[t1_b$Variable=="Parity, n (%) - 2", "Variable"]  <- "Parity, n (%) - 2 or more"
    t1_b[t1_b$Variable=="Parity, n (%) - 2", "Variable"]  <- "Parity, n (%) - 2 or more"
  
  # Save output table 
    write.table(t1_b[,c("Variable","No GDM_st","GDM_st", "START_PVal","No GDM_bib","GDM_bib", "BiB_PVal")] , "/Table1.tbl",sep="\t", row.names=F, quote=F)
 

# Table 2: Associations : GDM/AUC ~ PRS + RFs wAdj 

  # Prep tables for analysis
    tbl2_depvs_d <- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg","AUClwz", "Ogtt2hG", "FPG"), ]
    tbl2_indepvs_d <- cvt[cvt$VarGleShrt %in% c("T2dMah14PrsZ"),]
    
    tbl2_adjs_l <- list()
    tbl2_adjs_l[[1]] <- cvt[cvt$VarGleShrt %in% c("Age","BMI_v2", paste0("PC", 1:5), "Education2_n","BornInSaTF",   "Parity",  "ParHistDiabYN", "lowDQS"), ]
    names(tbl2_adjs_l)<-paste0("model", 1:length(tbl2_adjs_l))
    tbl2pref <- "Table2_GDMiadpsg_AssocwAdj"

  # Run association tests 
    #This section runs association tests on all the variables of interest in each study separately 
    tbl2o1l <- PGPmisc:::AssocInterMultiStudies(data=sbl, depVs_d=tbl2_depvs_d, indepVs_d=tbl2_indepvs_d, adjVs_l=tbl2_adjs_l,depVs_GleName_Coli=1, indepVs_GleName_Coli=1,adjVs_GleName_Coli=1, depVs_ColNames_Coli=2:3, indepVs_ColNames_Coli=2:3, depVs_ModelTypes_Coli=7, adjVs_ColNames_Coli=2:3,  saveShortOutput=F, saveLongOutput=T, OutputPath=o2_p, OutFileNamePrefx=tbl2pref, OutFileNameSufx=paste0("_", lmd)) [[2]]

  # Horizontally merge START and BiB association test results
  
    tbl2o1l_st <- tbl2o1l[tbl2o1l$Study== "START", ]
    tbl2o1l_bib <- tbl2o1l[tbl2o1l$Study== "BiB", ]
    tbl2colsToKeep<- c("DepVarGle", "IndepVarIGle","IndepVarSub_gleVars", "Model", "Ndiai")
    
    colnames(tbl2o1l_st)[which(!colnames(tbl2o1l_st) %in% tbl2colsToKeep)] <- paste0(colnames(tbl2o1l_st)[which(!colnames    (tbl2o1l_st) %in% tbl2colsToKeep)], "_st")
    colnames(tbl2o1l_bib)[which(!colnames(tbl2o1l_bib) %in% tbl2colsToKeep)] <- paste0(colnames(tbl2o1l_bib)[which  (!colnames  (tbl2o1l_bib) %in% tbl2colsToKeep)], "_bib")
  
    tbl2o2_stbib <- merge(tbl2o1l_st, tbl2o1l_bib, tbl2colsToKeep, suf=c("_st","_bib"), all=T)
    #visually check the merging 
    checkmt<-unique(tbl2o2_stbib[tbl2o2_stbib$IndepVarSub_gleVars != "(Intercept)", c("IndepVarSub_gleVars",     "IndepVarSub_bib", "IndepVarSub_st")])
    checkmt
    # Visually, merging is ok ! 

  # Remove rows of (Intercept) and PCA, and  lowDQS rows in tbl2o2_stbib before meta-A 
  
    # lowDQS is removed because it is only adjusted for in START
    # PCAs are removed because they are different between studies 
    tbl2o2_stbib_2 <- tbl2o2_stbib[-which(tbl2o2_stbib$IndepVarSub_gleVars %in% c("(Intercept)", paste0("PC", 1:5),   "lowDQS")),]


  # Run meta-analysis FE
    #used fixed effects models

    tbl2o2_l<-list()
    metaAMeth <- "FE"
    for( i in 1:nrow(tbl2o2_stbib_2)){
        tbl2o2_l[[i]] <- list()
    	if(!is.na(tbl2o2_stbib_2[i, "Beta_st"]) & !is.na(tbl2o2_stbib_2[i, "Beta_bib"])){
    
    	esci2 <- escalc(yi=as.numeric(tbl2o2_stbib_2[i,paste0("Beta", c("_st", "_bib"))]), sei=as.numeric(tbl2o2_stbib_2[i,    paste0("SE", c("_st", "_bib"))]), measure="SMD")
    
    	mai2 <- rma(yi, vi, method=metaAMeth, data=esci2)
    
    	tbl2o2_l[[i]] <- mai2
    
    	tbl2o2_stbib_2[i,"MetaA_Beta"]  <- mai2$"beta"
        tbl2o2_stbib_2[i,"MetaA_SE"]    <- mai2$"se"
        tbl2o2_stbib_2[i,"MetaA_Beta95%CiLB"] <- mai2$"ci.lb"
        tbl2o2_stbib_2[i,"MetaA_Beta95%CiUB"] <- mai2$"ci.ub"
        tbl2o2_stbib_2[i,"MetaA_P"]     <- mai2$"pval"
        tbl2o2_stbib_2[i,"MetaA_I2"]    <- mai2$"I2"
        tbl2o2_stbib_2[i,"MetaA_QE"]    <- mai2$"QE"
        tbl2o2_stbib_2[i,"MetaA_QEp"]   <- mai2$"QEp"
        tbl2o2_stbib_2[i,"MetaA_Model"] <- metaAMeth
        	
    	if(!is.na(tbl2o2_stbib_2[i, "OR_st"])){
          tbl2o2_stbib_2[i,"MetaA_OR"]         <- exp(mai2$"beta")
    	  tbl2o2_stbib_2[i,"MetaA_95%CiLB"]         <- exp(mai2$"ci.lb")
    	  tbl2o2_stbib_2[i,"MetaA_95%CiUB"]         <- exp(mai2$"ci.ub")
        }
      }
    }


  # Add back dietq_lo to the START model infos 
    
    dietq_lo was only adjuted for in START, and hence association results between dietq_lo was not included in the     eta-analysis. I want to add it back to the  meta-A table 
    
    # prep the dietqlo data before I add it back 
    t_dql <- tbl2o1l_st[tbl2o1l_st$IndepVarSub_st ==  "dietq_lo",]
    
    bl2o2_stbib_3 <- merge(tbl2o2_stbib_2,st_dql, all=T)
  

  # Format output table

    # use Beta [95%CIL - 95%CIU] format
    # round P values
    # order the rows inthe  table
    # replace NA[NA-NA] by NA
    # ... etc 
    
    tbl2o2_stbib2 <- tbl2o2_stbib_3[,tbl2colsToKeep]
    tbl2o2_stbib2$"Beta (SE)_st" <- paste0(round(tbl2o2_stbib_3$"Beta_st",3), " (", round(tbl2o2_stbib_3$"SE_st",3), ")")
    
    tbl2o2_stbib2$"Beta [CI]_st" <- paste0(round(tbl2o2_stbib_3$"Beta_st",3), " [", round(tbl2o2_stbib_3$"Beta95%CiLB_st",3),    " - ",  round(tbl2o2_stbib_3$"Beta95%CiUB_st",3), "]")
    
    tbl2o2_stbib2$"OR [95%CI]_st" <- paste0(round(tbl2o2_stbib_3$"OR_st",2), " [", round(tbl2o2_stbib_3$"95%CiLB_st" ,2), "     - ", round(tbl2o2_stbib_3$"95%CiUB_st",2),"]")
    tbl2o2_stbib2$"P_st" <- PGPmisc:::formatPvals(tbl2o2_stbib_3$"P_st")
    
    
    tbl2o2_stbib2$"Beta (SE)_bib" <- paste0(round(tbl2o2_stbib_3$"Beta_bib",3), " (", round(tbl2o2_stbib_3$"SE_bib",3), ")")
    
    tbl2o2_stbib2$"Beta [CI]_bib" <- paste0(round(tbl2o2_stbib_3$"Beta_bib",3), " [", round(tbl2o2_stbib_3$"Beta95%CiLB_bib",    3)," - ",  round(tbl2o2_stbib_3$"Beta95%CiUB_bib",3), "]")
    
    tbl2o2_stbib2$"OR [95%CI]_bib" <- paste0(round(tbl2o2_stbib_3$"OR_bib",2), " [", round(tbl2o2_stbib_3$"95%CiLB_bib" ,2),     " - ", round(tbl2o2_stbib_3$"95%CiUB_bib",2),"]")
    #some rows have NA P values for BiB because dietq_lo is not adjuted for in this study   
    w<-which(is.na(tbl2o2_stbib_3$"P_bib"))
    tbl2o2_stbib2[-w, "P_bib"] <- PGPmisc:::formatPvals(tbl2o2_stbib_3[-w,"P_bib"])
    
    
    tbl2o2_stbib2[-w,"Beta (SE)_MetaA"] <- paste0(round(tbl2o2_stbib_3[-w,"MetaA_Beta"],3), " (", round(tbl2o2_stbib_3[-w,    "MetaA_SE"],3), ")")
    
    tbl2o2_stbib2[-w,"Beta [CI]_MetaA"] <- paste0(round(tbl2o2_stbib_3[-w,"MetaA_Beta"],3), " [", round(tbl2o2_stbib_3[-w,    "MetaA_Beta95%CiLB"],3), " - ",  round(tbl2o2_stbib_3[-w,"MetaA_Beta95%CiUB"],3), "]")
    
    
    
    tbl2o2_stbib2[-w,"OR [95%CI]_MetaA"] <- paste0(round(tbl2o2_stbib_3[-w,"MetaA_OR"],2), " [", round(tbl2o2_stbib_3[-w,    "MetaA_95%CiLB"],2), " - ", round(tbl2o2_stbib_3[-w,"MetaA_95%CiUB"],2),"]")
    
    tbl2o2_stbib2[-w,"MetaA_P"] <- PGPmisc:::formatPvals(tbl2o2_stbib_3[-w,"MetaA_P"])
    tbl2o2_stbib2[-w,"I2_MetaA"] <- round(tbl2o2_stbib_3[-w,"MetaA_I2"],0)
    tbl2o2_stbib2[-w,"QEp_MetaA"] <- PGPmisc:::formatPvals(tbl2o2_stbib_3[-w,"MetaA_QEp"])
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"OR [95%CI]_MetaA"=="NA [NA - NA]"),"OR [95%CI]_MetaA"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"OR [95%CI]_st"=="NA [NA - NA]"),"OR [95%CI]_st"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"OR [95%CI]_bib"=="NA [NA - NA]"),"OR [95%CI]_bib"]<-NA
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta (SE)_bib"=="NA (NA)"),"Beta (SE)_bib"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta (SE)_st"=="NA (NA)"),"Beta (SE)_st"]<-NA
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta [CI]_bib"=="NA [NA - NA]"),"Beta [CI]_bib"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta [CI]_st"=="NA [NA - NA]"),"Beta [CI]_st"]<-NA
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta [CI]_MetaA"=="NA [NA - NA]"),"Beta [CI]_MetaA"]<-NA
    
    
    int2<- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg", "AUClwz", "Ogtt2hG", "FPG"), c(1:3, which(colnames(cvt) ==     "VarGleShrt")) ]
    
    colnames(int2)[1]<- "DepVarGle"
    int2[int2$VarGleShrt=="GDM_iadpsg",           "AssocTableDepVDescCol"] <- "GDM (iadpsg)"
    int2[int2$VarGleShrt=="GDM",           "AssocTableDepVDescCol"] <- "GDM (SA specific)"
    int2[int2$VarGleShrt=="AUClwz",        "AssocTableDepVDescCol"] <- "AUC glucose"
    int2[int2$VarGleShrt=="Ogtt2hG",       "AssocTableDepVDescCol"] <- "2h postload glucose"
    int2[int2$VarGleShrt=="FPG",       "AssocTableDepVDescCol"] <- "Fasting glucose"
    
    int2[int2$VarGleShrt=="GDM_iadpsg",     "DepvOrderinT2"] <- 4
    int2[int2$VarGleShrt=="AUClwz",  "DepvOrderinT2"] <- 3
    int2[int2$VarGleShrt=="Ogtt2hG", "DepvOrderinT2"] <- 2
    int2[int2$VarGleShrt=="FPG", "DepvOrderinT2"] <- 1
    int2[int2$VarGleShrt=="GDM",     "DepvOrderinT2"] <- 5
    
    int2 <-  int2[,c("DepVarGle","AssocTableDepVDescCol","DepvOrderinT2", "VarGleShrt")] 
    colnames(int2)[colnames(int2) == "VarGleShrt" ] <-  "DepVarGleShrt"
    
    indepint2<- cvt[cvt$VarGleShrt %in% c("T2dMah14PrsZ","T2dMah14PrsZtop10", "Age","BMI_v2", "Education2_n","BornInSaTF",     "Parity","ParHistDiabYN", "lowDQS"), c("Variable_Gle","VarName_START","VarName_BiB", "VarGleShrt")]
    colnames(indepint2)[1]<- "IndepVarSub_gleVars"
    
    indepint2[indepint2$VarGleShrt=="Age",           "AssocTableIndepVDescCol"] <- "Age (per year)"
    indepint2[indepint2$VarGleShrt=="BMI_v2",        "AssocTableIndepVDescCol"] <- "BMI (per kg/m2 point)"
    indepint2[indepint2$VarGleShrt=="Education2_n",  "AssocTableIndepVDescCol"] <- "Education level (per level)"
    indepint2[indepint2$VarGleShrt=="BornInSaTF",      "AssocTableIndepVDescCol"] <- "Born in South Asia (Yes/No)"
    indepint2[indepint2$VarGleShrt=="Parity",        "AssocTableIndepVDescCol"] <- "Parity (per unit increase)"
    indepint2[indepint2$VarGleShrt=="ParHistDiabYN", "AssocTableIndepVDescCol"] <- "Parental history of T2D (Yes/No)"
    indepint2[indepint2$VarGleShrt=="lowDQS",        "AssocTableIndepVDescCol"] <- "Low diet quality (Yes/No)"
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZ",   "AssocTableIndepVDescCol"] <- "PRS (per unit increase)"
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZtop10",   "AssocTableIndepVDescCol"] <- "PRS (Top 10% vs bottom 90%)"
    
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZ",   "IndepVvOrdinT2"] <- 1
    indepint2[indepint2$VarGleShrt=="Age",           "IndepVvOrdinT2"] <- 2
    indepint2[indepint2$VarGleShrt=="BMI_v2",        "IndepVvOrdinT2"] <- 3
    indepint2[indepint2$VarGleShrt=="Education2_n",  "IndepVvOrdinT2"] <- 7
    indepint2[indepint2$VarGleShrt=="BornInSaTF",      "IndepVvOrdinT2"] <- 4
    indepint2[indepint2$VarGleShrt=="Parity",        "IndepVvOrdinT2"] <- 6
    indepint2[indepint2$VarGleShrt=="ParHistDiabYN", "IndepVvOrdinT2"] <- 5
    indepint2[indepint2$VarGleShrt=="lowDQS",        "IndepVvOrdinT2"] <- 8
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZtop10",   "IndepVvOrdinT2"] <- 9
    
    indepint2 <-  indepint2[,c("IndepVarSub_gleVars", "AssocTableIndepVDescCol", "IndepVvOrdinT2", "VarGleShrt")] 
    colnames(indepint2)[colnames(indepint2) == "VarGleShrt" ] <-  "IndepVarGleShrt"
    
    colnames(tbl2o2_stbib2) [colnames(tbl2o2_stbib2) %in%  colnames(int2)]
    colnames(tbl2o2_stbib2) [colnames(tbl2o2_stbib2) %in%  colnames(indepint2)]
    
    tbl2o2_stbib2 <- merge(tbl2o2_stbib2, int2, all=T) 
    tbl2o2_stbib2 <- merge(tbl2o2_stbib2, indepint2, all=T)
    tbl2o2_stbib2 <- tbl2o2_stbib2[order(tbl2o2_stbib2$"DepvOrderinT2", tbl2o2_stbib2$"IndepVvOrdinT2"), ]   
    
    t2p1<- tbl2o2_stbib2[tbl2o2_stbib2$DepVarGleShrt %in% c("AUClwz", "Ogtt2hG", "FPG"), c("AssocTableDepVDescCol",     "AssocTableIndepVDescCol", "Beta [CI]_st", "P_st", "Beta [CI]_bib", "P_bib", "Beta [CI]_MetaA", "MetaA_P", "I2_MetaA",     "QEp_MetaA")]
    
    t2p2<- tbl2o2_stbib2[tbl2o2_stbib2$DepVarGleShrt %in% c("GDM_iadpsg"),c("AssocTableDepVDescCol",     "AssocTableIndepVDescCol", "OR [95%CI]_st", "P_st", "OR [95%CI]_bib", "P_bib", "OR [95%CI]_MetaA", "MetaA_P",     "I2_MetaA", "QEp_MetaA")]
    
    colnames(t2p1)[which(colnames(t2p1)%in% c("Beta [CI]_st","Beta [CI]_bib","Beta [CI]_MetaA"))] <- c("Beta/OR [95%CI]     _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
    colnames(t2p2)[which(colnames(t2p2)%in%c("OR [95%CI]_st","OR [95%CI]_bib","OR [95%CI]_MetaA"))] <- c("Beta/OR [95%CI]     _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
    colnames(t2p1)[!colnames(t2p1) %in% colnames(t2p2)]
    tbl2o2_stbib3 <- rbind(t2p1, t2p2)
  # Save output
    write.table(tbl2o2_stbib3, "Table2.tsv", row.names=F, quote=F, sep="\t")


# Table 3: PARs 

  # Calculate PARs in each study
    sbl[[1]]$T2dMah14PrsT12vT3 <- as.numeric(sbl[[1]]$T2dMah14PrsT12vT3)
    sbl[[1]]$BornInSAtf <- as.numeric(sbl[[1]]$BornInSAtf)
    sbl[[2]]$BornInSA_AL <- as.numeric(sbl[[2]]$BornInSA_AL)
    
    parVars_d <- cvt[cvt$"VarGleShrt" %in% c("GDM_iadpsg", "ParHistDiabYN", "BornInSaTF", "bmi_s23vle23_v2", "parity_2cat",    "Education_2c","age_s29vle2931", "age_s29vle32", "lowDQS", "BMI_v2", "Age", "T2dMah14PrsZ", "Parity", "Education2_n",    "T2dMah14PrsT12vT3", "T2dMah14PrsZtop10", paste0("PC", 1:5)),c( "Variable_Gle", "VarName_START", "VarName_BiB",      "VarGleShrt")]
    
    par_t <- as.data.frame(matrix(ncol=9, nrow=1))
    par_t <- par_t[-1,]
    colnames(par_t) <- c( "AF", "SE","z" ,"P","Study","Var","Formula", "Lower limit", "Upper limit")
    
    parInfo_d <- as.data.frame(cbind(c("anygdm ~ T2dMah14PrsT12vT3 + agemom + prepregbmi + BornInSAtf + mbeyself + famhxdmi     + dietq_lo",
    "anygdm ~ T2dMah14PrsZtop10 + agemom + prepregbmi + BornInSAtf + mbeyself + famhxdmi + dietq_lo","anygdm ~ T2dMah14PrsZ     + agemom + prepregbmi + BornInSAtf + mbeyself + famhxdmi + parity + dietq_lo", " anygdm ~ T2dMah14PrsZ+ prepregbmi +     BornInSAtf + mbeyself + famhxdmi + dietq_lo + age_s29vle2931", "anygdm ~ T2dMah14PrsZ+ prepregbmi + BornInSAtf +     mbeyself + famhxdmi + dietq_lo + age_s29vle32", "anygdm ~ T2dMah14PrsZ + agemom + BornInSAtf + mbeyself + famhxdmi +     dietq_lo + mbBMIs23vle23", "anygdm ~ T2dMah14PrsZ + agemom + prepregbmi + BornInSAtf + mbeyself2cse2v3p + famhxdmi +     dietq_lo", "anygdm ~ T2dMah14PrsZ + agemom + prepregbmi + BornInSAtf + mbeyself + famhxdmi + dietq_lo", "anygdm ~     T2dMah14PrsZ + agemom + prepregbmi + BornInSAtf + mbeyself + famhxdmi + dietq_lo", "GDMsaBiBwBfnT2D ~ T2dMah14PrsT12vT3     + agemy_mbqall + mms0mbkbmi + BornInSA_AL + edu0mumede4c1234 + ParHistDiabYNMaxValBFB", "GDMsaBiBwBfnT2D ~     T2dMah14PrsZtop10 + agemy_mbqall + mms0mbkbmi + BornInSA_AL + edu0mumede4c1234 + ParHistDiabYNMaxValBFB","     GDMsaBiBwBfnT2D ~ T2dMah14PrsZ + agemy_mbqall + mms0mbkbmi + BornInSA_AL + edu0mumede4c1234 + ParHistDiabYNMaxValBFB","     GDMsaBiBwBfnT2D ~ T2dMah14PrsZ+ mms0mbkbmi + BornInSA_AL + edu0mumede4c1234 + ParHistDiabYNMaxValBFB + age_s29vle2931",     "GDMsaBiBwBfnT2D ~ T2dMah14PrsZ+ mms0mbkbmi + BornInSA_AL + edu0mumede4c1234 + ParHistDiabYNMaxValBFB + age_s29vle32",     "GDMsaBiBwBfnT2D ~ T2dMah14PrsZ + agemy_mbqall + BornInSA_AL + edu0mumede4c1234 + ParHistDiabYNMaxValBFB +     BMIs23vle23v2", "GDMsaBiBwBfnT2D ~ T2dMah14PrsZ + agemy_mbqall + mms0mbkbmi + BornInSA_AL + edu0mumede2cse3v4 +     ParHistDiabYNMaxValBFB", "GDMsaBiBwBfnT2D ~ T2dMah14PrsZ + agemy_mbqall + mms0mbkbmi + BornInSA_AL + edu0mumede4c1234 +     ParHistDiabYNMaxValBFB"),
    c("T2dMah14PrsT12vT3", "T2dMah14PrsZtop10", "famhxdmi","age_s29vle2931", "age_s29vle32", "mbBMIs23vle23",     "mbeyself2cse2v3p", "BornInSAtf", "dietq_lo", "T2dMah14PrsT12vT3", "T2dMah14PrsZtop10", "ParHistDiabYNMaxValBFB",    "age_s29vle2931", "age_s29vle32", "BMIs23vle23v2", "edu0mumede2cse3v4", "BornInSA_AL"), 
    c("PRS T12vsT3", "PRS_Top10vB90", "FAMH", "Age <31","Age>32" ,"BMI s23vle23", "Education", "Born in SA", "dietqLo", "PRS     T12vsT3", "PRS_Top10vB90", "FAMH", "Age <31","Age>32", "BMI s23vle23", "Education", "Born in SA" ),
    c(rep("START",9 ), rep("BiB", 8))), stringsAsFactors=F)
    colnames(parInfo_d) <- c("PafFrmla","PafVar","PafVarGle", "Study") 
    
    parVars_d <- cvt[cvt$"VarGleShrt" %in% c("GDM_iadpsg", "ParHistDiabYN", "BornInSaTF", "bmi_s23vle23_v2", "parity_2cat",    "Education_2c","age_s29vle2931", "age_s29vle32", "lowDQS", "BMI_v2", "Age", "T2dMah14PrsZ", "Parity", "Education2_n",    "T2dMah14PrsT12vT3", "T2dMah14PrsZtop10", paste0("PC", 1:5)),c( "Variable_Gle", "VarName_START", "VarName_BiB",      "VarGleShrt")]
    
    par_t <- as.data.frame(matrix(ncol=9, nrow=1))
    par_t <- par_t[-1,]
    colnames(par_t) <- c( "AF", "SE","z" ,"P","Study","Var","Formula", "Lower limit", "Upper limit")
    
      
    for (std in 1:2){
    
      if(std ==1) parInfo_dsub <- parInfo_d[parInfo_d$Study=="START",]
      if(std ==2) parInfo_dsub <- parInfo_d[parInfo_d$Study=="BiB",]
      
      for (j in 1:nrow(parInfo_dsub)){
        parm_s1 <- glm(as.formula(parInfo_dsub[j,1] ), data=sbl[[std]], family=binomial)
    	afmod <- summary(AFglm(object=parm_s1, exposure=parInfo_dsub[j,2] , data=sbl[[std]]))
        AFj <- as.data.frame(afmod$AF)
        colnames(AFj) <- c("AF", "SE", "z", "P")
        AFj$Study<- parInfo_dsub[j,4]
        AFj$Var <- parInfo_dsub[j,3]
        AFj$Formula <- paste0(afmod$"formula"[2], afmod$"formula"[1],afmod$"formula"[3])
        AFj$"Lower limit" <- afmod$"confidence.interval" [1]
        AFj$"Upper limit" <- afmod$"confidence.interval" [2]
        par_t <- rbind(par_t, AFj)
      }
      
      
      
    }





  # Merge START and BiB association results
    s <- par_t[par_t$Study== "START", ]
    b <- par_t[par_t$Study== "BiB", ]
    colnames(s) <- paste0(colnames(s), "_st")
    colnames(b) <- paste0(colnames(b), "_bib")
    colnames(s) [which(colnames(s) %in% c("Var_st"))] <- c("Var")
    colnames(b) [which(colnames(b) %in% c("Var_bib"))] <- c("Var")
    m <- merge(s, b, by=c("Var"), all=T)
    m[,! colnames(m) %in% c("Formula_st","Formula_bib")]

  # Run Meta-Analysis (fixed effects)

    metaAF<-list()
    metaAMeth <- "FE"
    for( i in 1:nrow(m)){
        metaAF[[i]] <- list()
    	if(!is.na(m[i, "AF_st"]) & !is.na(m[i, "AF_bib"])){
    
    	esci <- escalc(yi=as.numeric(m[i,paste0("AF", c("_st", "_bib"))]), sei=as.numeric(m[i,paste0("SE", c("_st",     "_bib"))]), measure="SMD")
    
    	mai <- rma(yi, vi, method=metaAMeth, data=esci)
    
    	metaAF[[i]] <- mai
    
    	m[i,"MetaA_AF"]  <- mai$"beta"
        m[i,"MetaA_SE"]    <- mai$"se"
        m[i,"MetaA_P"]     <- mai$"pval"
        m[i,"MetaA_I2"]    <- mai$"I2"
        m[i,"MetaA_QE"]    <- mai$"QE"
        m[i,"MetaA_QEp"]   <- mai$"QEp"
        m[i,"MetaA_Model"] <- metaAMeth
    	m[i,"MetaA_CILB"] <- mai$ci.lb
    	m[i,"MetaA_CIUB"] <- mai$ci.ub
    
        }
      }



  # Format final output table

    m2 <- m[,c("Var", "MetaA_Model")]
    m2$"AF(SE)_st" <- paste0(round(m$"AF_st"*100,1), " (", round(m$"SE_st"*100,1), ")")
    
    m2$"AF[CI]_st" <- paste0(round(m$"AF_st"*100,1), " [", round(m$"Lower limit_st"*100,1), " - ", round(m$"Upper     limit_st"*100,1), "]")
    
    m2$"P_st" <- PGPmisc:::formatPvals(m$"P_st")
    m2$"AF(SE)_bib" <- paste0(round(m$"AF_bib"*100,1), " (", round(m$"SE_bib"*100,1), ")")
    
    m2$"AF[CI]_bib" <- paste0(round(m$"AF_bib"*100,1), " [", round(m$"Lower limit_bib"*100,1), " - ", round(m$"Upper     limit_bib"*100,1), "]")
    
    m2[which(!is.na(m$P_bib)), "P_bib"] <- PGPmisc:::formatPvals(m[which(!is.na(m$P_bib)), "P_bib"])
    
    m2$"AF(SE)_MetaA" <- paste0(round(m$"MetaA_AF"*100,1), " (", round(m$"MetaA_SE"*100,1), ")")
    
    m2$"AF[CI]_MetaA" <- paste0(round(m$"MetaA_AF"*100,1), " [", round(m$"MetaA_CILB"*100,1), " - ", round    (m$"MetaA_CIUB"*100,1), "]")
    
    
    m2[which(!is.na(m$MetaA_P)), "MetaA_P"] <- PGPmisc:::formatPvals(m[which(!is.na(m$MetaA_P)), "MetaA_P"])
    m2$"I2_MetaA" <- round(m$"MetaA_I2",0)
    m2[which(!is.na(m$"MetaA_QEp")), "MetaA_QEp"] <- PGPmisc:::formatPvals(m[which(!is.na(m$"MetaA_QEp")), "MetaA_QEp"])
    
    m2[m2$"AF[CI]_bib"=="NA [NA - NA]","AF[CI]_bib"]<-NA
    m2[m2$"AF[CI]_st"=="NA [NA - NA]","AF[CI]_st"]<-NA
    m2[m2$"AF[CI]_MetaA"=="NA [NA - NA]","AF[CI]_MetaA"]<-NA
    
    m2$MetaA_Model <- NULL
    m2$order <- c(1:4, 9, 5:8) 
    m3 <- m2[order(m2$order),]
    m3$"Independent Variable" <- c("Age (29-31 vs. <29 yrs)", "Age (>32yr vs. <29 yrs)", "Body Mass Index (≥ 23 vs. < 23)",   "Born in South Asia (Yes vs. No)", "Education (Post-secondary vs. less)", "Parental history of T2D (Yes vs. No)", "PRS   (Tertile 1+2 vs. 3)", "PRS (Top 10% vs. Bottom 90%)", "Low Diet Quality (Yes vs. No)")

  # Save outputs

    write.table(m, "Table3.tsv", row.names=F, quote=F, sep="\t")
  
    write.table(m3[,c("Independent Variable", "AF[CI]_st", "P_st", "AF[CI]_bib", "P_bib", "AF[CI]_MetaA", "MetaA_P", "I2_MetaA", "MetaA_QEp")] , "/Table3_FullFormatted.tsv", row.names=F, quote=F, sep="\t")





# Table 4: Interactions PRS * RFs with adjustments  
  #Remark: Steps are similar to Table 2
  # Prep tables needed for the analysis
    
    tbl3_depvs_d <- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg", "AUClwz", "Ogtt2hG", "FPG", "Ogtt1hG"), ]
    tbl3_indepvs_d <- cvt[cvt$VarGleShrt =="T2dMah14PrsZ",  ]
    tbl3_intervs_d <- cvt[cvt$VarGleShrt %in% c("Age","BMI_v2", "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "lowDQS"),  ]
    
    tbl3_adjs_l <- list()
    tbl3_adjs_l[[1]] <- cvt[cvt$VarGleShrt %in% c("Age","BMI_v2", paste0("PC", 1:5), "Education2_n","BornInSaTF",   "Parity", "ParHistDiabYN", "lowDQS"), ]
    names(tbl3_adjs_l)<-paste0("model", 1:length(tbl3_adjs_l))
    tbl3pref <- "Table4_AssocInteractwAdj"

  #  Run association tests 

    tbl3o1l <- PGPmisc:::AssocInterMultiStudies(data=sbl, depVs_d=tbl3_depvs_d, indepVs_d=tbl3_indepvs_d, adjVs_l=tbl3_adjs_l,depVs_GleName_Coli=1, indepVs_GleName_Coli=1,adjVs_GleName_Coli=1, depVs_ColNames_Coli=2:3, indepVs_ColNames_Coli=2:3, depVs_ModelTypes_Coli=7, adjVs_ColNames_Coli=2:3, interVs_d=tbl3_intervs_d, interVs_GleName_Coli=1, interVs_ColNames_Coli=2:3,   saveShortOutput=T, saveLongOutput=T, OutputPath=o2_p, OutFileNamePrefx=tbl3pref, OutFileNameSufx=paste0("_", lmd))[[1]]



  # Horizontally merge START and BiB association test results

    tbl3o1l_st <- tbl3o1l[tbl3o1l$Study== "START", ]
    tbl3o1l_bib <- tbl3o1l[tbl3o1l$Study== "BiB", ]
    tbl3colsToKeep<- c("DepVarGle", "IndepVarIGle", "InterVarIGle", "IndepVarSub_gleVars")
    
    colnames(tbl3o1l_st)[which(!colnames(tbl3o1l_st) %in% tbl3colsToKeep)] <- paste0(colnames(tbl3o1l_st)[which(!colnames    (tbl3o1l_st) %in% tbl3colsToKeep)], "_st")
    
    colnames(tbl3o1l_bib)[which(!colnames(tbl3o1l_bib) %in% tbl3colsToKeep)] <- paste0(colnames(tbl3o1l_bib)[which  (!colnames  (tbl3o1l_bib) %in% tbl3colsToKeep)], "_bib")
  
    tbl3o2_stbib <- merge(tbl3o1l_st, tbl3o1l_bib, tbl3colsToKeep, suf=c("_st","_bib"), all=T)
    
    # Visually check the merging 
    checkmt<-unique(tbl3o2_stbib[,c("DepVarGle","IndepVarIGle","InterVarIGle","IndepVarSub_gleVars","IndepVarSub_st",  "IndepVarSub_bib")])
    checkmt
    # merging is ok ! 




  # Remove lowDQS rows in tbl3o2_stbib before meta-A 
    # lowDQS is removed at this stage because it is only tested for interaction in START
    tbl3o2_stbib_2 <- tbl3o2_stbib[-which(tbl3o2_stbib$InterVarIGle == "Low diet quality score (vs Medium + High)"),]

  # Run Meta-Analysis (Fixed effects) 
    tbl3o2_l<-list()
    metaAMeth <- "FE"
    for( i in 1:nrow(tbl3o2_stbib_2)){
    tbl3o2_l[[i]] <- list()
	if(!is.na(tbl3o2_stbib_2[i, "Beta_st"]) & !is.na(tbl3o2_stbib_2[i, "Beta_bib"])){

	  esci2 <- escalc(yi=as.numeric(tbl3o2_stbib_2[i,paste0("Beta", c("_st", "_bib"))]), sei=as.numeric(tbl3o2_stbib_2[i,paste0("SE", c("_st", "_bib"))]), measure="SMD")

	  mai2 <- rma(yi, vi, method=metaAMeth, data=esci2)

	  tbl3o2_l[[i]] <- mai2
  
	  tbl3o2_stbib_2[i,"MetaA_Beta"]  <- mai2$"beta"
      tbl3o2_stbib_2[i,"MetaA_SE"]    <- mai2$"se"
  
      tbl3o2_stbib_2[i,"MetaA_Beta95%CiLB"]    <- mai2$"ci.lb"
      tbl3o2_stbib_2[i,"MetaA_Beta95%CiUB"]    <- mai2$"ci.ub"
  
  
      tbl3o2_stbib_2[i,"MetaA_P"]     <- mai2$"pval"
      tbl3o2_stbib_2[i,"MetaA_I2"]    <- mai2$"I2"
      tbl3o2_stbib_2[i,"MetaA_QE"]    <- mai2$"QE"
      tbl3o2_stbib_2[i,"MetaA_QEp"]   <- mai2$"QEp"
      tbl3o2_stbib_2[i,"MetaA_Model"] <- metaAMeth
      if(!is.na(tbl3o2_stbib_2[i, "OR_st"])){
        tbl3o2_stbib_2[i,"MetaA_OR"]         <- exp(mai2$"beta")
	    tbl3o2_stbib_2[i,"MetaA_95%CiLB"]         <- exp(mai2$"ci.lb")
	    tbl3o2_stbib_2[i,"MetaA_95%CiUB"]         <- exp(mai2$"ci.ub")
      }
    }
    }

  # Add back dietq_lo to the START model 
  
    # Rq: dietq_lo was only adjuted for in START, and hence association results between dietq_lo was not included in the meta-analysis. I want to add it back to the  meta-A table 
    # prep the dietqlo data before I add it back 
    st_dql <- tbl3o2_stbib[which(tbl3o2_stbib$InterVarIGle == "Low diet quality score (vs Medium + High)"),]
  
    tbl3o2_stbib_3 <- merge(tbl3o2_stbib_2,st_dql, all=T)
    head(tbl3o2_stbib_3)


  # Format Output table

    # Part 1 :  
      # add () and [], round values, replace NA[NA-NA] etc ...
      tbl3o2_stbib2 <- tbl3o2_stbib_3[,tbl3colsToKeep]
      tbl3o2_stbib2$"Beta (SE)_st" <- paste0(round(tbl3o2_stbib_3$"Beta_st",3), " (", round(tbl3o2_stbib_3$"SE_st",3), ")")
      
      tbl3o2_stbib2$"Beta [95%CI]_st" <- paste0(round(tbl3o2_stbib_3$"Beta_st",3), " [", round(tbl3o2_stbib_3$"Beta95%CiLB_st",3), " - ", round(tbl3o2_stbib_3$"Beta95%CiUB_st",3), "]")
      
      tbl3o2_stbib2$"OR [95%CI]_st" <- paste0(round(tbl3o2_stbib_3$"OR_st",2), " [", round(tbl3o2_stbib_3$"95%CiLB_st" ,2), "-", round(tbl3o2_stbib_3$"95%CiUB_st",2),"]")
      tbl3o2_stbib2$"P_st" <- PGPmisc:::formatPvals(tbl3o2_stbib_3$"P_st")
      
      tbl3o2_stbib2$"Beta (SE)_bib" <- paste0(round(tbl3o2_stbib_3$"Beta_bib",3), " (", round(tbl3o2_stbib_3$"SE_bib",3), ")")
      
      tbl3o2_stbib2$"Beta [95%CI]_bib" <- paste0(round(tbl3o2_stbib_3$"Beta_bib",3), " [", round(tbl3o2_stbib_3$"Beta95%CiLB_bib",3), " - ", round(tbl3o2_stbib_3$"Beta95%CiUB_bib",3), "]")
      
      tbl3o2_stbib2$"OR [95%CI]_bib" <- paste0(round(tbl3o2_stbib_3$"OR_bib",2), " [", round(tbl3o2_stbib_3$"95%CiLB_bib" ,2), "-", round(tbl3o2_stbib_3$"95%CiUB_bib",2),"]")
      #some rows have NA P values for BiB because dietq_lo is not adjuted for in this study   
      w<-which(is.na(tbl3o2_stbib_3$"P_bib"))
      tbl3o2_stbib2[-w, "P_bib"] <- PGPmisc:::formatPvals(tbl3o2_stbib_3[-w,"P_bib"])
      
      
      tbl3o2_stbib2[-w,"Beta (SE)_MetaA"] <- paste0(round(tbl3o2_stbib_3[-w,"MetaA_Beta"],3), " (", round(tbl3o2_stbib_3[-w,"MetaA_SE"],3), ")")
      
      tbl3o2_stbib2[-w,"Beta [95%CI]_MetaA"] <- paste0(round(tbl3o2_stbib_3[-w,"MetaA_Beta"],3), " [", round(tbl3o2_stbib_3[-w,"MetaA_Beta95%CiLB"],3), " - ", round(tbl3o2_stbib_3[-w,"MetaA_Beta95%CiUB"],3) ,"]")
      
      tbl3o2_stbib2[-w,"OR [95%CI]_MetaA"] <- paste0(round(tbl3o2_stbib_3[-w,"MetaA_OR"],2), " [", round(tbl3o2_stbib_3[-w,"MetaA_95%CiLB"],2), "-", round(tbl3o2_stbib_3[-w,"MetaA_95%CiUB"],2),"]")
      tbl3o2_stbib2[-w,"MetaA_P"] <- PGPmisc:::formatPvals(tbl3o2_stbib_3[-w,"MetaA_P"])
      tbl3o2_stbib2[-w,"I2_MetaA"] <- round(tbl3o2_stbib_3[-w,"MetaA_I2"],0)
      tbl3o2_stbib2[-w,"QEp_MetaA"] <- PGPmisc:::formatPvals(tbl3o2_stbib_3[-w,"MetaA_QEp"])
      
      tbl3o2_stbib2[which(tbl3o2_stbib2$"OR [95%CI]_MetaA"=="NA [NA-NA]"),"OR [95%CI]_MetaA"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"OR [95%CI]_st"=="NA [NA-NA]"),"OR [95%CI]_st"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"OR [95%CI]_bib"=="NA [NA-NA]"),"OR [95%CI]_bib"]<-NA
      
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta (SE)_bib"=="NA (NA)"),"Beta (SE)_bib"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta (SE)_st"=="NA (NA)"),"Beta (SE)_st"]<-NA
      
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta [95%CI]_bib"=="NA [NA - NA]"),"Beta [95%CI]_bib"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta [95%CI]_st"=="NA [NA - NA]"),"Beta [95%CI]_st"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta [95%CI]_bib"=="NA [NA - NA]"),"Beta [95%CI]_bib"]<-NA
  

  
    # Part 2
      int3<- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg", "AUClwz", "Ogtt2hG", "FPG", "GDM"), c("Variable_Gle", "VarName_START", "VarName_BiB", "TableLegend","VarGleShrt")]
      colnames(int3)[colnames(int3) == "Variable_Gle"]<- "DepVarGle"
      
      int3[int3$VarGleShrt=="GDM", 		"AssocTableDepVDescCol"] <- "GDM (South Asian cutoffs)"
      int3[int3$VarGleShrt=="GDM_iadpsg", 		"AssocTableDepVDescCol"] <- "GDM (IADPSG)"
      int3[int3$VarGleShrt=="Ogtt2hG", 		"AssocTableDepVDescCol"] <- "2h post-load glucose"
      int3[int3$VarGleShrt=="FPG", 		"AssocTableDepVDescCol"] <- "Fasting plasma glucose"
      int3[int3$VarGleShrt=="AUClwz", 		"AssocTableDepVDescCol"] <- "AUC glucose"
      
      int3[int3$VarGleShrt=="FPG", 		"DepvOrderinT3"] <- 1
      int3[int3$VarGleShrt=="Ogtt2hG", 	"DepvOrderinT3"] <- 2
      int3[int3$VarGleShrt=="AUClwz",  	"DepvOrderinT3"] <- 3
      int3[int3$VarGleShrt=="GDM_iadpsg",	"DepvOrderinT3"] <- 4
      int3[int3$VarGleShrt=="GDM", 		"DepvOrderinT3"] <- 5
      colnames(int3)[colnames(int3) == "VarGleShrt"]<- "DepVarGleShrt"
      int3 <-  int3[,c("DepVarGle", "DepvOrderinT3", "AssocTableDepVDescCol", "DepVarGleShrt")] 
      
      
      indepint3<- cvt[cvt$VarGleShrt %in% c("T2dMah14PrsZ", "Age","BMI_v2", "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "lowDQS"), c("Variable_Gle", "VarName_START", "VarName_BiB", "TableLegend","VarGleShrt")]
      colnames(indepint3)[colnames(indepint3) == "Variable_Gle"]<- "InterVarIGle"
      
      indepint3[indepint3$VarGleShrt=="Age",           "AssocTableIndepVDescCol"] <- "Age"
      indepint3[indepint3$VarGleShrt=="BMI_v2",        "AssocTableIndepVDescCol"] <- "BMI"
      indepint3[indepint3$VarGleShrt=="Education2_n",  "AssocTableIndepVDescCol"] <- "Education level"
      indepint3[indepint3$VarGleShrt=="BornInSaTF",      "AssocTableIndepVDescCol"] <- "Born in South Asia (Yes/No)"
      indepint3[indepint3$VarGleShrt=="Parity",        "AssocTableIndepVDescCol"] <- "Parity"
      indepint3[indepint3$VarGleShrt=="ParHistDiabYN", "AssocTableIndepVDescCol"] <- "Parental history of T2D (Yes/No)"
      indepint3[indepint3$VarGleShrt=="lowDQS",        "AssocTableIndepVDescCol"] <- "Low diet quality (Yes/No)"
      indepint3[indepint3$VarGleShrt=="T2dMah14PrsZ",   "AssocTableIndepVDescCol"] <- "PRS"
      
      # Set order of covars in table 
      indepint3[indepint3$VarGleShrt=="T2dMah14PrsZ",  "IndepVvOrdinT3"] <- 1
      indepint3[indepint3$VarGleShrt=="Age",           "IndepVvOrdinT3"] <- 2
      indepint3[indepint3$VarGleShrt=="BMI_v2",        "IndepVvOrdinT3"] <- 3
      indepint3[indepint3$VarGleShrt=="Education2_n",  "IndepVvOrdinT3"] <- 7
      indepint3[indepint3$VarGleShrt=="BornInSaTF",    "IndepVvOrdinT3"] <- 4
      indepint3[indepint3$VarGleShrt=="Parity",        "IndepVvOrdinT3"] <- 6
      indepint3[indepint3$VarGleShrt=="ParHistDiabYN", "IndepVvOrdinT3"] <- 5
      indepint3[indepint3$VarGleShrt=="lowDQS",        "IndepVvOrdinT3"] <- 8
      
      colnames(indepint3)[colnames(indepint3) == "VarGleShrt"]<- "InterVarGleShrt"
      
      indepint3 <-  indepint3[,c("InterVarIGle","AssocTableIndepVDescCol", "IndepVvOrdinT3", "InterVarGleShrt")] 

  
    # Part 3 
      names(table(indepint3$InterVarIGle))
      names(table(tbl3o2_stbib2$InterVarIGle))
      
      names(table(int3$"DepVarGle"))
      names(table(tbl3o2_stbib2$"DepVarGle"))
      
      colnames(tbl3o2_stbib2) [colnames(tbl3o2_stbib2) %in%  colnames(int3)]
      colnames(tbl3o2_stbib2) [colnames(tbl3o2_stbib2) %in%  colnames(indepint3)]
      
      tbl3o2_stbib2 <- merge(tbl3o2_stbib2, int3, all=T) 
      tbl3o2_stbib2 <- merge(tbl3o2_stbib2, indepint3, all=T)
      tbl3o2_stbib2 <- tbl3o2_stbib2[order(tbl3o2_stbib2$"DepvOrderinT3", tbl3o2_stbib2$"IndepVvOrdinT3"), ]   
      tbl3o2_stbib2$interactionTerm <- paste0( "PRS x ", tbl3o2_stbib2$AssocTableIndepVDescCol)
      
      t3p1<- tbl3o2_stbib2[tbl3o2_stbib2$"DepVarGleShrt" %in% c("AUClwz", "Ogtt2hG", "FPG"), c("DepVarGle", "interactionTerm", "Beta [95%CI]_st", "P_st", "Beta [95%CI]_bib", "P_bib", "Beta [95%CI]_MetaA", "MetaA_P", "I2_MetaA", "QEp_MetaA")]
      
      t3p2<- tbl3o2_stbib2[tbl3o2_stbib2$"DepVarGleShrt" %in% c("GDM_iadpsg", "GDM"),c("DepVarGle", "interactionTerm", "OR [95%CI]_st", "P_st", "OR [95%CI]_bib", "P_bib", "OR [95%CI]_MetaA", "MetaA_P", "I2_MetaA", "QEp_MetaA")]
      
      colnames(t3p1)[which(colnames(t3p1) %in% c("Beta [95%CI]_st", "Beta [95%CI]_bib", "Beta [95%CI]_MetaA"))]  <- c("Beta/OR [95%CI] _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
      
      colnames(t3p2)[which(colnames(t3p2) %in% c("OR [95%CI]_st","OR [95%CI]_bib","OR [95%CI]_MetaA"))] <- c("Beta/OR [95%CI] _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
      
      tbl3o2_stbib3 <- rbind(t3p1, t3p2)
      
      tbl3o2_stbib3[which(tbl3o2_stbib3$"Beta/OR [95%CI] _bib"=="NA [NA - NA]"),"Beta/OR [95%CI] _bib"]<-NA


  # Save outputs
    write.table(tbl3o2_stbib3, "Table4.tsv", row.names=F, quote=F, sep="\t")



# Figure 1 - GDM ~ PRS * BMI_tertiles

  bibintervars <- cvt[cvt[,"VarGleShrt"] %in% c("T2dMah14PrsZ","mms0mbkbmi", "Age","BMI_v2","T2dMah14PrsTer" ,paste0("PC", 1:5), "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "GDM_iadpsg"),"VarName_BiB"]
  bibinterd <- sbl[[2]][ , bibintervars] 
  bibinterd <- bibinterd[complete.cases(bibinterd),]
  bibinterd$BMI <- bibinterd$mms0mbkbmi # Here I just rename the variable to faciliated the naming on the plot 
  f1bib <- "drvgesdiab01 ~ BMI * T2dMah14PrsZ + agemy_mbqall "
  jpeg("/Figure1.jpeg", width = 6.5, height = 4, units = "in", res=300)
    m1bib<- glm(f1bib, data=bibinterd, family=binomial)
    interact_plot(m1bib, pred=T2dMah14PrsZ, modx="BMI", x.label="PRS (Z scores)", y.label="Predicted probability of GDM", interval=T, modx.values ="terciles")
    dev.off()

# Figure 2 - FPG ~ PRS * dietq_lo
  stintervars <- cvt[cvt[,"VarGleShrt"] %in% c("FPG","T2dMah14PrsZ","Age","BMI_v2", paste0("PC", 1:5), "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "lowDQS", "GDM_iadpsg"),"VarName_START" ]
  stinterd <- sbl[[1]][ , stintervars] 
  stinterd <- stinterd[complete.cases(stinterd),]

  fst <- "lbrfpgv ~ dietq_lo * T2dMah14PrsZ + agemom + prepregbmi + mbeyself + parity + PC1_v1.0 + PC2_v1.0 + PC3_v1.0 + PC4_v1.0 + PC5_v1.0 +  + BornInSAtf"
  mst<- lm(fst, data=stinterd)
  jpeg(paste(o2_p, "/Figure2.jpeg", sep=""), width = 5, height = 4, units = "in", res=300)
  interact_plot(mst, pred=T2dMah14PrsZ, modx=dietq_lo, interval = TRUE, int.width = 0.8)
  dev.off() 

# Supplementary File 1a (Table S1) : List of SNPs in the PRS

  # Miscellaneous: Make sure the catgorical variable is categorical
    sbl[[1]]$T2dMah14PrsTer <- as.character(sbl[[1]]$T2dMah14PrsTer)
    sbl[[2]]$T2dMah14PrsTer <- as.character(sbl[[2]]$T2dMah14PrsTer)

  # Open the list of SNPs included in the PRS
    common <- read.table(paste0(comsnpinfp_p, "/START_BiB_Mahajan2014_1KG_CommonSNPsOnly_infos.tsv"), head=T, string=F, sep="\t")
  
  # Subset SNPs iof interest, and clean data 
    stbl <- common[common$P<=0.0014,c("RSID_Mah", "CHR", "POS", "RiskAllele", "OtherAllele", "OR", "OR_95L", "OR_95U", "P")]
    colnames(stbl) <- c("SNP ID", "Chromosome", "Position", "Risk Allele", "Other Allele", "OR", "OR_95L", "OR_95U", "P")
    stbl$orci<- paste0( stbl$OR, " [", stbl$OR_95L,"-", stbl$OR_95U, "]")
    stbl2 <-stbl[order(stbl$Chromosome, stbl$Position),]
    stbl2$OR <- NULL
    stbl2$OR_95L <- NULL
    stbl2$OR_95U <- NULL
  # Save output 
    write.csv(stbl2, "SuplementaryFile1a.txt", row.names=F, quote=F)



# Supplementary File 1b (Table S2) : univariate assiocation results between PRS tertiles and traits of interest 
# Supplementary File 1c (Table S3) : GDM SA ~ PRS + Adj  

  # Prep tables for analysis
    tbl2_depvs_d <- cvt[cvt$VarGleShrt %in% c("GDM",), ]
    tbl2_indepvs_d <- cvt[cvt$VarGleShrt %in% c("T2dMah14PrsZ"),]
    
    tbl2_adjs_l <- list()
    tbl2_adjs_l[[1]] <- cvt[cvt$VarGleShrt %in% c("Age","BMI_v2", paste0("PC", 1:5), "Education2_n","BornInSaTF",   "Parity",  "ParHistDiabYN", "lowDQS"), ]
    names(tbl2_adjs_l)<-paste0("model", 1:length(tbl2_adjs_l))
    tbl2pref <- "Table2_GDMiadpsg_AssocwAdj"

  # Run association tests 
    #This section runs association tests on all the variables of interest in each study separately 
    tbl2o1l <- PGPmisc:::AssocInterMultiStudies(data=sbl, depVs_d=tbl2_depvs_d, indepVs_d=tbl2_indepvs_d, adjVs_l=tbl2_adjs_l,depVs_GleName_Coli=1, indepVs_GleName_Coli=1,adjVs_GleName_Coli=1, depVs_ColNames_Coli=2:3, indepVs_ColNames_Coli=2:3, depVs_ModelTypes_Coli=7, adjVs_ColNames_Coli=2:3,  saveShortOutput=F, saveLongOutput=T, OutputPath=o2_p, OutFileNamePrefx=tbl2pref, OutFileNameSufx=paste0("_", lmd)) [[2]]

  # Horizontally merge START and BiB association test results
  
    tbl2o1l_st <- tbl2o1l[tbl2o1l$Study== "START", ]
    tbl2o1l_bib <- tbl2o1l[tbl2o1l$Study== "BiB", ]
    tbl2colsToKeep<- c("DepVarGle", "IndepVarIGle","IndepVarSub_gleVars", "Model", "Ndiai")
    
    colnames(tbl2o1l_st)[which(!colnames(tbl2o1l_st) %in% tbl2colsToKeep)] <- paste0(colnames(tbl2o1l_st)[which(!colnames    (tbl2o1l_st) %in% tbl2colsToKeep)], "_st")
    colnames(tbl2o1l_bib)[which(!colnames(tbl2o1l_bib) %in% tbl2colsToKeep)] <- paste0(colnames(tbl2o1l_bib)[which  (!colnames  (tbl2o1l_bib) %in% tbl2colsToKeep)], "_bib")
  
    tbl2o2_stbib <- merge(tbl2o1l_st, tbl2o1l_bib, tbl2colsToKeep, suf=c("_st","_bib"), all=T)
    #visually check the merging 
    checkmt<-unique(tbl2o2_stbib[tbl2o2_stbib$IndepVarSub_gleVars != "(Intercept)", c("IndepVarSub_gleVars",     "IndepVarSub_bib", "IndepVarSub_st")])
    checkmt
    # Visually, merging is ok ! 

  # Remove rows of (Intercept) and PCA, and  lowDQS rows in tbl2o2_stbib before meta-A 
  
    # lowDQS is removed because it is only adjusted for in START
    # PCAs are removed because they are different between studies 
    tbl2o2_stbib_2 <- tbl2o2_stbib[-which(tbl2o2_stbib$IndepVarSub_gleVars %in% c("(Intercept)", paste0("PC", 1:5),   "lowDQS")),]


  # Run meta-analysis FE
    #used fixed effects models

    tbl2o2_l<-list()
    metaAMeth <- "FE"
    for( i in 1:nrow(tbl2o2_stbib_2)){
        tbl2o2_l[[i]] <- list()
    	if(!is.na(tbl2o2_stbib_2[i, "Beta_st"]) & !is.na(tbl2o2_stbib_2[i, "Beta_bib"])){
    
    	esci2 <- escalc(yi=as.numeric(tbl2o2_stbib_2[i,paste0("Beta", c("_st", "_bib"))]), sei=as.numeric(tbl2o2_stbib_2[i,    paste0("SE", c("_st", "_bib"))]), measure="SMD")
    
    	mai2 <- rma(yi, vi, method=metaAMeth, data=esci2)
    
    	tbl2o2_l[[i]] <- mai2
    
    	tbl2o2_stbib_2[i,"MetaA_Beta"]  <- mai2$"beta"
        tbl2o2_stbib_2[i,"MetaA_SE"]    <- mai2$"se"
        tbl2o2_stbib_2[i,"MetaA_Beta95%CiLB"] <- mai2$"ci.lb"
        tbl2o2_stbib_2[i,"MetaA_Beta95%CiUB"] <- mai2$"ci.ub"
        tbl2o2_stbib_2[i,"MetaA_P"]     <- mai2$"pval"
        tbl2o2_stbib_2[i,"MetaA_I2"]    <- mai2$"I2"
        tbl2o2_stbib_2[i,"MetaA_QE"]    <- mai2$"QE"
        tbl2o2_stbib_2[i,"MetaA_QEp"]   <- mai2$"QEp"
        tbl2o2_stbib_2[i,"MetaA_Model"] <- metaAMeth
        	
    	if(!is.na(tbl2o2_stbib_2[i, "OR_st"])){
          tbl2o2_stbib_2[i,"MetaA_OR"]         <- exp(mai2$"beta")
    	  tbl2o2_stbib_2[i,"MetaA_95%CiLB"]         <- exp(mai2$"ci.lb")
    	  tbl2o2_stbib_2[i,"MetaA_95%CiUB"]         <- exp(mai2$"ci.ub")
        }
      }
    }


  # Add back dietq_lo to the START model infos 
    
    dietq_lo was only adjuted for in START, and hence association results between dietq_lo was not included in the     eta-analysis. I want to add it back to the  meta-A table 
    
    # prep the dietqlo data before I add it back 
    t_dql <- tbl2o1l_st[tbl2o1l_st$IndepVarSub_st ==  "dietq_lo",]
    
    bl2o2_stbib_3 <- merge(tbl2o2_stbib_2,st_dql, all=T)
  

  # Format output table

    # use Beta [95%CIL - 95%CIU] format
    # round P values
    # order the rows inthe  table
    # replace NA[NA-NA] by NA
    # ... etc 
    
    tbl2o2_stbib2 <- tbl2o2_stbib_3[,tbl2colsToKeep]
    tbl2o2_stbib2$"Beta (SE)_st" <- paste0(round(tbl2o2_stbib_3$"Beta_st",3), " (", round(tbl2o2_stbib_3$"SE_st",3), ")")
    
    tbl2o2_stbib2$"Beta [CI]_st" <- paste0(round(tbl2o2_stbib_3$"Beta_st",3), " [", round(tbl2o2_stbib_3$"Beta95%CiLB_st",3),    " - ",  round(tbl2o2_stbib_3$"Beta95%CiUB_st",3), "]")
    
    tbl2o2_stbib2$"OR [95%CI]_st" <- paste0(round(tbl2o2_stbib_3$"OR_st",2), " [", round(tbl2o2_stbib_3$"95%CiLB_st" ,2), "     - ", round(tbl2o2_stbib_3$"95%CiUB_st",2),"]")
    tbl2o2_stbib2$"P_st" <- PGPmisc:::formatPvals(tbl2o2_stbib_3$"P_st")
    
    
    tbl2o2_stbib2$"Beta (SE)_bib" <- paste0(round(tbl2o2_stbib_3$"Beta_bib",3), " (", round(tbl2o2_stbib_3$"SE_bib",3), ")")
    
    tbl2o2_stbib2$"Beta [CI]_bib" <- paste0(round(tbl2o2_stbib_3$"Beta_bib",3), " [", round(tbl2o2_stbib_3$"Beta95%CiLB_bib",    3)," - ",  round(tbl2o2_stbib_3$"Beta95%CiUB_bib",3), "]")
    
    tbl2o2_stbib2$"OR [95%CI]_bib" <- paste0(round(tbl2o2_stbib_3$"OR_bib",2), " [", round(tbl2o2_stbib_3$"95%CiLB_bib" ,2),     " - ", round(tbl2o2_stbib_3$"95%CiUB_bib",2),"]")
    #some rows have NA P values for BiB because dietq_lo is not adjuted for in this study   
    w<-which(is.na(tbl2o2_stbib_3$"P_bib"))
    tbl2o2_stbib2[-w, "P_bib"] <- PGPmisc:::formatPvals(tbl2o2_stbib_3[-w,"P_bib"])
    
    
    tbl2o2_stbib2[-w,"Beta (SE)_MetaA"] <- paste0(round(tbl2o2_stbib_3[-w,"MetaA_Beta"],3), " (", round(tbl2o2_stbib_3[-w,    "MetaA_SE"],3), ")")
    
    tbl2o2_stbib2[-w,"Beta [CI]_MetaA"] <- paste0(round(tbl2o2_stbib_3[-w,"MetaA_Beta"],3), " [", round(tbl2o2_stbib_3[-w,    "MetaA_Beta95%CiLB"],3), " - ",  round(tbl2o2_stbib_3[-w,"MetaA_Beta95%CiUB"],3), "]")
    
    
    
    tbl2o2_stbib2[-w,"OR [95%CI]_MetaA"] <- paste0(round(tbl2o2_stbib_3[-w,"MetaA_OR"],2), " [", round(tbl2o2_stbib_3[-w,    "MetaA_95%CiLB"],2), " - ", round(tbl2o2_stbib_3[-w,"MetaA_95%CiUB"],2),"]")
    
    tbl2o2_stbib2[-w,"MetaA_P"] <- PGPmisc:::formatPvals(tbl2o2_stbib_3[-w,"MetaA_P"])
    tbl2o2_stbib2[-w,"I2_MetaA"] <- round(tbl2o2_stbib_3[-w,"MetaA_I2"],0)
    tbl2o2_stbib2[-w,"QEp_MetaA"] <- PGPmisc:::formatPvals(tbl2o2_stbib_3[-w,"MetaA_QEp"])
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"OR [95%CI]_MetaA"=="NA [NA - NA]"),"OR [95%CI]_MetaA"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"OR [95%CI]_st"=="NA [NA - NA]"),"OR [95%CI]_st"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"OR [95%CI]_bib"=="NA [NA - NA]"),"OR [95%CI]_bib"]<-NA
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta (SE)_bib"=="NA (NA)"),"Beta (SE)_bib"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta (SE)_st"=="NA (NA)"),"Beta (SE)_st"]<-NA
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta [CI]_bib"=="NA [NA - NA]"),"Beta [CI]_bib"]<-NA
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta [CI]_st"=="NA [NA - NA]"),"Beta [CI]_st"]<-NA
    
    tbl2o2_stbib2[which(tbl2o2_stbib2$"Beta [CI]_MetaA"=="NA [NA - NA]"),"Beta [CI]_MetaA"]<-NA
    
    
    int2<- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg", "AUClwz", "Ogtt2hG", "FPG"), c(1:3, which(colnames(cvt) ==     "VarGleShrt")) ]
    
    colnames(int2)[1]<- "DepVarGle"
    int2[int2$VarGleShrt=="GDM_iadpsg",           "AssocTableDepVDescCol"] <- "GDM (iadpsg)"
    int2[int2$VarGleShrt=="GDM",           "AssocTableDepVDescCol"] <- "GDM (SA specific)"
    int2[int2$VarGleShrt=="AUClwz",        "AssocTableDepVDescCol"] <- "AUC glucose"
    int2[int2$VarGleShrt=="Ogtt2hG",       "AssocTableDepVDescCol"] <- "2h postload glucose"
    int2[int2$VarGleShrt=="FPG",       "AssocTableDepVDescCol"] <- "Fasting glucose"
    
    int2[int2$VarGleShrt=="GDM_iadpsg",     "DepvOrderinT2"] <- 4
    int2[int2$VarGleShrt=="AUClwz",  "DepvOrderinT2"] <- 3
    int2[int2$VarGleShrt=="Ogtt2hG", "DepvOrderinT2"] <- 2
    int2[int2$VarGleShrt=="FPG", "DepvOrderinT2"] <- 1
    int2[int2$VarGleShrt=="GDM",     "DepvOrderinT2"] <- 5
    
    int2 <-  int2[,c("DepVarGle","AssocTableDepVDescCol","DepvOrderinT2", "VarGleShrt")] 
    colnames(int2)[colnames(int2) == "VarGleShrt" ] <-  "DepVarGleShrt"
    
    indepint2<- cvt[cvt$VarGleShrt %in% c("T2dMah14PrsZ","T2dMah14PrsZtop10", "Age","BMI_v2", "Education2_n","BornInSaTF",     "Parity","ParHistDiabYN", "lowDQS"), c("Variable_Gle","VarName_START","VarName_BiB", "VarGleShrt")]
    colnames(indepint2)[1]<- "IndepVarSub_gleVars"
    
    indepint2[indepint2$VarGleShrt=="Age",           "AssocTableIndepVDescCol"] <- "Age (per year)"
    indepint2[indepint2$VarGleShrt=="BMI_v2",        "AssocTableIndepVDescCol"] <- "BMI (per kg/m2 point)"
    indepint2[indepint2$VarGleShrt=="Education2_n",  "AssocTableIndepVDescCol"] <- "Education level (per level)"
    indepint2[indepint2$VarGleShrt=="BornInSaTF",      "AssocTableIndepVDescCol"] <- "Born in South Asia (Yes/No)"
    indepint2[indepint2$VarGleShrt=="Parity",        "AssocTableIndepVDescCol"] <- "Parity (per unit increase)"
    indepint2[indepint2$VarGleShrt=="ParHistDiabYN", "AssocTableIndepVDescCol"] <- "Parental history of T2D (Yes/No)"
    indepint2[indepint2$VarGleShrt=="lowDQS",        "AssocTableIndepVDescCol"] <- "Low diet quality (Yes/No)"
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZ",   "AssocTableIndepVDescCol"] <- "PRS (per unit increase)"
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZtop10",   "AssocTableIndepVDescCol"] <- "PRS (Top 10% vs bottom 90%)"
    
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZ",   "IndepVvOrdinT2"] <- 1
    indepint2[indepint2$VarGleShrt=="Age",           "IndepVvOrdinT2"] <- 2
    indepint2[indepint2$VarGleShrt=="BMI_v2",        "IndepVvOrdinT2"] <- 3
    indepint2[indepint2$VarGleShrt=="Education2_n",  "IndepVvOrdinT2"] <- 7
    indepint2[indepint2$VarGleShrt=="BornInSaTF",      "IndepVvOrdinT2"] <- 4
    indepint2[indepint2$VarGleShrt=="Parity",        "IndepVvOrdinT2"] <- 6
    indepint2[indepint2$VarGleShrt=="ParHistDiabYN", "IndepVvOrdinT2"] <- 5
    indepint2[indepint2$VarGleShrt=="lowDQS",        "IndepVvOrdinT2"] <- 8
    indepint2[indepint2$VarGleShrt=="T2dMah14PrsZtop10",   "IndepVvOrdinT2"] <- 9
    
    indepint2 <-  indepint2[,c("IndepVarSub_gleVars", "AssocTableIndepVDescCol", "IndepVvOrdinT2", "VarGleShrt")] 
    colnames(indepint2)[colnames(indepint2) == "VarGleShrt" ] <-  "IndepVarGleShrt"
    
    colnames(tbl2o2_stbib2) [colnames(tbl2o2_stbib2) %in%  colnames(int2)]
    colnames(tbl2o2_stbib2) [colnames(tbl2o2_stbib2) %in%  colnames(indepint2)]
    
    tbl2o2_stbib2 <- merge(tbl2o2_stbib2, int2, all=T) 
    tbl2o2_stbib2 <- merge(tbl2o2_stbib2, indepint2, all=T)
    tbl2o2_stbib2 <- tbl2o2_stbib2[order(tbl2o2_stbib2$"DepvOrderinT2", tbl2o2_stbib2$"IndepVvOrdinT2"), ]   
    
    t2p1<- tbl2o2_stbib2[tbl2o2_stbib2$DepVarGleShrt %in% c("AUClwz", "Ogtt2hG", "FPG"), c("AssocTableDepVDescCol",     "AssocTableIndepVDescCol", "Beta [CI]_st", "P_st", "Beta [CI]_bib", "P_bib", "Beta [CI]_MetaA", "MetaA_P", "I2_MetaA",     "QEp_MetaA")]
    
    t2p2<- tbl2o2_stbib2[tbl2o2_stbib2$DepVarGleShrt %in% c("GDM_iadpsg"),c("AssocTableDepVDescCol",     "AssocTableIndepVDescCol", "OR [95%CI]_st", "P_st", "OR [95%CI]_bib", "P_bib", "OR [95%CI]_MetaA", "MetaA_P",     "I2_MetaA", "QEp_MetaA")]
    
    colnames(t2p1)[which(colnames(t2p1)%in% c("Beta [CI]_st","Beta [CI]_bib","Beta [CI]_MetaA"))] <- c("Beta/OR [95%CI]     _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
    colnames(t2p2)[which(colnames(t2p2)%in%c("OR [95%CI]_st","OR [95%CI]_bib","OR [95%CI]_MetaA"))] <- c("Beta/OR [95%CI]     _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
    colnames(t2p1)[!colnames(t2p1) %in% colnames(t2p2)]
    tbl2o2_stbib3 <- rbind(t2p1, t2p2)
  # Save output
    write.csv(tbl2o2_stbib3, "SuplementaryFile1d.csv", row.names=F, quote=F)
  

# Supplementary File 1d (Table S4) : multivariate association results between PRS tertiles and traits of interest 

  # Prep tables for analysis
   
    tbls3_depvs_d <- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg", "AUClwz", "Ogtt2hG", "FPG", "GDM"), ]
    tbls3_indepvs_d <- cvt[cvt$VarGleShrt =="T2dMah14PrsTer",  ]
    
    tbls3_adjs_l <- list()
    tbls3_adjs_l[[1]] <- cvt[cvt$VarGleShrt %in% c("Age","BMI_v2", paste0("PC", 1:5), "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "lowDQS"), ]
    names(tbls3_adjs_l)<-paste0("model", 1:length(tbls3_adjs_l))
    tbls3pref <- "TablS3_AssocwAdj"

  # Run associations 
    
    tbls3o1l <- PGPmisc:::AssocInterMultiStudies(data=sbl, depVs_d=tbls3_depvs_d, indepVs_d=tbls3_indepvs_d, adjVs_l=tbls3_adjs_l,depVs_GleName_Coli=1, indepVs_GleName_Coli=1,adjVs_GleName_Coli=1, depVs_ColNames_Coli=2:3, indepVs_ColNames_Coli=2:3, depVs_ModelTypes_Coli=which(colnames(tbls3_depvs_d)== "Model_Type"), adjVs_ColNames_Coli=2:3,  saveShortOutput=F, saveLongOutput=F) [[1]]
    
  # Horizontally merge START and BiB association tests
    
    tbls3o1l_st <- tbls3o1l[tbls3o1l$Study== "START", ]
    tbls3o1l_bib <- tbls3o1l[tbls3o1l$Study== "BiB", ]
    tbls3colsToKeep<- c("DepVarGle", "IndepVarIGle","IndepVarSub", "IndepVarSub_gleVars", "Model", "Ndiai")
    
    colnames(tbls3o1l_st)[which(!colnames(tbls3o1l_st) %in% tbls3colsToKeep)] <- paste0(colnames(tbls3o1l_st)[which(!colnames(tbls3o1l_st) %in% tbls3colsToKeep)], "_st")
    colnames(tbls3o1l_bib)[which(!colnames(tbls3o1l_bib) %in% tbls3colsToKeep)] <- paste0(colnames(tbls3o1l_bib)[which(!colnames(tbls3o1l_bib) %in% tbls3colsToKeep)], "_bib")
   
    tbls3o2_stbib <- merge(tbls3o1l_st, tbls3o1l_bib, tbls3colsToKeep, suf=c("_st","_bib"), all=T)
    #visually check the merging 
    checkmt<-tbls3o2_stbib[, c("IndepVarSub_gleVars", "IndepVarSub")]
    checkmt
    #merging is ok ! 
    

    
  # Meta-A FE
    
    tbls3o2_l<-list()
    metaAMeth <- "FE"
    for( i in 1:nrow(tbls3o2_stbib)){
        tbls3o2_l[[i]] <- list()
    	if(!is.na(tbls3o2_stbib[i, "Beta_st"]) & !is.na(tbls3o2_stbib[i, "Beta_bib"])){
    
    	esci2 <- escalc(yi=as.numeric(tbls3o2_stbib[i,paste0("Beta", c("_st", "_bib"))]), sei=as.numeric(tbls3o2_stbib[i,paste0("SE", c("_st", "_bib"))]), measure="SMD")
    
    	mai2 <- rma(yi, vi, method=metaAMeth, data=esci2)
    
    	tbls3o2_l[[i]] <- mai2
    
    	tbls3o2_stbib[i,"MetaA_Beta"]  <- mai2$"beta"
        tbls3o2_stbib[i,"MetaA_SE"]    <- mai2$"se"
    	tbls3o2_stbib[i,"MetaA_Beta95%CiLB"]    <- mai2$"ci.lb"
        tbls3o2_stbib[i,"MetaA_Beta95%CiUB"]    <- mai2$"ci.ub"
        tbls3o2_stbib[i,"MetaA_P"]     <- mai2$"pval"
        tbls3o2_stbib[i,"MetaA_I2"]    <- mai2$"I2"
        tbls3o2_stbib[i,"MetaA_QE"]    <- mai2$"QE"
        tbls3o2_stbib[i,"MetaA_QEp"]   <- mai2$"QEp"
        tbls3o2_stbib[i,"MetaA_Model"] <- metaAMeth
        if(!is.na(tbls3o2_stbib[i, "OR_st"])){
          tbls3o2_stbib[i,"MetaA_OR"]         <- exp(mai2$"beta")
    	  tbls3o2_stbib[i,"MetaA_95%CiLB"]         <- exp(mai2$"ci.lb")
    	  tbls3o2_stbib[i,"MetaA_95%CiUB"]         <- exp(mai2$"ci.ub")
        }
      }
    }

    
  # Format output table
    
    # Part 1 

      tbls3o2_stbib2 <- tbls3o2_stbib[,tbls3colsToKeep]
      
      tbls3o2_stbib2$"Beta (SE)_st" <- paste0(round(tbls3o2_stbib$"Beta_st",3), " (", round(tbls3o2_stbib$"SE_st",3), ")")
      
      tbls3o2_stbib2$"Beta [CI]_st" <- paste0(round(tbls3o2_stbib$"Beta_st",3), " [", round(tbls3o2_stbib$"Beta95%CiLB_st",3), " - ", round(tbls3o2_stbib$"Beta95%CiUB_st",3), "]")
      
      tbls3o2_stbib2$"OR [95%CI]_st" <- paste0(round(tbls3o2_stbib$"OR_st",2), " [", round(tbls3o2_stbib$"95%CiLB_st" ,2), "-", round(tbls3o2_stbib$"95%CiUB_st",2),"]")
      tbls3o2_stbib2$"P_st" <- PGPmisc:::formatPvals(tbls3o2_stbib$"P_st")
      tbls3o2_stbib2[, "P_bib"] <- PGPmisc:::formatPvals(tbls3o2_stbib[,"P_bib"])
      
      tbls3o2_stbib2$"Beta (SE)_bib" <- paste0(round(tbls3o2_stbib$"Beta_bib",3), " (", round(tbls3o2_stbib$"SE_bib",3), ")")
      
      tbls3o2_stbib2$"Beta [CI]_bib" <- paste0(round(tbls3o2_stbib$"Beta_bib",3), " [", round(tbls3o2_stbib$"Beta95%CiLB_bib",3), " - ", round(tbls3o2_stbib$"Beta95%CiUB_bib",3), "]")
      
      tbls3o2_stbib2$"OR [95%CI]_bib" <- paste0(round(tbls3o2_stbib$"OR_bib",2), " [", round(tbls3o2_stbib$"95%CiLB_bib" ,2), "-", round(tbls3o2_stbib$"95%CiUB_bib",2),"]")
      
      tbls3o2_stbib2[,"Beta (SE)_MetaA"] <- paste0(round(tbls3o2_stbib[,"MetaA_Beta"],3), " (", round(tbls3o2_stbib[,"MetaA_SE"],3), ")")
      tbls3o2_stbib2[,"Beta [CI]_MetaA"] <- paste0(round(tbls3o2_stbib[,"MetaA_Beta"],3), " [", round(tbls3o2_stbib[,"MetaA_Beta95%CiLB"],3), " - ", round(tbls3o2_stbib[,"MetaA_Beta95%CiUB"],3) , "]")
      
      tbls3o2_stbib2[,"OR [95%CI]_MetaA"] <- paste0(round(tbls3o2_stbib[,"MetaA_OR"],2), " [", round(tbls3o2_stbib[,"MetaA_95%CiLB"],2), "-", round(tbls3o2_stbib[,"MetaA_95%CiUB"],2),"]")
      tbls3o2_stbib2[,"MetaA_P"] <- PGPmisc:::formatPvals(tbls3o2_stbib[,"MetaA_P"])
      tbls3o2_stbib2[,"I2_MetaA"] <- round(tbls3o2_stbib[,"MetaA_I2"],0)
      tbls3o2_stbib2[,"QEp_MetaA"] <- PGPmisc:::formatPvals(tbls3o2_stbib[,"MetaA_QEp"])
      
      tbls3o2_stbib2[tbls3o2_stbib2=="NA [NA-NA]"]<-NA
      tbls3o2_stbib2[tbls3o2_stbib2=="NA [NA - NA]"]<-NA
      tbls3o2_stbib2[tbls3o2_stbib2=="NA (NA)"]<-NA

    
    # Part 2 

    
      ints3<- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg", "AUClwz", "Ogtt2hG", "FPG", "GDM"), c(colnames(cvt)[1:3], "VarGleShrt")]
      colnames(ints3)[1]<- "DepVarGle"
      ints3[ints3$VarGleShrt=="GDM_iadpsg",           "AssocTableDepVDescCol"] <- "GDM (IADPSG)"
      ints3[ints3$VarGleShrt=="AUClwz",        "AssocTableDepVDescCol"] <- "AUC glucose"
      ints3[ints3$VarGleShrt=="Ogtt2hG",       "AssocTableDepVDescCol"] <- "2h postload glucose"
      ints3[ints3$VarGleShrt=="FPG",       "AssocTableDepVDescCol"] <- "Fasting glucose"
      ints3[ints3$VarGleShrt=="GDM",       "AssocTableDepVDescCol"] <- "GDM (South Asian-specific definition)"
      
      ints3[ints3$VarGleShrt=="GDM_iadpsg",     "DepvOrderinST3"] <- 4
      ints3[ints3$VarGleShrt=="AUClwz",  "DepvOrderinST3"] <- 3
      ints3[ints3$VarGleShrt=="Ogtt2hG", "DepvOrderinST3"] <- 2
      ints3[ints3$VarGleShrt=="FPG", "DepvOrderinST3"] <- 1
      ints3[ints3$VarGleShrt=="GDM", "DepvOrderinST3"] <- 5

    
    # Part 3
      
      tbls3o2_stbib2 <- merge(tbls3o2_stbib2, ints3, all=T) 
      tbls3o2_stbib2 <- tbls3o2_stbib2[order(tbls3o2_stbib2$"DepvOrderinST3"), ]   
      
      t2p1 <- tbls3o2_stbib2[tbls3o2_stbib2$DepVarGle %in% c("AUC glucose (log transfromed, winsorized and standardized)", "OGTT - 2h Glucose", "Fasting Glucose"), c("AssocTableDepVDescCol", "IndepVarSub", "Beta [CI]_st", "P_st", "Beta [CI]_bib", "P_bib", "Beta [CI]_MetaA", "MetaA_P", "I2_MetaA", "QEp_MetaA")]
      
      t2p2<- tbls3o2_stbib2[tbls3o2_stbib2$DepVarGle %in% c("GDM (IADPSG)", "GDM (South Asian cutoffs)"),c("AssocTableDepVDescCol","IndepVarSub", "OR [95%CI]_st", "P_st", "OR [95%CI]_bib", "P_bib", "OR [95%CI]_MetaA", "MetaA_P", "I2_MetaA", "QEp_MetaA")]
      
      colnames(t2p1)[colnames(t2p1)=="Beta [CI]_st"] <- "Beta/OR [95%CI]_st"
      colnames(t2p1)[colnames(t2p1)=="Beta [CI]_bib"] <- "Beta/OR [95%CI]_bib"
      colnames(t2p1)[colnames(t2p1)=="Beta [CI]_MetaA"] <- "Beta/OR [95%CI]_MetaA"
      
      colnames(t2p2)[colnames(t2p2)=="OR [95%CI]_st"] <- "Beta/OR [95%CI]_st"
      colnames(t2p2)[colnames(t2p2)=="OR [95%CI]_bib"] <- "Beta/OR [95%CI]_bib"
      colnames(t2p2)[colnames(t2p2)=="OR [95%CI]_MetaA"] <- "Beta/OR [95%CI]_MetaA"
      
      tbls3o2_stbib3 <- rbind(t2p1, t2p2)
      
      tbls3o2_stbib3$IndepVarSub <- c("PRS (T1 vs. T2)", "PRS (T1 vs. T3)")
      tbls3o2_stbib3[tbls3o2_stbib3$IndepVarSub=="T2dMah14PrsTer2", "IndepVarSub"] <- "PRS (T2 vs. T1)"
      tbls3o2_stbib3[tbls3o2_stbib3$IndepVarSub=="T2dMah14PrsTer3", "IndepVarSub"] <- "PRS (T3 vs. T1)"
    

    
  # Save output

      write.csv(tbls3o2_stbib3, "SuplementaryFile1d.csv", row.names=F, quote=F)








   



# Suplementaruy file 1e (Table S5): GDM(South Asian definition)  ~ Interactions PRS * RFs with adjustments  
  # Remark: Steps are similar to Table 4
  # Prep tables needed for the analysis
    
    tbl3_depvs_d <- cvt[cvt$VarGleShrt %in% c("GDM"), ]
    tbl3_indepvs_d <- cvt[cvt$VarGleShrt =="T2dMah14PrsZ",  ]
    tbl3_intervs_d <- cvt[cvt$VarGleShrt %in% c("Age","BMI_v2", "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "lowDQS"),  ]
    
    tbl3_adjs_l <- list()
    tbl3_adjs_l[[1]] <- cvt[cvt$VarGleShrt %in% c("Age","BMI_v2", paste0("PC", 1:5), "Education2_n","BornInSaTF",   "Parity", "ParHistDiabYN", "lowDQS"), ]
    names(tbl3_adjs_l)<-paste0("model", 1:length(tbl3_adjs_l))
    tbl3pref <- "Table4_AssocInteractwAdj"

  #  Run association tests 

    tbl3o1l <- PGPmisc:::AssocInterMultiStudies(data=sbl, depVs_d=tbl3_depvs_d, indepVs_d=tbl3_indepvs_d, adjVs_l=tbl3_adjs_l,depVs_GleName_Coli=1, indepVs_GleName_Coli=1,adjVs_GleName_Coli=1, depVs_ColNames_Coli=2:3, indepVs_ColNames_Coli=2:3, depVs_ModelTypes_Coli=7, adjVs_ColNames_Coli=2:3, interVs_d=tbl3_intervs_d, interVs_GleName_Coli=1, interVs_ColNames_Coli=2:3,   saveShortOutput=T, saveLongOutput=T, OutputPath=o2_p, OutFileNamePrefx=tbl3pref, OutFileNameSufx=paste0("_", lmd))[[1]]



  # Horizontally merge START and BiB association test results

    tbl3o1l_st <- tbl3o1l[tbl3o1l$Study== "START", ]
    tbl3o1l_bib <- tbl3o1l[tbl3o1l$Study== "BiB", ]
    tbl3colsToKeep<- c("DepVarGle", "IndepVarIGle", "InterVarIGle", "IndepVarSub_gleVars")
    
    colnames(tbl3o1l_st)[which(!colnames(tbl3o1l_st) %in% tbl3colsToKeep)] <- paste0(colnames(tbl3o1l_st)[which(!colnames    (tbl3o1l_st) %in% tbl3colsToKeep)], "_st")
    
    colnames(tbl3o1l_bib)[which(!colnames(tbl3o1l_bib) %in% tbl3colsToKeep)] <- paste0(colnames(tbl3o1l_bib)[which  (!colnames  (tbl3o1l_bib) %in% tbl3colsToKeep)], "_bib")
  
    tbl3o2_stbib <- merge(tbl3o1l_st, tbl3o1l_bib, tbl3colsToKeep, suf=c("_st","_bib"), all=T)
    
    # Visually check the merging 
    checkmt<-unique(tbl3o2_stbib[,c("DepVarGle","IndepVarIGle","InterVarIGle","IndepVarSub_gleVars","IndepVarSub_st",  "IndepVarSub_bib")])
    checkmt
    # merging is ok ! 




  # Remove lowDQS rows in tbl3o2_stbib before meta-A 
    # lowDQS is removed at this stage because it is only tested for interaction in START
    tbl3o2_stbib_2 <- tbl3o2_stbib[-which(tbl3o2_stbib$InterVarIGle == "Low diet quality score (vs Medium + High)"),]

  # Run Meta-Analysis (Fixed effects) 
    tbl3o2_l<-list()
    metaAMeth <- "FE"
    for( i in 1:nrow(tbl3o2_stbib_2)){
    tbl3o2_l[[i]] <- list()
	if(!is.na(tbl3o2_stbib_2[i, "Beta_st"]) & !is.na(tbl3o2_stbib_2[i, "Beta_bib"])){

	  esci2 <- escalc(yi=as.numeric(tbl3o2_stbib_2[i,paste0("Beta", c("_st", "_bib"))]), sei=as.numeric(tbl3o2_stbib_2[i,paste0("SE", c("_st", "_bib"))]), measure="SMD")

	  mai2 <- rma(yi, vi, method=metaAMeth, data=esci2)

	  tbl3o2_l[[i]] <- mai2
  
	  tbl3o2_stbib_2[i,"MetaA_Beta"]  <- mai2$"beta"
      tbl3o2_stbib_2[i,"MetaA_SE"]    <- mai2$"se"
  
      tbl3o2_stbib_2[i,"MetaA_Beta95%CiLB"]    <- mai2$"ci.lb"
      tbl3o2_stbib_2[i,"MetaA_Beta95%CiUB"]    <- mai2$"ci.ub"
  
  
      tbl3o2_stbib_2[i,"MetaA_P"]     <- mai2$"pval"
      tbl3o2_stbib_2[i,"MetaA_I2"]    <- mai2$"I2"
      tbl3o2_stbib_2[i,"MetaA_QE"]    <- mai2$"QE"
      tbl3o2_stbib_2[i,"MetaA_QEp"]   <- mai2$"QEp"
      tbl3o2_stbib_2[i,"MetaA_Model"] <- metaAMeth
      if(!is.na(tbl3o2_stbib_2[i, "OR_st"])){
        tbl3o2_stbib_2[i,"MetaA_OR"]         <- exp(mai2$"beta")
	    tbl3o2_stbib_2[i,"MetaA_95%CiLB"]         <- exp(mai2$"ci.lb")
	    tbl3o2_stbib_2[i,"MetaA_95%CiUB"]         <- exp(mai2$"ci.ub")
      }
    }
    }

  # Add back dietq_lo to the START model 
  
    # Rq: dietq_lo was only adjuted for in START, and hence association results between dietq_lo was not included in the meta-analysis. I want to add it back to the  meta-A table 
    # prep the dietqlo data before I add it back 
    st_dql <- tbl3o2_stbib[which(tbl3o2_stbib$InterVarIGle == "Low diet quality score (vs Medium + High)"),]
  
    tbl3o2_stbib_3 <- merge(tbl3o2_stbib_2,st_dql, all=T)
    head(tbl3o2_stbib_3)


  # Format Output table

    # Part 1 :  
      # add () and [], round values, replace NA[NA-NA] etc ...
      tbl3o2_stbib2 <- tbl3o2_stbib_3[,tbl3colsToKeep]
      tbl3o2_stbib2$"Beta (SE)_st" <- paste0(round(tbl3o2_stbib_3$"Beta_st",3), " (", round(tbl3o2_stbib_3$"SE_st",3), ")")
      
      tbl3o2_stbib2$"Beta [95%CI]_st" <- paste0(round(tbl3o2_stbib_3$"Beta_st",3), " [", round(tbl3o2_stbib_3$"Beta95%CiLB_st",3), " - ", round(tbl3o2_stbib_3$"Beta95%CiUB_st",3), "]")
      
      tbl3o2_stbib2$"OR [95%CI]_st" <- paste0(round(tbl3o2_stbib_3$"OR_st",2), " [", round(tbl3o2_stbib_3$"95%CiLB_st" ,2), "-", round(tbl3o2_stbib_3$"95%CiUB_st",2),"]")
      tbl3o2_stbib2$"P_st" <- PGPmisc:::formatPvals(tbl3o2_stbib_3$"P_st")
      
      tbl3o2_stbib2$"Beta (SE)_bib" <- paste0(round(tbl3o2_stbib_3$"Beta_bib",3), " (", round(tbl3o2_stbib_3$"SE_bib",3), ")")
      
      tbl3o2_stbib2$"Beta [95%CI]_bib" <- paste0(round(tbl3o2_stbib_3$"Beta_bib",3), " [", round(tbl3o2_stbib_3$"Beta95%CiLB_bib",3), " - ", round(tbl3o2_stbib_3$"Beta95%CiUB_bib",3), "]")
      
      tbl3o2_stbib2$"OR [95%CI]_bib" <- paste0(round(tbl3o2_stbib_3$"OR_bib",2), " [", round(tbl3o2_stbib_3$"95%CiLB_bib" ,2), "-", round(tbl3o2_stbib_3$"95%CiUB_bib",2),"]")
      #some rows have NA P values for BiB because dietq_lo is not adjuted for in this study   
      w<-which(is.na(tbl3o2_stbib_3$"P_bib"))
      tbl3o2_stbib2[-w, "P_bib"] <- PGPmisc:::formatPvals(tbl3o2_stbib_3[-w,"P_bib"])
      
      
      tbl3o2_stbib2[-w,"Beta (SE)_MetaA"] <- paste0(round(tbl3o2_stbib_3[-w,"MetaA_Beta"],3), " (", round(tbl3o2_stbib_3[-w,"MetaA_SE"],3), ")")
      
      tbl3o2_stbib2[-w,"Beta [95%CI]_MetaA"] <- paste0(round(tbl3o2_stbib_3[-w,"MetaA_Beta"],3), " [", round(tbl3o2_stbib_3[-w,"MetaA_Beta95%CiLB"],3), " - ", round(tbl3o2_stbib_3[-w,"MetaA_Beta95%CiUB"],3) ,"]")
      
      tbl3o2_stbib2[-w,"OR [95%CI]_MetaA"] <- paste0(round(tbl3o2_stbib_3[-w,"MetaA_OR"],2), " [", round(tbl3o2_stbib_3[-w,"MetaA_95%CiLB"],2), "-", round(tbl3o2_stbib_3[-w,"MetaA_95%CiUB"],2),"]")
      tbl3o2_stbib2[-w,"MetaA_P"] <- PGPmisc:::formatPvals(tbl3o2_stbib_3[-w,"MetaA_P"])
      tbl3o2_stbib2[-w,"I2_MetaA"] <- round(tbl3o2_stbib_3[-w,"MetaA_I2"],0)
      tbl3o2_stbib2[-w,"QEp_MetaA"] <- PGPmisc:::formatPvals(tbl3o2_stbib_3[-w,"MetaA_QEp"])
      
      tbl3o2_stbib2[which(tbl3o2_stbib2$"OR [95%CI]_MetaA"=="NA [NA-NA]"),"OR [95%CI]_MetaA"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"OR [95%CI]_st"=="NA [NA-NA]"),"OR [95%CI]_st"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"OR [95%CI]_bib"=="NA [NA-NA]"),"OR [95%CI]_bib"]<-NA
      
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta (SE)_bib"=="NA (NA)"),"Beta (SE)_bib"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta (SE)_st"=="NA (NA)"),"Beta (SE)_st"]<-NA
      
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta [95%CI]_bib"=="NA [NA - NA]"),"Beta [95%CI]_bib"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta [95%CI]_st"=="NA [NA - NA]"),"Beta [95%CI]_st"]<-NA
      tbl3o2_stbib2[which(tbl3o2_stbib2$"Beta [95%CI]_bib"=="NA [NA - NA]"),"Beta [95%CI]_bib"]<-NA
  

  
    # Part 2
      int3<- cvt[cvt$VarGleShrt %in% c("GDM_iadpsg", "AUClwz", "Ogtt2hG", "FPG", "GDM"), c("Variable_Gle", "VarName_START", "VarName_BiB", "TableLegend","VarGleShrt")]
      colnames(int3)[colnames(int3) == "Variable_Gle"]<- "DepVarGle"
      
      int3[int3$VarGleShrt=="GDM", 		"AssocTableDepVDescCol"] <- "GDM (South Asian cutoffs)"
      int3[int3$VarGleShrt=="GDM_iadpsg", 		"AssocTableDepVDescCol"] <- "GDM (IADPSG)"
      int3[int3$VarGleShrt=="Ogtt2hG", 		"AssocTableDepVDescCol"] <- "2h post-load glucose"
      int3[int3$VarGleShrt=="FPG", 		"AssocTableDepVDescCol"] <- "Fasting plasma glucose"
      int3[int3$VarGleShrt=="AUClwz", 		"AssocTableDepVDescCol"] <- "AUC glucose"
      
      int3[int3$VarGleShrt=="FPG", 		"DepvOrderinT3"] <- 1
      int3[int3$VarGleShrt=="Ogtt2hG", 	"DepvOrderinT3"] <- 2
      int3[int3$VarGleShrt=="AUClwz",  	"DepvOrderinT3"] <- 3
      int3[int3$VarGleShrt=="GDM_iadpsg",	"DepvOrderinT3"] <- 4
      int3[int3$VarGleShrt=="GDM", 		"DepvOrderinT3"] <- 5
      colnames(int3)[colnames(int3) == "VarGleShrt"]<- "DepVarGleShrt"
      int3 <-  int3[,c("DepVarGle", "DepvOrderinT3", "AssocTableDepVDescCol", "DepVarGleShrt")] 
      
      
      indepint3<- cvt[cvt$VarGleShrt %in% c("T2dMah14PrsZ", "Age","BMI_v2", "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "lowDQS"), c("Variable_Gle", "VarName_START", "VarName_BiB", "TableLegend","VarGleShrt")]
      colnames(indepint3)[colnames(indepint3) == "Variable_Gle"]<- "InterVarIGle"
      
      indepint3[indepint3$VarGleShrt=="Age",           "AssocTableIndepVDescCol"] <- "Age"
      indepint3[indepint3$VarGleShrt=="BMI_v2",        "AssocTableIndepVDescCol"] <- "BMI"
      indepint3[indepint3$VarGleShrt=="Education2_n",  "AssocTableIndepVDescCol"] <- "Education level"
      indepint3[indepint3$VarGleShrt=="BornInSaTF",      "AssocTableIndepVDescCol"] <- "Born in South Asia (Yes/No)"
      indepint3[indepint3$VarGleShrt=="Parity",        "AssocTableIndepVDescCol"] <- "Parity"
      indepint3[indepint3$VarGleShrt=="ParHistDiabYN", "AssocTableIndepVDescCol"] <- "Parental history of T2D (Yes/No)"
      indepint3[indepint3$VarGleShrt=="lowDQS",        "AssocTableIndepVDescCol"] <- "Low diet quality (Yes/No)"
      indepint3[indepint3$VarGleShrt=="T2dMah14PrsZ",   "AssocTableIndepVDescCol"] <- "PRS"
      
      # Set order of covars in table 
      indepint3[indepint3$VarGleShrt=="T2dMah14PrsZ",  "IndepVvOrdinT3"] <- 1
      indepint3[indepint3$VarGleShrt=="Age",           "IndepVvOrdinT3"] <- 2
      indepint3[indepint3$VarGleShrt=="BMI_v2",        "IndepVvOrdinT3"] <- 3
      indepint3[indepint3$VarGleShrt=="Education2_n",  "IndepVvOrdinT3"] <- 7
      indepint3[indepint3$VarGleShrt=="BornInSaTF",    "IndepVvOrdinT3"] <- 4
      indepint3[indepint3$VarGleShrt=="Parity",        "IndepVvOrdinT3"] <- 6
      indepint3[indepint3$VarGleShrt=="ParHistDiabYN", "IndepVvOrdinT3"] <- 5
      indepint3[indepint3$VarGleShrt=="lowDQS",        "IndepVvOrdinT3"] <- 8
      
      colnames(indepint3)[colnames(indepint3) == "VarGleShrt"]<- "InterVarGleShrt"
      
      indepint3 <-  indepint3[,c("InterVarIGle","AssocTableIndepVDescCol", "IndepVvOrdinT3", "InterVarGleShrt")] 

  
    # Part 3 
      names(table(indepint3$InterVarIGle))
      names(table(tbl3o2_stbib2$InterVarIGle))
      
      names(table(int3$"DepVarGle"))
      names(table(tbl3o2_stbib2$"DepVarGle"))
      
      colnames(tbl3o2_stbib2) [colnames(tbl3o2_stbib2) %in%  colnames(int3)]
      colnames(tbl3o2_stbib2) [colnames(tbl3o2_stbib2) %in%  colnames(indepint3)]
      
      tbl3o2_stbib2 <- merge(tbl3o2_stbib2, int3, all=T) 
      tbl3o2_stbib2 <- merge(tbl3o2_stbib2, indepint3, all=T)
      tbl3o2_stbib2 <- tbl3o2_stbib2[order(tbl3o2_stbib2$"DepvOrderinT3", tbl3o2_stbib2$"IndepVvOrdinT3"), ]   
      tbl3o2_stbib2$interactionTerm <- paste0( "PRS x ", tbl3o2_stbib2$AssocTableIndepVDescCol)
      
      t3p1<- tbl3o2_stbib2[tbl3o2_stbib2$"DepVarGleShrt" %in% c("AUClwz", "Ogtt2hG", "FPG"), c("DepVarGle", "interactionTerm", "Beta [95%CI]_st", "P_st", "Beta [95%CI]_bib", "P_bib", "Beta [95%CI]_MetaA", "MetaA_P", "I2_MetaA", "QEp_MetaA")]
      
      t3p2<- tbl3o2_stbib2[tbl3o2_stbib2$"DepVarGleShrt" %in% c("GDM_iadpsg", "GDM"),c("DepVarGle", "interactionTerm", "OR [95%CI]_st", "P_st", "OR [95%CI]_bib", "P_bib", "OR [95%CI]_MetaA", "MetaA_P", "I2_MetaA", "QEp_MetaA")]
      
      colnames(t3p1)[which(colnames(t3p1) %in% c("Beta [95%CI]_st", "Beta [95%CI]_bib", "Beta [95%CI]_MetaA"))]  <- c("Beta/OR [95%CI] _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
      
      colnames(t3p2)[which(colnames(t3p2) %in% c("OR [95%CI]_st","OR [95%CI]_bib","OR [95%CI]_MetaA"))] <- c("Beta/OR [95%CI] _st", "Beta/OR [95%CI] _bib", "Beta/OR [95%CI] _MetaA")
      
      tbl3o2_stbib3 <- rbind(t3p1, t3p2)
      
      tbl3o2_stbib3[which(tbl3o2_stbib3$"Beta/OR [95%CI] _bib"=="NA [NA - NA]"),"Beta/OR [95%CI] _bib"]<-NA


  # Save outputs
    write.table(tbl3o2_stbib3, "Table4.tsv", row.names=F, quote=F, sep="\t")




# Suplementaruy file 1f (Table S6): Subgroup analysis for the interactions PRS x BMI and PRS x lowDQ
  # Remark: the table was manually filled in (in excel) based on the following results 
  # GDM (IADPSG) ~ PRS By BMI categories 
    # in  BiB
      f1<-"drvgesdiab01 ~ T2dMah14PrsZ + agemy_mbqall + BornInSA_AL + edu0mumede4c1234 + eclregpart + PC1 + PC2 + PC3 + PC4 + PC5"
      
      d1 <- bibinterd[bibinterd$BMI_Ter_v3==1,]
      m<-glm(f1, data=d1, family=binomial)
      sm <- summary(m)
      ttt<-as.data.frame(t(cbind(round(exp(sm$coef),2)[,1],sm$coef[,4], round(exp(confint(glm(f1, data=d1, family=binomial))),2))[2,]))
      colnames(ttt)[1:2] <- c("OR","P")
      ttt$N<-length(m$residuals)
      
      
      d2 <- bibinterd[bibinterd$BMI_Ter_v3==2,]
      m<-glm(f1, data=d2, family=binomial)
      sm <- summary(m)
      tttb<-as.data.frame(t(cbind(round(exp(sm$coef),2)[,1],sm$coef[,4], round(exp(confint(glm(f1, data=d2, family=binomial))),2))[2,]))
      colnames(tttb)[1:2] <- c("OR","P")
      tttb$N<-length(m$residuals)
      ttt<- rbind(ttt, tttb)
      
      d3 <- bibinterd[bibinterd$BMI_Ter_v3==3,]
      m<-glm(f1, data=d3, family=binomial)
      sm <- summary(m)
      tttb<-as.data.frame(t(cbind(round(exp(sm$coef),2)[,1],sm$coef[,4], round(exp(confint(glm(f1, data=d3, family=binomial))),2))[2,]))
      colnames(tttb)[1:2] <- c("OR","P")
      tttb$N<-length(m$residuals)
      ttt<- rbind(ttt, tttb)

    # in START 
  
      bibintervars <- cvt[cvt[,"VarGleShrt"] %in% c("T2dMah14PrsZ","Age","BMI_v2","BMI_Ter_v3", "T2dMah14PrsTer" ,paste0("PC", 1:5), "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "GDM_iadpsg"),"VarName_BiB"]
      
      stintervars <- cvt[cvt[,"VarGleShrt"] %in% c("T2dMah14PrsZ","Age","BMI_v2","BMI_Ter_v3", "T2dMah14PrsTer", paste0("PC", 1:5), "Education2_n","BornInSaTF", "Parity","ParHistDiabYN", "GDM_iadpsg"),"VarName_START" ]
      
      bibinterd <- sbl[[2]][ , bibintervars] 
      stinterd <- sbl[[1]][ , stintervars] 
      
      bibinterd <- bibinterd[complete.cases(bibinterd),]
      stinterd <- stinterd[complete.cases(stinterd),]
      
      
      f1<-"anygdm ~ T2dMah14PrsZ + agemom + BornInSAtf + mbeyself + parity + PC1_v1.0 + PC2_v1.0 + PC3_v1.0 + PC4_v1.0 + PC5_v1.0"
      
      d1s <- stinterd[stinterd$BMI_Ter_v3==1,]
      m<-glm(f1, data=d1s, family=binomial)
      
      sm <- summary(m)
      ttts<-as.data.frame(t(cbind(round(exp(sm$coef),2)[,1],sm$coef[,4], round(exp(confint(glm(f1, data=d1s, family=binomial))),2))[2,]))
      colnames(ttts)[1:2] <- c("OR","P")
      ttts$N<-length(m$residuals)
      
      
      d2s <- stinterd[stinterd$BMI_Ter_v3==2,]
      m<-glm(f1, data=d2s, family=binomial)
      sm <- summary(m)
      tttbs<-as.data.frame(t(cbind(round(exp(sm$coef),2)[,1],sm$coef[,4], round(exp(confint(glm(f1, data=d2s, family=binomial))),2))[2,]))
      colnames(tttbs)[1:2] <- c("OR","P")
      tttbs$N<-length(m$residuals)
      ttts<- rbind(ttts, tttbs)
      
      d3s <- stinterd[stinterd$BMI_Ter_v3==3,]
      m<-glm(f1, data=d3s, family=binomial)
      sm <- summary(m)
      tttbs<-as.data.frame(t(cbind(round(exp(sm$coef),2)[,1],sm$coef[,4], round(exp(confint(glm(f1, data=d3s, family=binomial))),2))[2,]))
      colnames(tttbs)[1:2] <- c("OR","P")
      tttbs$N<-length(m$residuals)
      ttts<- rbind(ttts, tttbs)
      
      colnames(ttt) <-  paste0(colnames(ttt), "_bib")
      colnames(ttts) <-  paste0(colnames(ttt), "_st")
      tttf <- cbind(ttt, ttts)



  # FPG ~ PRS By LowDietQuality categories Subgroup analysis 
    # ran in START only since the variable dietq_lo is not available in BiB

    fst <- "lbrfpgv ~ T2dMah14PrsZ + agemom + prepregbmi + mbeyself + parity + PC1_v1.0 + PC2_v1.0 + PC3_v1.0 + PC4_v1.0 + PC5_v1.0 + BornInSAtf"
     
    mst<- lm(fst, data=stinterd)
    d1 <- stinterd[stinterd$dietq_lo==0,]
    summary(lm(fst, data=d1))
    length(lm(fst, data=d1)$residuals)
    round(confint(lm(fst, data=d1)),2)
    
    
    
    d2 <- stinterd[stinterd$dietq_lo==1,]
    summary(lm(fst, data=d2))
    length(lm(fst, data=d2)$residuals)
    round(confint(lm(fst, data=d2)),2)

  # FPG ~ PRS by Born in South Asia categories
    # BiB 


      fbib <- "gttfastglu ~ T2dMah14PrsZ + agemy_mbqall + mms0mbkbmi + edu0mumede4c1234 + eclregpart + PC1 + PC2 + PC3 + PC4 + PC5 "
      
      summary(lm(fbib, data=bibinterd[!as.logical(bibinterd$BornInSA_AL),]))
      round(confint (lm(fbib, data=bibinterd[!as.logical(bibinterd$BornInSA_AL),])),2)
      
      
      summary(lm(fbib, data=bibinterd[as.logical(bibinterd$BornInSA_AL),]))
      round(confint(lm(fbib, data=bibinterd[as.logical(bibinterd$BornInSA_AL),])),2)
      

    # START 
      fst <- "lbrfpgv ~ T2dMah14PrsZ + agemom + prepregbmi + mbeyself + parity + PC1_v1.0 + PC2_v1.0 + PC3_v1.0 + PC4_v1.0 + PC5_v1.0 + dietq_lo"
      
      summary(lm(fst, data=stinterd[stinterd$BornInSAtf==0,]))
      round(confint(lm(fst, data=stinterd[stinterd$BornInSAtf==0,])),2)
      
      
      summary(lm(fst, data=stinterd[stinterd$BornInSAtf==1,]))
      round(confint(lm(fst, data=stinterd[stinterd$BornInSAtf==1,])),2)


  # FPG ~ PRS by BMI categories 

    # BiB 

      f1<-"gttfastglu ~ T2dMah14PrsZ + agemy_mbqall + BornInSA_AL + edu0mumede4c1234 + eclregpart + PC1 + PC2 + PC3 + PC4 + PC5"
      
      d1 <- bibinterd[bibinterd$BMI_Ter_v3==1,]
      m<-lm(f1, data=d1)
      sm <- summary(m)
      round(sm$coef,2)
      length(lm(f1, data=d1)$residuals)
      round(confint(lm(f1, data=d1)),2)
      
      d2 <- bibinterd[bibinterd$BMI_Ter_v3==2,]
      summary(lm(f1, data=d2))
      length(lm(f1, data=d2)$residuals)
      round(confint(lm(f1, data=d2)),2)
      
      
      d3 <- bibinterd[bibinterd$BMI_Ter_v3==3,]
      summary(lm(f1, data=d3))
      length(lm(f1, data=d3)$residuals)
      round(confint(lm(f1, data=d3)),2)
      
    # START 
      
      f1<-"lbrfpgv ~ T2dMah14PrsZ + agemom + BornInSAtf + mbeyself + parity + PC1_v1.0 + PC2_v1.0 + PC3_v1.0 + PC4_v1.0 + PC5_v1.0 + dietq_lo"
      
      d1 <- stinterd[stinterd$BMI_Ter_v3==1,]
      summary(lm(f1 , data=d1))
      length(lm(f1, data=d1)$residuals)
      round(confint(lm(f1, data=d1)),2)
      
      d2 <- stinterd[stinterd$BMI_Ter_v3==2,]
      summary(lm(f1 , data=d2))
      length(lm(f1, data=d2)$residuals)
      round(confint(lm(f1, data=d2)),2)
      
      d3 <- stinterd[stinterd$BMI_Ter_v3==3,]
      summary(lm(f1 , data=d3))
      length(lm(f1, data=d3)$residuals)
      round(confint(lm(f1, data=d3)),2)








  
# In text information  
  # N and % 
    # Here I caluclate some of the numbers included in the main text 
    
    # % GDM (IADPSG)
    round(prop.table(table(sbl[[1]]$anygdm))*100,1)
    round(prop.table(table(sbl[[2]]$drvgesdiab01))*100,1)
    
    # BIB South asians specific definition  
    round(prop.table(table(sbl[[1]]$GdmSaBiBp01))*100,1)
    round(prop.table(table(sbl[[2]]$GDMsaBiBwBfnT2D))*100,1)
    
    # Porportion of south asians 
    round(prop.table(table(sbl[[1]]$CtryOrgin3c))*100,1)
    round(prop.table(table(sbl[[2]]$CountryOrigin3c))*100,1)
       
    # Born in South Asia  
    round(prop.table(table(sbl[[1]]$BornInSAtf ))*100,1)
    round(prop.table(table(sbl[[2]]$BornInSA_AL ))*100,1)
    
    # Year spent in recruitment country 
    round(mean(sbl[[1]][sbl[[1]]$BornInSAtf, "YrsInCA"], na.rm=T),1)
    round(mean(sbl[[2]][sbl[[2]]$BornInSA_AL, "YrsInUK"], na.rm=T),1)
    
    # N pregnancies 
    #primiparous 
    round(prop.table(table(sbl[[1]]$parity3c0v1v2p))[1]*100,1)
    round(prop.table(table(sbl[[2]]$eclregpart3c0v1v2p))[1]*100,1)
    #1 previous pregnancies 
    round(prop.table(table(sbl[[1]]$parity3c0v1v2p))[2]*100,1)
    round(prop.table(table(sbl[[2]]$eclregpart3c0v1v2p))[2]*100,1)
    #2 or more  
    round(prop.table(table(sbl[[1]]$parity3c0v1v2p))[3]*100,1)
    round(prop.table(table(sbl[[2]]$eclregpart3c0v1v2p))[3]*100,1)
  
    # Vegetarians 
    round(prop.table(table(sbl[[1]]$Vegetarian ))*100,1)
    round(prop.table(table(sbl[[2]]$Vegetarian ))*100,1)
    
    round(prop.table(table(sbl[[1]]$mbeyself2cse2v3p ))*100,1)
    round(prop.table(table(sbl[[2]]$Vegetarian ))*100,1) 
  
    # PRS range (standardized PRS)
    round(range(st1$"T2dMah14PrsZ", na.rm=T),2)
    round(range(bib1$"T2dMah14PrsZ", na.rm=T),2)
    
    # Education (2 categories variable)
    round(prop.table(table(sbl[[1]]$mbeyself2cse2v3p ))*100,2)
    round(prop.table(table(sbl[[2]]$edu0mumede2cse3v4 ))*100,2)



  # Power analysis of the interactions 
    # Keep only interaction results of the rows of interest 
      interpower <- tbl3o2_stbib_3[which((tbl3o2_stbib_3$"P_st" < 0.05 | tbl3o2_stbib_3$"P_bib" < 0.05) ),] 
    
      interpower[interpower$P_st <0.05, "CorX1X2Y"] <- interpower[interpower$P_st <0.05, "Beta_st"]
      interpower[!is.na(interpower$P_bib) & interpower$P_bib <0.05, "CorX1X2Y"] <- interpower[!is.na(interpower$P_bib) & interpower$P_bib <0.05, "Beta_bib"]
    
    #visually check phenotypes 
      interpower[, c("DepVarGle", "IndepVarIGle", "InterVarIGle")]
      interpower$interactionTermE <- c("AUC glucose", rep("Fasting glucose", 3), "GDM (IADPSG criteria)")
      interpower$SourceStudyE <- c(rep("START", 4),"BiB")
      interpower$TargetStudyE <- c(rep("BiB", 4),"START")
    
    #fill in Beta interactions
    interpower$InteractionCorE <-NA
    interpower$InteractionCorE[1:4] <- interpower$Beta_st[1:4]
    interpower$InteractionCorE[5] <- interpower$Beta_bib[5]
  
    #Visually check formulas 
    interpower$Formula_st 
    interpower$Formula_bib
    
    
    interpower$CorX1Y <- NA
    interpower$CorX2Y <- NA
    interpower$CorX1Y <- NA
    interpower$N <- NA
  
    set.seed(13)
    # 1 - FPG ~ PRS * dietqlo
    #2 are continuous (y and x1) vs 1 binary (x2) 

  
    interpower$CorX1Y[4] <- cor(st1$lbrfpgv ,st1$T2dMah14PrsZ, use="complete.obs")
    interpower$CorX2Y[4] <- cor(st1$lbrfpgv ,st1$dietq_lo, use="complete.obs")
    interpower$CorX1X2[4] <- cor(st1$T2dMah14PrsZ ,st1$dietq_lo, use="complete.obs")
    interpower$N[4] <- sum(complete.cases(st1[,c("lbrfpgv", "T2dMah14PrsZ", "dietq_lo" )]))
        
    # Smallest detectable effect 
      power_PRSdqlo <- power_interaction(
      alpha = 0.05, # p-value
      N = interpower$N[4] , # sample size
      r.x1.y = interpower$CorX1Y[4], # correlation between x1 and y
      r.x2.y = interpower$CorX2Y[4], # correlation between x2 and y
      r.x1.x2 = interpower$CorX1X2[4], # correlation between x1 and x2
      r.x1x2.y = seq(.1,.2,by=.005), # correlation between x1x2 and y
      k.x2 = 2,
      cl = 35, # number of cores
      n.iter = 10000,
      adjust.correlations = FALSE
      ) 
      #
      power_PRSdqlo
    
    # powet to detect our effect 
      power_PRSdqlo1 <- power_interaction(
        alpha = 0.05, # p-value
        N = interpower$N[4] , # sample size
        r.x1.y = interpower$CorX1Y[4], # correlation between x1 and y
        r.x2.y = interpower$CorX2Y[4], # correlation between x2 and y
        r.x1.x2 = interpower$CorX1X2[4], # correlation between x1 and x2
        r.x1x2.y = interpower$CorX1X2Y[4], # correlation between x1x2 and y
        k.x2 = 2,
        cl = 35, # number of cores
        n.iter = 10000,
        adjust.correlations = FALSE
        ) 
      
      interpower$Power[4] <- power_PRSdqlo1[2]  									
  
    # 2 - GDM ~ PRS * BMI :  
    # 3 - Auc ~ PRS * BMI (Data not Shown)
    # all three variables are continuous in both studies 

    interpower$CorX1Y[1] <- cor(st1$AUClwz  ,st1$T2dMah14PrsZ, use="complete.obs")
    interpower$CorX2Y[1] <- cor(st1$AUClwz  ,st1$prepregbmi, use="complete.obs")
    interpower$CorX1X2[1] <- cor(st1$T2dMah14PrsZ  ,st1$prepregbmi, use="complete.obs")
    interpower$N[1] <- sum(complete.cases(bib1[,c("AUClwz", "T2dMah14PrsZ", "mms0mbkbmi"  )]))
    power1<- power_interaction(n.iter = 10000,
    						   alpha = 0.05,
    						   interpower$N[1],
    						   r.x1.y = interpower$CorX1Y[1],
    						   r.x2.y =  interpower$CorX2Y[1], # correlation between x2 and y
    						   r.x1.x2 = interpower$CorX1X2[1],# correlation between x1 and x2
    						   r.x1x2.y = interpower$CorX1X2Y[1],  # correlation between x1x2 and y
    						   cl = 35# number of clusters for parallel analyses
    						   )
    interpower$Power[1] <- power1[2]

  
    # 4 - FPG ~ PRS * BMI (Data not Shown)
    #all three variables are continuous in both studies 
    interpower$CorX1Y[2] <- cor(st1$lbrfpgv ,st1$T2dMah14PrsZ, use="complete.obs")
    interpower$CorX2Y[2] <- cor(st1$lbrfpgv ,st1$prepregbmi, use="complete.obs")
    interpower$CorX1X2[2] <- cor(st1$T2dMah14PrsZ ,st1$prepregbmi, use="complete.obs")
    interpower$N[2] <- sum(complete.cases(bib1[,c("gttfastglu", "T2dMah14PrsZ", "mms0mbkbmi" )]))
    						
    power2<- power_interaction(n.iter = 10000,
    										alpha = 0.05,
    										interpower$N[2],
    										r.x1.y = interpower$CorX1Y[2],
    										r.x2.y =  interpower$CorX2Y[2],
    										r.x1.x2 = interpower$CorX1X2[2], # correlation between x1 and x2
    										r.x1x2.y = interpower$CorX1X2Y[2], # correlation between x1x2 and y
    										cl = 35# number of clusters for parallel analyses
    										)	
    interpower$Power[2] <- power2[2]

  
    # 5 - FPG ~ PRS * Born in South Asia
    # 2 are continuous vs 1 binary

    interpower$CorX1Y[3] <- cor(st1$lbrfpgv ,st1$T2dMah14PrsZ, use="complete.obs")
    interpower$CorX2Y[3] <- cor(st1$lbrfpgv ,st1$BornInSAtf, use="complete.obs")
    interpower$CorX1X2[3] <- cor(st1$T2dMah14PrsZ ,st1$BornInSAtf, use="complete.obs")
    interpower$N[3] <- sum(complete.cases(bib1[,c("gttfastglu", "T2dMah14PrsZ", "BornInSA_AL" )]))
    				
    
    power3 <- power_interaction(n.iter = 10000,
    							alpha = 0.05,
    							interpower$N[3],
    							r.x1.y = interpower$CorX1Y[3],
    							r.x2.y =  interpower$CorX2Y[3],
    							r.x1.x2 = interpower$CorX1X2[3],
    							r.x1x2.y = interpower$CorX1X2Y[3],
    							k.x2 = 2, # Born in south Asia is Binary
    							cl = 35) # number of clusters for parallel analyses
    										
    interpower$Power[3] <- power3[2]

  

    # 2 are continuous (x1 and x2) vs 1 binary (y) 

    interpower$CorX1Y[5] <- cor(bib1$drvgesdiab01 ,bib1$T2dMah14PrsZ, use="complete.obs")
    interpower$CorX2Y[5] <- cor(bib1$drvgesdiab01 ,bib1$mms0mbkbmi, use="complete.obs")
    interpower$CorX1X2[5] <- cor(bib1$T2dMah14PrsZ ,bib1$mms0mbkbmi, use="complete.obs")
    interpower$N[5] <- sum(complete.cases(st1[,c("anygdm", "T2dMah14PrsZ", "prepregbmi")]))
    power5b <- power_interaction(n.iter = 10000,
      alpha = 0.05,
      interpower$N[5],
      r.x1.y = interpower$CorX1Y[5],
      r.x2.y =  interpower$CorX2Y[5],
      r.x1.x2 = interpower$CorX1X2[5], 
      r.x1x2.y = interpower$CorX1X2Y[5], 
      k.y = 2,
      cl = 35 )
    
    interpower$Power[5] <- power5[2]



