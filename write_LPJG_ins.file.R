# Tom Pugh
# outsourced from main code by Annemarie Eckes-Shephard May 2021
# convenience function that contains all code bits to create .ins files for LPJ-GUESS-HYD
# output_fol : specify output folder location as string
# basePFT: number between 1 and 5, to reflect PFTs:  'TeBE'(1),'TeBS'(2),'IBS'(3),'TrBE'(4),'TrBR'(5)
# traits_LPJG: list object of traits derived from function  lpjg_traits_conv()
# insfile_template: contains all instructions you want for LPJG to run except for the below and
# the PFT headers in lpjg_traits_conv.R
write_LPJG_ins.file <- function(output_fol,basePFT,traits_LPJG,insfile_template="global_cf_base_Tom.ins",lpjg_version=4.1){
  basePFT_names <- c('TeBE','TeBS','IBS','TrBE','TrBR')
  
  if(lpjg_version!=4.1){
    error("this output file creator can only write insfiles for 4.1 versions of LPJ-GUESS. 
          Pleasemodify the headers (i.e. TeBE_header in lpjg_traits_conv.R")
  }
# Set the name for the output file  
if (basePFT==1) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TeBE.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TeBE.csv",sep="")
} else if (basePFT==2) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TeBS.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TeBS.csv",sep="")
} else if (basePFT==3) {
  LPJG_outfile <- paste(output_fol,"LPJG_PFT_insfile_IBS.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_IBS.csv",sep="")
} else if (basePFT==4) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TrBE.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TrBE.csv",sep="")
} else if (basePFT==5) {
  LPJG_outfile <- paste(output_fol,"/LPJG_PFT_insfile_TrBR.ins",sep="")
  LPJG_summaryfile <- paste(output_fol,"/LPJG_PFT_summary_TrBR.csv",sep="")
} else {
  stop("basePFT must be equal to 1, 2, 3, 4 or 5")
}

# Write out to LPJ-GUESS instruction file
PFTfile <- file(LPJG_outfile)
for (nn in 1:length(traits_LPJG$Ks)) {
  if (nn>1) {
    PFTfile <- file(LPJG_outfile,open="append")
  }

  # 'TeBE'(1),'TeBS'(2),'IBS'(3),'TrBE'(4),'TrBR'(5)
  Line1 <- paste("pft \"PFT",nn,"\" (",sep="")
  if (basePFT==1) {
    Line2 <- TeBE_header
    Line10 <- paste("\t kL_max 7.24") # mean value from Daijun's literature screening,
    # generated in extract_Kleaf_instances_for_trait_network.R
    #previous value (temporarily used): 9  from Sack et Holbrook 2006 Fig 2, mean for Temp woody angiosperms
    traits_LPJG$Kleaf <- 7.24
  } else if (basePFT==2) {
    Line2 <- TeBS_header
    Line10 <- paste("\t kL_max 11.6")# mean value from Daijun's literature screening,
    # generated in extract_Kleaf_instances_for_trait_network.R
    # prev. value: 9, Fig 2, mean for Temp woody angiosperms
    # make sure KL_max is correctly assigned a constant value, also in this output:
    traits_LPJG$Kleaf <- 11.6
  } else if (basePFT==3) {
    Line2 <- IBS_header
    Line10 <- paste("\t kL_max 11.6")# mean value from Daijun's literature screening,
    # generated in extract_Kleaf_instances_for_trait_network.R
    #previous value (temporarily used): 9, Fig 2, mean for Temp woody angiosperms
    traits_LPJG$Kleaf <- 11.6
  } else if (basePFT==4) {
    Line2 <- TrBE_header
    Line10 <- paste("\t kL_max 7.24")# mean value from Daijun's literature screening,
    # generated in extract_Kleaf_instances_for_trait_network.R
    # prev. value (temporarily used): 13 .  Fig 2, mean for Trop. woody angiosperms
    traits_LPJG$Kleaf <- 7.24
  } else if (basePFT==5) {
    Line2 <- TrBR_header
    Line10 <- paste("\t kL_max 17.3")# mean value from Daijun's literature screening,
    # generated in extract_Kleaf_instances_for_trait_network.R
    #previous value (temporarily used): Fig 2, mean for Trop. woody angiosperms
    traits_LPJG$Kleaf <- 17.3
  }
  Line3 <- "\t !Hydraulics"
  Line4 <- paste("\t isohydricity ",traits_LPJG$lambda[nn],sep="")
  Line5 <- paste("\t delta_psi_max ",traits_LPJG$DeltaPsiWW[nn],sep="")
  Line6 <- paste("\t cav_slope ",traits_LPJG$polyslope[nn],sep="") #done
  Line7 <- paste("\t psi50_xylem ",traits_LPJG$P50[nn],sep="")
  Line8 <- paste("\t ks_max ",traits_LPJG$Ks[nn],sep="")
  Line9 <- paste("\t kr_max ",11.2e-4,sep="") # LPJ-GUESS default from Hickler et al. (2006)
  #Line10 <- paste("\t kL_max ",traits_LPJG$Kleaf[nn],sep="")
  Line11 <- paste("\t wooddens ",traits_LPJG$WD[nn],sep="") # done
  Line12 <- paste("\t k_latosa ",traits_LPJG$LS[nn],sep="")
  Line13 <- paste("\t sla ",traits_LPJG$SLA[nn],sep="") # done
  Line14 <- paste("\t cton_leaf_min ",traits_LPJG$CtoNmin_LPJG[nn],sep="")
  if (basePFT==1 | basePFT==4) {
    Line15 <- paste("\t leaflong ",traits_LPJG$leaflong[nn],sep="")
    Line16 <- paste("\t turnover_leaf ",1/traits_LPJG$leaflong[nn],sep="")
  } else {
    #Use LPJ-GUESS standard values for deciduous
    Line15 <- paste("\t leaflong 0.5")
    Line16 <- paste("\t turnover_leaf 1.0")
  } 
  
  writeLines(c(Line1,Line2,Line3,Line4,Line5,Line6,Line7,Line8,Line9,Line10,Line11,Line12,Line13,Line14,Line15,Line16,"",")",""),PFTfile)
  close(PFTfile)
}

# Write out to a set of LPJ-GUESS instruction files, 1 per PFT
for (nn in 1:length(traits_LPJG$Ks)) {
  LPJG_outfile_pft <- paste(LPJG_outfile,".PFT",nn,sep="")
  file.copy(insfile_template,LPJG_outfile_pft,overwrite=T)
  PFTfile <- file(LPJG_outfile_pft,open="append")
  
  Line1 <- paste("pft \"PFT",nn,"\" (",sep="")
  if (basePFT==1) {
    Line2 <- TeBE_header
  } else if (basePFT==2) {
    Line2 <- TeBS_header
  } else if (basePFT==3) {
    Line2 <- IBS_header
  } else if (basePFT==4) {
    Line2 <- TrBE_header
  } else if (basePFT==5) {
    Line2 <- TrBR_header
  }
  Line3 <- "\t !Hydraulics"
  Line4 <- paste("\t isohydricity ",traits_LPJG$lambda[nn],sep="")
  Line5 <- paste("\t delta_psi_max ",traits_LPJG$DeltaPsiWW[nn],sep="")
  Line6 <- paste("\t cav_slope ",traits_LPJG$polyslope[nn],sep="")
  Line7 <- paste("\t psi50_xylem ",traits_LPJG$P50[nn],sep="")
  Line8 <- paste("\t ks_max ",traits_LPJG$Ks[nn],sep="")
  Line9 <- paste("\t kr_max ",11.2e-4,sep="") # LPJ-GUESS default from Hickler et al. (2006)
  #Line10 <- paste("\t kL_max ",traits_LPJG$Kleaf[nn],sep="")
  Line11 <- paste("\t wooddens ",traits_LPJG$WD[nn],sep="")
  Line12 <- paste("\t k_latosa ",traits_LPJG$LS[nn],sep="")
  Line13 <- paste("\t sla ",traits_LPJG$SLA[nn],sep="")
  Line14 <- paste("\t cton_leaf_min ",traits_LPJG$CtoNmin_LPJG[nn],sep="")
  if (basePFT==1 | basePFT==4) {
    Line15 <- paste("\t leaflong ",traits_LPJG$leaflong[nn],sep="")
    Line16 <- paste("\t turnover_leaf ",1/traits_LPJG$leaflong[nn],sep="")
  } else {
    #Use LPJ-GUESS standard values for deciduous
    Line15 <- paste("\t leaflong 0.5")
    Line16 <- paste("\t turnover_leaf 1.0")
  } 
  
  writeLines(c(Line1,Line2,Line3,Line4,Line5,Line6,Line7,Line8,Line9,Line10,Line11,Line12,Line13,Line14,Line15,Line16,"",")",""),PFTfile)
  close(PFTfile)
}

# Write out a trait values table by PFT to be used for post-processing of LPJ-GUESS output
write.table(traits_LPJG,file=LPJG_summaryfile,sep=",",row.names = F)

print(paste0('LPJ-GUESS .ins files for PFT ',basePFT_names[basePFT], ' created in folder ', output_fol))
}