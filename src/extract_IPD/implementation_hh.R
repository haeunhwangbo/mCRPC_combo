library(tidyverse)
setwd("/Users/haeunhwangbo/Dropbox/PalmerLab/collaborations/Stephanie_Huang/mCRPC_combo/")
source(file.path(getwd(), 'src', 'extract_IPD', 'IPD_Functions.R'))
dir = '/Users/haeunhwangbo/Dropbox/PalmerLab/collaborations/Stephanie_Huang/mCRPC_combo/data/raw/trials/'
comb_prefix = "Prostate_Talazoparib-Enzalutamide_TALAPARO2_rPFS_HRD"
mono_prefix = "Prostate_Enzalutamide_TALAPARO2_rPFS_HRD"
atrisk_file = "Prostate_Talaparib_TALAPARO2_rPFS_at-risk"


#DF <- read.csv(file.path("at_risk_table.csv"))
#AR <-filter(DF, Slide == 'C', Arm == 'TC3/IC3')
DF <- read.csv(paste0(dir, atrisk_file, '.csv'))

for(arm in c('mono', 'comb')){
  if((arm == "mono") && (mono_prefix != '')){
    AR <- data.frame(interval = seq(1, dim(DF)[1]), 
                     trisk = DF$time, nrisk = DF$control, TE = DF$control[1],
                     Arm = "control")
    fileprefix = mono_prefix
  }
  else{
    AR <- data.frame(interval = seq(1, dim(DF)[1]), 
                     trisk = DF$time, nrisk = DF$treatment, TE = DF$treatment[1],
                     Arm = "treatment")
    fileprefix = comb_prefix
  }
  
  TSD <- read.csv(paste0(dir, fileprefix, ".csv"))
  digizeit <- DIGI.CLEANUP(TSD)
  
  
  # if survival is in percentage, convert it to 0-1
  if (max(digizeit$S) > 2){
    digizeit$S = digizeit$S / 100
  }
  # set upperbound to 1
  digizeit$S[which(digizeit$S > 1)] <- 1
  
  pub.risk <- K.COORDINATES(AR, digizeit)
  
  IPD <- GENERATEINDIVIDUALDATA(tot.events = unique(AR$TE), 
                                arm.id = unique(AR$Arm), 
                                digizeit = digizeit, 
                                pub.risk = pub.risk)
  IPD <- data.frame(Time=IPD$Time, Event=IPD$Event, 
                    Arm=IPD$Arm)
  nm <- paste(unique(AR$Slide),  paste(unique(AR$Arm), unique(AR$Subpop), sep='_'),sep='_')
  write.csv(IPD, paste0(dir, fileprefix, "_indiv.csv"), row.names = FALSE)
}


#IPD$Time <- as.numeric(as.character(IPD$Time))
#IPD$Event <- as.numeric(as.character(IPD$Event))
#survfit(Surv(Time, Event) ~ 1, data=IPD)