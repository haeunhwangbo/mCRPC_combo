library(tidyverse)
setwd("/Users/haeunhwangbo/Dropbox/PalmerLab/collaborations/Stephanie_Huang/mCRPC_combo/")
source(file.path(getwd(), 'src', 'extract_IPD', 'IPD_Functions.R'))
dir = '/Users/haeunhwangbo/Dropbox/PalmerLab/collaborations/Stephanie_Huang/mCRPC_combo/data/raw/trials/'
digitized_file_prefix = "Prostate_Olaparib_Bono2020_rPFS"
atrisk_file_prefix = "Prostate_Olaparib_Bono2020_rPFS_at-risk"

DF <- read.csv(paste0(dir, atrisk_file_prefix, '.csv'))

AR <- data.frame(interval = seq(1, dim(DF)[1]), 
                 trisk = DF$time, nrisk = DF$treatment, TE = DF$treatment[1],
                 Arm = "treatment")

  
TSD <- read.csv(paste0(dir, digitized_file_prefix, ".csv"))
digizeit <- DIGI.CLEANUP(TSD)
  
  
# if survival is in percentage, convert it to 0-1
if (max(digizeit$S) > 2){
  digizeit$S = digizeit$S / 100
}
# set upper bound to 1
digizeit$S[which(digizeit$S > 1)] <- 1

pub.risk <- K.COORDINATES(AR, digizeit)

IPD <- GENERATEINDIVIDUALDATA(tot.events = unique(AR$TE), 
                              arm.id = unique(AR$Arm), 
                              digizeit = digizeit, 
                              pub.risk = pub.risk)
IPD <- data.frame(Time=IPD$Time, Event=IPD$Event, 
                  Arm=IPD$Arm)
nm <- paste(unique(AR$Slide),  paste(unique(AR$Arm), unique(AR$Subpop), sep='_'),sep='_')
write.csv(IPD, paste0(dir, digitized_file_prefix, "_indiv.csv"), row.names = FALSE)