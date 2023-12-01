library("dplyr")
library("readxl")
source("../R/01_clean_lipid_list.R")
source("../R/02_lipid_miner.R")


# load the lipid list from MS-DIAL output
ll<-lipid.list.from.msdial("../../AgLung_hIPF_lipids_for_stats.xlsx")
print(ll)

# use clean lipid list function then run the miner
ll.clean<-clean.lipid.list(ll)
print(ll.clean)
result <- lipid.miner(ll.clean)
