### integration test for the lipid.list.from.msdial function
### make sure the results from lipid miner are as expected


suppressMessages(library("dplyr"))
suppressMessages(library("readxl"))
suppressMessages(library("devtools"))
suppressMessages(devtools::load_all("../"))


# load the lipid list from MS-DIAL output
ll<-lipid.list.from.msdial("../../AgLung_hIPF_lipids_for_stats.xlsx")
print("---------- load lipid list from xlsx ----------")
#print(ll)

# use clean lipid list function then run the miner
ll.clean<-clean.lipid.list(ll)
print("---------- clean the lipid list ----------")
#print(ll.clean)
print("---------- run the lipid miner ----------")
lipid.miner(ll.clean, name = "test")
# results should now be in test.intact, test.chain, and test.allchain
print("intact")
print(head(test.intact))
print("chain")
print(head(test.chain))
print("allchain")
print(head(test.allchains))

# try to run one of the tests querying test against the UniverseExample
Universe <- UniverseExample
cleaned.universeExample <- clean.lipid.list(Universe)
lipid.miner(cleaned.universeExample, name = "Universe")
print("---------- UniverseExamples miner results ----------")
print("intact")
print(head(Universe.intact))
print("chain")
print(head(Universe.chain))
print("allchain")
print(head(Universe.allchains))

print("---------- run stats ----------")
# run the Fisher (and create results objects)
intact_fisher_cat_result <- intact.fisher(test.intact$Category, Universe.intact$Category)
intact_fisher_main_result <- intact.fisher(test.intact$`Main class`, Universe.intact$`Main class`)
intact_fisher_sub_result <- intact.fisher(test.intact$`Sub class`, Universe.intact$`Sub class`)
chain_fisher_result <- chain.fisher(test.chain, Universe.chain)
allchains_fisher_result <- allchains.fisher(test.allchains, Universe.allchains)