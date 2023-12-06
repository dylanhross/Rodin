### integration test for the lipid.list.from.msdial function
### make sure the rest of the analysis can work from this
### modified/reformatted input


suppressMessages(library("dplyr"))
suppressMessages(library("readxl"))
suppressMessages(library("devtools"))
suppressMessages(devtools::load_all("../"))


# load the lipid list from MS-DIAL output
ll<-lipid.list.from.msdial("../../AgLung_hIPF_lipids_for_stats.xlsx")
print("---------- load lipid list from xlsx ----------")
print(ll)

# use clean lipid list function then run the miner
ll.clean<-clean.lipid.list(ll)
print("---------- clean the lipid list ----------")
print(ll.clean)
print("---------- run the lipid miner ----------")
lipid.miner(ll.clean, name = "test")
# results should now be in test.intact, test.chain, and test.allchain
print(test.intact)


### the rest are parts adapted from the example script

# try to run one of the tests querying test against the UniverseExample
Universe <- UniverseExample
cleaned.universeExample <- clean.lipid.list(Universe)
lipid.miner(cleaned.universeExample, name = "Universe")

# run the Fisher (and create results objects)
intact_fisher_cat_result <- intact.fisher(test.intact$Category, Universe.intact$Category)
intact_fisher_main_result <- intact.fisher(test.intact$`Main class`, Universe.intact$`Main class`)
intact_fisher_sub_result <- intact.fisher(test.intact$`Sub class`, Universe.intact$`Sub class`)
chain_fisher_result <- chain.fisher(test.chain, Universe.chain)
allchains_fisher_result <- allchains.fisher(test.allchains, Universe.allchains)

#reformat the tables to add the type of test
intact_fisher_cat_result<- cbind(Type=(rep("Intact")),intact_fisher_cat_result)
intact_fisher_main_result<- cbind(Type=(rep("Main Class")),intact_fisher_main_result)
intact_fisher_sub_result<- cbind(Type=(rep("Sub Class")),intact_fisher_sub_result)
chain_fisher_result<- cbind(Type=(rep("Chain characteristics")),chain_fisher_result)
allchains_fisher_result<- cbind(Type=(rep("Specific chain")),allchains_fisher_result)

# some graphics 
#category stack barplot
grid.arrange(intact.cat.stack(Universe.intact)+ggtitle("Universe"),intact.cat.stack(test.intact)+ggtitle("test"),ncol=2)

#Main class stack barplot
grid.arrange(intact.main.stack(Universe.intact)+ggtitle("Universe"),intact.main.stack(test.intact)+ggtitle("test"),ncol=2)

#Sub class stack barplot
grid.arrange(intact.sub.stack(Universe.intact)+ggtitle("Universe"),intact.sub.stack(test.intact)+ggtitle("test"),ncol=2)
