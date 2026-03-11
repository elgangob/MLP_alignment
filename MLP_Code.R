#===================================
#Packages
#===================================

library(tidyverse)

library(seqinr)

#==================================================================================================================
#Data Importing 

#================================
#Import FASTA
#================================

seqdump1 <- read.fasta(file = "Documents/Stevenson Project/seqdump1.txt")

seqdump2 <- read.fasta(file = "Documents/Stevenson Project/seqdump2.txt")

seqdump3 <- read.fasta(file = "Documents/Stevenson Project/seqdump3.txt")

seqdump4 <- read.fasta(file = "Documents/Stevenson Project/seqdump4.txt")

prot <- read.fasta(file = "Documents/Stevenson Project/All Mlp proteins.txt")

#================================
#Getting FASTA in DF format
#================================

df_seqdump1 <- as.data.frame(t(data.frame(Annot = getAnnot(seqdump1))))
rownames(df_seqdump1) <- 1:length(df_seqdump1[,1])
colnames(df_seqdump1) <- "Seq"

df_seqdump2 <- as.data.frame(t(data.frame(Annot = getAnnot(seqdump2))))
rownames(df_seqdump2) <- 1:length(df_seqdump2[,1])
colnames(df_seqdump2) <- "Seq"

df_seqdump3 <- as.data.frame(t(data.frame(Annot = getAnnot(seqdump3))))
rownames(df_seqdump3) <- 1:length(df_seqdump3[,1])
colnames(df_seqdump3) <- "Seq"

df_seqdump4 <- as.data.frame(t(data.frame(Annot = getAnnot(seqdump4))))
rownames(df_seqdump4) <- 1:length(df_seqdump4[,1])
colnames(df_seqdump4) <- "Seq"

df_prot <- as.data.frame(t(data.frame(Annot = getAnnot(prot))))
rownames(df_prot) <- 1:length(df_prot[,1])
colnames(df_prot) <- "Seq"

#Small fix for prot seq
df_prot[166:167,] <- c(">Bdut_CR2A_4", ">Bdut_CR2A_5")

#==================================================================================================================
#Data Filtering 

#================================
#Breaking down annotation for protein seq
#================================

df_prot_sep <- df_prot %>% 
  separate(col = Seq, into = c("Species", "Strain", "Plasmid"), sep = "_", remove = FALSE)

#================================
#Writing out species name
#================================

df_prot_sep <- df_prot_sep %>% 
  mutate(
    Name = case_when(
      str_starts(Species, ">Bbur") ~ "burgdorferi",
      str_starts(Species, ">Bmayo") ~ "mayonii",
      str_starts(Species, ">Bbis") ~ "bissettiae",
      str_starts(Species, ">Bbav") ~ "bavariensis",
      str_starts(Species, ">Bafz") ~ "afzelii",
      str_starts(Species, ">Bval") ~ "valaisiana",
      str_starts(Species, ">Bturd") ~ "turdi",
      str_starts(Species, ">Bjap") ~ "japonica",
      str_starts(Species, ">Bturc") ~ "turcica",
      str_starts(Species, ">Bmiya") ~ "miyamotoi",
      str_starts(Species, ">Bherm") ~ "hermsii",
      str_starts(Species, ">Bniet") ~ "nietonii",
      str_starts(Species, ">Bpark") ~ "parkeri",
      str_starts(Species, ">Bcroc") ~ "crocidurae",
      str_starts(Species, ">Bdut") ~ "duttonii",
      str_starts(Species, ">Brec") ~ "recurrentis"
    )
  )

#==================================================================================================================
#First Iteration with seqdump1

#================================
#Function for Relating Prot to DNA
#================================

df_prot_sep <- df_prot_sep %>% 
  mutate(DNA_seq = NA)
df_prot_sep_trial <- df_prot_sep

for (i in 1:length(df_seqdump1[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    if (str_detect(df_seqdump1[i,1], df_prot_sep[j,3]) & 
        str_detect(df_seqdump1[i,1], paste("32-", df_prot_sep[j,4], sep = "")) & 
        str_detect(df_seqdump1[i,1], df_prot_sep[j,5]) & 
        is.na(df_prot_sep[j,6])) {
      df_prot_sep[j,6] <- df_seqdump1[i,1]
    }
  }
}

#==================================================================================================================
#Second Iteration with seqdump2

#================================
#Function for Relating Prot to DNA
#================================

for (i in 1:length(df_seqdump2[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    if (str_detect(df_seqdump2[i,1], df_prot_sep[j,3]) & 
        str_detect(df_seqdump2[i,1], paste("32-", df_prot_sep[j,4], sep = "")) & 
        str_detect(df_seqdump2[i,1], df_prot_sep[j,5]) & 
        is.na(df_prot_sep[j,6])) {
      df_prot_sep[j,6] <- df_seqdump2[i,1]
    }
  }
}

#==================================================================================================================
#Third Iteration with seqdump3

#================================
#Function for Relating Prot to DNA
#================================

for (i in 1:length(df_seqdump3[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    if (str_detect(df_seqdump3[i,1], df_prot_sep[j,3]) & 
        str_detect(df_seqdump3[i,1], df_prot_sep[j,4]) & 
        str_detect(df_seqdump3[i,1], df_prot_sep[j,5]) & 
        is.na(df_prot_sep[j,6])) {
      df_prot_sep[j,6] <- df_seqdump3[i,1]
    }
  }
}

#==================================================================================================================
#Fourth Iteration with seqdump4

#================================
#Function for Relating Prot to DNA
#================================

for (i in 1:length(df_seqdump4[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    if (str_detect(df_seqdump4[i,1], df_prot_sep[j,3]) & 
        str_detect(df_seqdump4[i,1], df_prot_sep[j,4]) & 
        str_detect(df_seqdump4[i,1], df_prot_sep[j,5]) & 
        is.na(df_prot_sep[j,6])) {
      df_prot_sep[j,6] <- df_seqdump4[i,1]
    }
  }
}
