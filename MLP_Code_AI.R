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

seqdump5 <- read.fasta(file = "Documents/Stevenson Project/seqdump5.txt")

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

df_seqdump5 <- as.data.frame(t(data.frame(Annot = getAnnot(seqdump5))))
rownames(df_seqdump5) <- 1:length(df_seqdump5[,1])
colnames(df_seqdump5) <- "Seq"

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

#================================
# Creating a Plasmid Search Key
#================================

df_prot_sep <- df_prot_sep %>% 
  mutate(
    Plasmid_Search = case_when(
      # 1. Lyme Disease Group (cp32 plasmids)
      Name %in% c("burgdorferi", "mayonii", "bissettiae", "bavariensis", 
                  "afzelii", "valaisiana", "turdi", "japonica") ~ paste("cp32-", Plasmid, sep = ""),
      
      # 2. Borrelia parkeri (SLO) mapping
      Strain == "SLO" & Plasmid == "1" ~ "lp18",
      Strain == "SLO" & Plasmid == "2" ~ "cp28",
      Strain == "SLO" & Plasmid == "3" ~ "cp29",
      Strain == "SLO" & Plasmid == "4" ~ "cp58",
      Strain == "SLO" & Plasmid == "5" ~ "cp58",
      Strain == "SLO" & Plasmid == "6" ~ "lp34-1",
      Strain == "SLO" & Plasmid == "7" ~ "lp60",
      Strain == "SLO" & Plasmid == "8" ~ "lp32",
      Strain == "SLO" & Plasmid == "9" ~ "lp32",
      
      # 3. Borrelia hermsii (DAH) mapping
      Strain == "DAH" & Plasmid == "1" ~ "cp30-2",
      Strain == "DAH" & Plasmid == "2" ~ "cp55",
      Strain == "DAH" & Plasmid == "3" ~ "cp90",
      Strain == "DAH" & Plasmid == "4" ~ "cp90",
      Strain == "DAH" & Plasmid == "5" ~ "lp59",
      Strain == "DAH" & Plasmid == "6" ~ "lp59",
      
      # 4. Borrelia turcica (IST7) mapping
      Strain == "IST7" & Plasmid %in% c("1", "2") ~ "lp129",
      
      # 5. SAFEGUARD FALLBACK
      # For DOU, FR64b, CR2A, etc., where NCBI uses "Contig" or "unnamed".
      # This prevents str_detect from falsely matching the "1" in the NCBI ID.
      TRUE ~ paste0("MANUAL_CHECK_", Plasmid) 
    )
  ) %>%
  # We initialize the DNA_seq column here too to keep things clean
  mutate(DNA_seq = NA)

#==================================================================================================================
# Iterations for Seqdumps
#==================================================================================================================

# First Iteration with seqdump1
for (i in 1:length(df_seqdump1[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    
    current_strain <- df_prot_sep$Strain[j]
    search_term <- df_prot_sep$Plasmid_Search[j]
    
    if (!is.na(search_term) && !str_detect(search_term, "MANUAL_CHECK")) {
      # Added paste0(search_term, "\\b") to prevent "cp32-1" from matching "cp32-10"
      if (str_detect(df_seqdump1[i,1], current_strain) & 
          str_detect(df_seqdump1[i,1], paste0(search_term, "\\b")) & 
          is.na(df_prot_sep$DNA_seq[j])) {
        
        df_prot_sep$DNA_seq[j] <- df_seqdump1[i,1]
      }
    }
  }
}

#==================================================================================================================
# Second Iteration with seqdump2

for (i in 1:length(df_seqdump2[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    
    current_strain <- df_prot_sep$Strain[j]
    search_term <- df_prot_sep$Plasmid_Search[j]
    
    if (!is.na(search_term) && !str_detect(search_term, "MANUAL_CHECK")) {
      if (str_detect(df_seqdump2[i,1], current_strain) & 
          str_detect(df_seqdump2[i,1], paste0(search_term, "\\b")) & 
          is.na(df_prot_sep$DNA_seq[j])) {
        
        df_prot_sep$DNA_seq[j] <- df_seqdump2[i,1]
      }
    }
  }
}

#==================================================================================================================
# Third Iteration with seqdump3

for (i in 1:length(df_seqdump3[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    
    current_strain <- df_prot_sep$Strain[j]
    search_term <- df_prot_sep$Plasmid_Search[j]
    
    if (!is.na(search_term) && !str_detect(search_term, "MANUAL_CHECK")) {
      if (str_detect(df_seqdump3[i,1], current_strain) & 
          str_detect(df_seqdump3[i,1], paste0(search_term, "\\b")) & 
          is.na(df_prot_sep$DNA_seq[j])) {
        
        df_prot_sep$DNA_seq[j] <- df_seqdump3[i,1]
      }
    }
  }
}

#==================================================================================================================
# Fourth Iteration with seqdump4

for (i in 1:length(df_seqdump4[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    
    current_strain <- df_prot_sep$Strain[j]
    search_term <- df_prot_sep$Plasmid_Search[j]
    
    if (!is.na(search_term) && !str_detect(search_term, "MANUAL_CHECK")) {
      if (str_detect(df_seqdump4[i,1], current_strain) & 
          str_detect(df_seqdump4[i,1], paste0(search_term, "\\b")) & 
          is.na(df_prot_sep$DNA_seq[j])) {
        
        df_prot_sep$DNA_seq[j] <- df_seqdump4[i,1]
      }
    }
  }
}


#==================================================================================================================
# Fifth Iteration with seqdump5

for (i in 1:length(df_seqdump5[,1])) {
  for (j in 1:length(df_prot_sep[,1])) {
    
    current_strain <- df_prot_sep$Strain[j]
    search_term <- df_prot_sep$Plasmid_Search[j]
    
    if (!is.na(search_term) && !str_detect(search_term, "MANUAL_CHECK")) {
      if (str_detect(df_seqdump5[i,1], current_strain) & 
          str_detect(df_seqdump5[i,1], paste0(search_term, "\\b")) & 
          is.na(df_prot_sep$DNA_seq[j])) {
        
        df_prot_sep$DNA_seq[j] <- df_seqdump5[i,1]
      }
    }
  }
}