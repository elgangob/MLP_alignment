#Normal Sequence

# 1. Load the required libraries
library(seqinr)
library(ggplot2)

# 2. Load your alignment
alignment <- read.fasta("Documents/Stevenson Project/Aligned Only/Normal/All Mlp DNA_aligned_only.fas", 
                        forceDNAtolower = FALSE)

# Read the .rec file, completely ignoring the broken headers (header = FALSE, skip = 1)
rec_table <- read.delim("Coding/3seq build 170612/3seq_Normal.3s.rec", 
                        header = FALSE, 
                        skip = 1,
                        stringsAsFactors = FALSE)

# 3. Sort the table by the corrected p-values (which live in Column 11)
rec_table <- rec_table[order(as.numeric(rec_table[, 11])), ]

# 4. Scrub invisible spaces using raw column numbers
names(alignment) <- trimws(names(alignment))
rec_table[, 1] <- trimws(rec_table[, 1]) # Parent 1
rec_table[, 2] <- trimws(rec_table[, 2]) # Parent 2
rec_table[, 3] <- trimws(rec_table[, 3]) # Child

# 5. Set parameters
num_plots_to_make <- 20 
window_size <- 50
step_size <- 10

# 6. Open the PDF
pdf("Documents/Stevenson Project/Mlp_Recombination_Plot_Normal.pdf", width = 10, height = 6)

# 7. The Automation Loop
for (row in 1:num_plots_to_make) {
  
  # Pull the names and p-value using raw column numbers!
  parent1_name <- rec_table[row, 1]
  parent2_name <- rec_table[row, 2]
  child_name <- rec_table[row, 3]
  p_value <- signif(as.numeric(rec_table[row, 11]), digits = 3) 
  
  # Safety check
  if(parent1_name %in% names(alignment) && 
     parent2_name %in% names(alignment) && 
     child_name %in% names(alignment)) {
    
    child_seq <- alignment[[child_name]]
    parent1_seq <- alignment[[parent1_name]]
    parent2_seq <- alignment[[parent2_name]]
    
    positions <- c()
    sim_to_p1 <- c()
    sim_to_p2 <- c()
    
    for (i in seq(1, length(child_seq) - window_size, by = step_size)) {
      window_indices <- i:(i + window_size - 1)
      
      window_child <- child_seq[window_indices]
      window_p1 <- parent1_seq[window_indices]
      window_p2 <- parent2_seq[window_indices]
      
      match_p1 <- sum(window_child == window_p1) / window_size
      match_p2 <- sum(window_child == window_p2) / window_size
      
      positions <- c(positions, i + (window_size / 2))
      sim_to_p1 <- c(sim_to_p1, match_p1 * 100)
      sim_to_p2 <- c(sim_to_p2, match_p2 * 100)
    }
    
    plot_data <- data.frame(
      Position = rep(positions, 2),
      Similarity = c(sim_to_p1, sim_to_p2),
      Parent = rep(c(paste("Parent 1:", parent1_name), 
                     paste("Parent 2:", parent2_name)), each = length(positions))
    )
    
    p <- ggplot(plot_data, aes(x = Position, y = Similarity, color = Parent)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      labs(
        title = paste("Recombinant:", child_name),
        subtitle = paste("p-value =", p_value),
        x = "Nucleotide Position",
        y = "Sequence Identity (%)"
      ) +
      scale_color_manual(values = c("blue", "red"))
    
    print(p)
    
  } else {
    message(paste("Skipping row", row, "- Sequence names not found in alignment."))
  }
}

# 8. Close and save the PDF document
dev.off()

print("All plots successfully saved to Mlp_Recombination_Plots.pdf!")

#====================================================================================
#Long Sequence

# 2. Load your alignment
alignment <- read.fasta("Documents/Stevenson Project/Aligned Only/Long/All Mlp DNA_long_aligned.fas", 
                        forceDNAtolower = FALSE)

# Read the .rec file, completely ignoring the broken headers (header = FALSE, skip = 1)
rec_table <- read.delim("Coding/3seq build 170612/3seq_Long.3s.rec", 
                        header = FALSE, 
                        skip = 1,
                        stringsAsFactors = FALSE)

# 3. Sort the table by the corrected p-values (which live in Column 11)
rec_table <- rec_table[order(as.numeric(rec_table[, 11])), ]

# 4. Scrub invisible spaces using raw column numbers
names(alignment) <- trimws(names(alignment))
rec_table[, 1] <- trimws(rec_table[, 1]) # Parent 1
rec_table[, 2] <- trimws(rec_table[, 2]) # Parent 2
rec_table[, 3] <- trimws(rec_table[, 3]) # Child

# 5. Set parameters
num_plots_to_make <- 20 
window_size <- 50
step_size <- 10

# 6. Open the PDF
pdf("Documents/Stevenson Project/Mlp_Recombination_Plots_Long.pdf", width = 10, height = 6)

# 7. The Automation Loop
for (row in 1:num_plots_to_make) {
  
  # Pull the names and p-value using raw column numbers!
  parent1_name <- rec_table[row, 1]
  parent2_name <- rec_table[row, 2]
  child_name <- rec_table[row, 3]
  p_value <- signif(as.numeric(rec_table[row, 11]), digits = 3) 
  
  # Safety check
  if(parent1_name %in% names(alignment) && 
     parent2_name %in% names(alignment) && 
     child_name %in% names(alignment)) {
    
    child_seq <- alignment[[child_name]]
    parent1_seq <- alignment[[parent1_name]]
    parent2_seq <- alignment[[parent2_name]]
    
    positions <- c()
    sim_to_p1 <- c()
    sim_to_p2 <- c()
    
    for (i in seq(1, length(child_seq) - window_size, by = step_size)) {
      window_indices <- i:(i + window_size - 1)
      
      window_child <- child_seq[window_indices]
      window_p1 <- parent1_seq[window_indices]
      window_p2 <- parent2_seq[window_indices]
      
      match_p1 <- sum(window_child == window_p1) / window_size
      match_p2 <- sum(window_child == window_p2) / window_size
      
      positions <- c(positions, i + (window_size / 2))
      sim_to_p1 <- c(sim_to_p1, match_p1 * 100)
      sim_to_p2 <- c(sim_to_p2, match_p2 * 100)
    }
    
    plot_data <- data.frame(
      Position = rep(positions, 2),
      Similarity = c(sim_to_p1, sim_to_p2),
      Parent = rep(c(paste("Parent 1:", parent1_name), 
                     paste("Parent 2:", parent2_name)), each = length(positions))
    )
    
    p <- ggplot(plot_data, aes(x = Position, y = Similarity, color = Parent)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      labs(
        title = paste("Recombinant:", child_name),
        subtitle = paste("p-value =", p_value),
        x = "Nucleotide Position",
        y = "Sequence Identity (%)"
      ) +
      scale_color_manual(values = c("blue", "red"))
    
    print(p)
    
  } else {
    message(paste("Skipping row", row, "- Sequence names not found in alignment."))
  }
}

# 8. Close and save the PDF document
dev.off()

print("All plots successfully saved to Mlp_Recombination_Plots.pdf!")

#====================================================================================
#Short Sequence

# 2. Load your alignment
alignment <- read.fasta("Documents/Stevenson Project/Aligned Only/Short/All Mlp DNA_short_aligned.fas", 
                        forceDNAtolower = FALSE)

# Read the .rec file, completely ignoring the broken headers (header = FALSE, skip = 1)
rec_table <- read.delim("Coding/3seq build 170612/3seq_Short.3s.rec", 
                        header = FALSE, 
                        skip = 1,
                        stringsAsFactors = FALSE)

# 3. Sort the table by the corrected p-values (which live in Column 11)
rec_table <- rec_table[order(as.numeric(rec_table[, 11])), ]

# 4. Scrub invisible spaces using raw column numbers
names(alignment) <- trimws(names(alignment))
rec_table[, 1] <- trimws(rec_table[, 1]) # Parent 1
rec_table[, 2] <- trimws(rec_table[, 2]) # Parent 2
rec_table[, 3] <- trimws(rec_table[, 3]) # Child

# 5. Set parameters
num_plots_to_make <- 20 
window_size <- 50
step_size <- 10

# 6. Open the PDF
pdf("Documents/Stevenson Project/Mlp_Recombination_Plots_Short.pdf", width = 10, height = 6)

# 7. The Automation Loop
for (row in 1:num_plots_to_make) {
  
  # Pull the names and p-value using raw column numbers!
  parent1_name <- rec_table[row, 1]
  parent2_name <- rec_table[row, 2]
  child_name <- rec_table[row, 3]
  p_value <- signif(as.numeric(rec_table[row, 11]), digits = 3) 
  
  # Safety check
  if(parent1_name %in% names(alignment) && 
     parent2_name %in% names(alignment) && 
     child_name %in% names(alignment)) {
    
    child_seq <- alignment[[child_name]]
    parent1_seq <- alignment[[parent1_name]]
    parent2_seq <- alignment[[parent2_name]]
    
    positions <- c()
    sim_to_p1 <- c()
    sim_to_p2 <- c()
    
    for (i in seq(1, length(child_seq) - window_size, by = step_size)) {
      window_indices <- i:(i + window_size - 1)
      
      window_child <- child_seq[window_indices]
      window_p1 <- parent1_seq[window_indices]
      window_p2 <- parent2_seq[window_indices]
      
      match_p1 <- sum(window_child == window_p1) / window_size
      match_p2 <- sum(window_child == window_p2) / window_size
      
      positions <- c(positions, i + (window_size / 2))
      sim_to_p1 <- c(sim_to_p1, match_p1 * 100)
      sim_to_p2 <- c(sim_to_p2, match_p2 * 100)
    }
    
    plot_data <- data.frame(
      Position = rep(positions, 2),
      Similarity = c(sim_to_p1, sim_to_p2),
      Parent = rep(c(paste("Parent 1:", parent1_name), 
                     paste("Parent 2:", parent2_name)), each = length(positions))
    )
    
    p <- ggplot(plot_data, aes(x = Position, y = Similarity, color = Parent)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      labs(
        title = paste("Recombinant:", child_name),
        subtitle = paste("p-value =", p_value),
        x = "Nucleotide Position",
        y = "Sequence Identity (%)"
      ) +
      scale_color_manual(values = c("blue", "red"))
    
    print(p)
    
  } else {
    message(paste("Skipping row", row, "- Sequence names not found in alignment."))
  }
}

# 8. Close and save the PDF document
dev.off()

print("All plots successfully saved to Mlp_Recombination_Plots.pdf!")