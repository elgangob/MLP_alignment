# 1. Load the visualization and data wrangling libraries
library(ggplot2)
library(dplyr)

# 2. Define the exact file path to your ClonalFrameML text output
base_dir <- "~/Documents/Stevenson_Project/Aligned_+_Filter_(For_Phylogenetic)/Short Aligned + Filtered/ClonalFrameML/"
import_file <- file.path(base_dir, "ClonalFrameML_Out.importation_status.txt")

# 3. Read the data
imports <- read.table(import_file, header = TRUE)

# 4. Data Wrangling: Calculate Recombination per Branch
recomb_data <- imports %>%
  # Calculate the physical size of each swapped cassette
  mutate(Import_Length = End - Beg) %>%
  # Group the data by each specific branch/strain
  group_by(Node) %>%
  # Sum up total bp and count how many separate swaps happened
  summarise(Total_Recombined_bp = sum(Import_Length),
            Event_Count = n()) %>%
  ungroup() %>%
  # Sort from highest recombination to lowest for a clean aesthetic
  arrange(Total_Recombined_bp)

# Lock in the sorting order so ggplot doesn't scramble it alphabetically
recomb_data$Node <- factor(recomb_data$Node, levels = recomb_data$Node)

# 5. Create the customized Shelter Island highlight logic for the y-axis
# This creates a vector of colors/fonts that spots B31 and 297 automatically
axis_colors <- ifelse(grepl("B31|297", levels(recomb_data$Node)), "firebrick", "black")
axis_fonts  <- ifelse(grepl("B31|297", levels(recomb_data$Node)), "bold", "plain")

# 6. Generate the Aesthetic Lollipop Chart
final_plot <- ggplot(recomb_data, aes(x = Node, y = Total_Recombined_bp)) +
  # Draw the 'sticks' of the lollipop
  geom_segment(aes(x = Node, xend = Node, y = 0, yend = Total_Recombined_bp), 
               color = "gray70", linewidth = 0.8) +
  # Draw the 'candy' on the ends, sized by number of events and colored by impact
  geom_point(aes(size = Event_Count, color = Total_Recombined_bp), alpha = 0.9) +
  # Apply a sleek, colorblind-friendly gradient
  scale_color_viridis_c(option = "plasma", guide = "none") +
  # Flip the chart sideways so strain names are readable
  coord_flip() +
  # Clean up the labels
  labs(title = "Recombination Impact per Borrelia Branch",
       subtitle = "Total imported DNA base pairs across Mlp short sequences",
       x = "Strain / Ancestral Node",
       y = "Total Recombined DNA (bp)",
       size = "Number of\nImport Events") +
  # Apply a minimalist publication theme
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray30", size = 12, margin = margin(b = 15)),
    # Apply the custom colors/fonts to the y-axis names
    axis.text.y = element_text(size = 9, color = axis_colors, face = axis_fonts),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    panel.grid.major.y = element_blank(), # Removes horizontal lines for a cleaner look
    legend.position = "right"
  )

# 7. Print the plot to your RStudio Viewer
print(final_plot)

# 8. Save a high-resolution PDF
# 8. Save a high-resolution, elongated PDF
ggsave(
  filename = file.path(base_dir, "Mlp_Recombination_Lollipop_Chart.pdf"), 
  plot = final_plot, 
  width = 10, 
  height = 40,        
  limitsize = FALSE,  # Tells ggplot not to panic about the massive size
  dpi = 300
)

#==========================================================================
library(ape)

# 1. Define your file paths
base_dir <- "~/Documents/Stevenson_Project/Aligned_+_Filter_(For_Phylogenetic)/Short Aligned + Filtered/ClonalFrameML"
tree_file <- file.path(base_dir, "ClonalFrameML_Out.labelled_tree.newick")
output_pdf <- file.path(base_dir, "Large_Phylogeny_Tree.pdf")

# 2. Read the tree
my_tree <- read.tree(tree_file)

# 3. Open a massive PDF canvas
# Height is set to 40 inches to give all ~160+ sequences plenty of vertical breathing room
pdf(file = output_pdf, width = 12, height = 40)

# 4. Adjust the margins (bottom, left, top, right)
# We make the right margin very large (8) so long strain names don't get cut off the page
par(mar = c(2, 1, 2, 8)) 

# 5. Plot the tree
plot(my_tree, 
     type = "phylogram",       # Standard square-branch tree
     cex = 0.8,                # Controls the text size of the strain names (adjust if too big/small)
     edge.color = "black",     # Branch color
     edge.width = 1.5,         # Makes the branches slightly thicker and easier to read
     label.offset = 0.005)     # Adds a tiny gap between the physical branch and the text

# Add a scale bar at the bottom to show mutational distance
add.scale.bar(cex = 0.8)

# 6. Close the PDF writer (This actually saves the file!)
dev.off()

print("PDF successfully generated! Check your folder.")