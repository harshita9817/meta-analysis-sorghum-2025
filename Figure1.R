# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(colorspace)
library(patchwork)

#----------------------------#
# PANEL A: Trait Counts
#----------------------------#

# Load new and old trait datasets
newpaper <- read.csv("new_data_updated.csv")
oldpaper <- read.csv("old_data_updated.csv")

# Extract phenotype classes from first row (excluding Genotype)
old_traits <- oldpaper[1, -1]
new_traits <- newpaper[1, -1]

# Create trait dataframes
old_df <- data.frame(
  Trait = colnames(old_traits),
  PhenotypeClass = as.character(old_traits[1, ]),
  Source = "Old"
)

new_df <- data.frame(
  Trait = colnames(new_traits),
  PhenotypeClass = as.character(new_traits[1, ]),
  Source = "New"
)

# Combine and count
all_traits <- bind_rows(old_df, new_df)
trait_counts <- all_traits %>%
  filter(!is.na(PhenotypeClass)) %>%
  group_by(PhenotypeClass, Source) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(ClassSource = paste0(PhenotypeClass, "_", Source))

# Color map (Old = dark, New = light per class)
classes <- unique(trait_counts$PhenotypeClass)
base_colors <- qualitative_hcl(length(classes), palette = "Dynamic")
names(base_colors) <- classes
color_map <- unlist(lapply(classes, function(cls) {
  base_col <- base_colors[cls]
  c(darken(base_col, 0.3), lighten(base_col, 0.3))
}))
names(color_map) <- c(rbind(paste0(classes, "_Old"), paste0(classes, "_New")))

# Panel A plot
panel_a_plot <- ggplot(trait_counts, aes(x = PhenotypeClass, y = Count, fill = ClassSource)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3.2, color = "black") +
  scale_fill_manual(values = color_map) +
  labs(title = "A", x = "Phenotype Class", y = "Number of Traits") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

#----------------------------#
# PANEL B: Correlation Heatmap
#----------------------------#

# Load merged trait data with class in first row
merged <- read.csv("Merged_data_updated.csv")

# Extract trait classes and numeric values
trait_classes <- as.character(merged[1, -1])
names(trait_classes) <- colnames(merged)[-1]
trait_data <- merged[-1, ]
rownames(trait_data) <- trait_data$Genotype
trait_data <- trait_data[, -1]
trait_data_numeric <- trait_data %>% mutate(across(everything(), ~as.numeric(as.character(.))))
trait_data_numeric <- trait_data_numeric[, colSums(!is.na(trait_data_numeric)) > 0]
trait_data_numeric <- trait_data_numeric[, apply(trait_data_numeric, 2, function(x) sd(x, na.rm = TRUE) > 0)]

# Correlation matrix
cor_matrix <- cor(trait_data_numeric, use = "pairwise.complete.obs")
trait_classes <- trait_classes[colnames(cor_matrix)]

# Annotation and color setup
annotation <- data.frame(Class = trait_classes)
rownames(annotation) <- names(trait_classes)
ordered_classes <- sort(unique(trait_classes))
class_colors <- qualitative_hcl(length(ordered_classes), palette = "Dynamic")
names(class_colors) <- ordered_classes
annotation_colors <- list(Class = class_colors)

# Create heatmap
heatmap_obj <- pheatmap(
  mat = cor_matrix,
  annotation_col = annotation,
  annotation_row = annotation,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  border_color = NA,
  fontsize = 12,
  main = "B",
  silent = TRUE
)
panel_b_plot <- heatmap_obj[[4]]

#----------------------------#
# PANEL C: Trait h² Boxplot
#----------------------------#

heritability <- read.csv("Vg_Ve_Heritability_with_PhenoClass.csv")
heritability$PhenotypeClass <- factor(heritability$PhenotypeClass, levels = sort(unique(heritability$PhenotypeClass)))

panel_c_plot <- ggplot(heritability, aes(x = PhenotypeClass, y = h2, fill = PhenotypeClass)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 1.8, color = "black", alpha = 0.9) +
  geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.5) +
  scale_fill_manual(values = class_colors) +
  labs(title = "C", x = "Phenotype Class", y = expression(h^2)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(face = "bold", hjust = 0, size = 16),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

#----------------------------#
# PANEL D: PC h² Barplot
#----------------------------#

pcheritability <- read.csv("pc_heritability_output.csv")
pcheritability$Trait <- factor(pcheritability$Trait, levels = rev(pcheritability$Trait))

panel_d_plot <- ggplot(pcheritability, aes(x = h2, y = Trait)) +
  geom_bar(stat = "identity", fill = "#4682B4", color = "black", width = 0.6) +
  geom_text(aes(label = round(h2, 2)), hjust = -0.1, size = 3.5, color = "black") +
  labs(title = "D", x = expression(h^2), y = "Principal Components") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(face = "bold", hjust = 0, size = 16),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  xlim(0, 1.05)

#----------------------------#
# COMBINE PANELS A–D
#----------------------------#

combined_plot <- (panel_a_plot | panel_b_plot) /
  (panel_c_plot | panel_d_plot)

# Save combined figure
ggsave("Figure1_AllPanels.svg", plot = combined_plot, width = 14, height = 12, dpi = 300)
ggsave("Figure1_AllPanels.png", plot = combined_plot, width = 14, height = 12, dpi = 300)

