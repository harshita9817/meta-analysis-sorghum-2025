# =======================
# Panel A — GWAS RMIP Manhattan (across PCs)
# =======================

# ---- Packages ----
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)

# ---- Params ----
rmip_thr   <- 0.5
pcs        <- paste0("PC", 1:10)

# ---- Helpers ----
std_cols <- function(df) {
  names(df) <- trimws(names(df))
  # CHROM
  if (!"CHROM" %in% names(df)) {
    cand <- intersect(names(df), c("Chromosome","CHR","Chr","chrom","Chrom"))
    if (length(cand)) names(df)[names(df) == cand[1]] <- "CHROM"
  }
  # POS
  if (!"POS" %in% names(df)) {
    cand <- intersect(names(df), c("Position.","Position","BP","pos"))
    if (length(cand)) names(df)[names(df) == cand[1]] <- "POS"
  }
  # RMIP
  if (!"RMIP" %in% names(df)) {
    cand <- intersect(names(df), c("Rmip","rmip"))
    if (length(cand)) names(df)[names(df) == cand[1]] <- "RMIP"
  }
  df
}

find_pc_file <- function(pc) {
  cands <- c(
    paste0("Z", pc, "_signals.csv"),
    paste0("Z", pc, "signals.csv"),
    paste0(pc, "_signals.csv")
  )
  existing <- cands[file.exists(cands)]
  if (length(existing)) existing[1] else NA_character_
}

# ---- Step A1: Load & stack GWAS results for PCs ----
all_gwas_list <- lapply(pcs, function(pc) {
  f <- find_pc_file(pc)
  if (is.na(f)) return(NULL)
  df <- fread(f, data.table = FALSE)
  df <- std_cols(df)
  df$Trait <- pc
  df
})

all_gwas <- bind_rows(all_gwas_list)
if (nrow(all_gwas) == 0) stop("No GWAS PC signal files found. Check filenames like 'ZPC1_signals.csv'.")

# ---- Step A2: Clean CHROM/POS, keep chr 1–10 ----
all_gwas <- all_gwas %>%
  mutate(
    CHROM = gsub("^Chr", "", as.character(CHROM)),
    CHROM = gsub("^0",   "", CHROM),
    CHROM = factor(CHROM, levels = as.character(1:10), ordered = TRUE),
    POS   = as.numeric(POS)
  ) %>%
  filter(CHROM %in% as.character(1:10), !is.na(POS)) %>%
  arrange(CHROM, POS)

# ---- Step A3: Keep only PCs with any RMIP > threshold ----
pcs_with_hits <- all_gwas %>%
  group_by(Trait) %>%
  summarise(any_hit = any(RMIP > rmip_thr, na.rm = TRUE), .groups = "drop") %>%
  filter(any_hit) %>%
  pull(Trait)

filtered_gwas <- all_gwas %>% filter(Trait %in% pcs_with_hits)
if (nrow(filtered_gwas) == 0) stop("No PCs have RMIP > 0.5. Try lowering the threshold or check your data.")

# ---- Step A4: Chromosome offsets & BPcum ----
chrom_lengths <- filtered_gwas %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(POS, na.rm = TRUE), .groups = "drop")

chrom_offsets <- chrom_lengths %>%
  mutate(offset = dplyr::lag(cumsum(chr_len), default = 0)) %>%
  dplyr::select(CHROM, offset)

filtered_gwas <- filtered_gwas %>%
  left_join(chrom_offsets, by = "CHROM") %>%
  mutate(BPcum = POS + offset) %>%
  arrange(CHROM, POS)

# ---- Step A5: Axis centers ----
axis.set <- filtered_gwas %>%
  group_by(CHROM) %>%
  summarise(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop")

# ---- Step A6: Colors (fixed mapping PC1..PC10; subset to present PCs) ----
pc_colors_full <- setNames(brewer.pal(10, "Paired"), paste0("PC", 1:10))
filtered_gwas$Trait <- factor(filtered_gwas$Trait, levels = paste0("PC", 1:10))
pc_levels_present <- levels(droplevels(filtered_gwas$Trait))
pc_colors <- pc_colors_full[pc_levels_present]

# (save PC1 color for later panels, if needed)
pc1_col <- unname(pc_colors_full["PC1"])

# ---- Step A7: Plot Panel A ----
panelA <- ggplot(filtered_gwas, aes(x = BPcum, y = RMIP, color = Trait)) +
  geom_point(size = 3.4, alpha = 0.85) +
  geom_hline(yintercept = rmip_thr, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = pc_colors, name = "PC Component", breaks = pc_levels_present) +
  scale_x_continuous(breaks = axis.set$center, labels = as.character(axis.set$CHROM)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
  labs(x = "Chromosome", y = "RMIP", title = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background  = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    axis.line         = element_line(color = "black"),
    axis.text.x       = element_text(color = "black"),
    axis.text.y       = element_text(color = "black"),
    legend.title      = element_text(face = "bold"),
    legend.key        = element_blank(),
    legend.position   = "bottom"
  )

# ---- Step A8: Save (white background) ----
ggsave("PanelA_GWAS_PC_RMIP.png", panelA, width = 12, height = 4.8, dpi = 400, bg = "white")
########################################
# =======================
# Panel B — PC1 only (uses PC1 color from Panel A)
# =======================

library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)

# If pc1_col wasn't created in Panel A for some reason, define it here:
if (!exists("pc1_col")) {
  pc_colors_full <- setNames(brewer.pal(10, "Paired"), paste0("PC", 1:10))
  pc1_col <- unname(pc_colors_full["PC1"])
}

# ------- Load PC1 data -------
pc1_df <- fread("GWAS_PC1.csv", data.table = FALSE)

# Harmonize columns
names(pc1_df) <- trimws(names(pc1_df))
if (!"CHROM" %in% names(pc1_df)) {
  cand <- intersect(names(pc1_df), c("Chromosome","CHR","Chr","chrom","Chrom"))
  if (length(cand)) names(pc1_df)[names(pc1_df) == cand[1]] <- "CHROM"
}
if (!"POS" %in% names(pc1_df)) {
  cand <- intersect(names(pc1_df), c("Position.","Position","BP","pos"))
  if (length(cand)) names(pc1_df)[names(pc1_df) == cand[1]] <- "POS"
}
# Compute FDR if missing
if (!"FDR" %in% names(pc1_df)) {
  if (!"P.value" %in% names(pc1_df)) stop("GWAS_PC1.csv must have FDR or P.value.")
  pc1_df$FDR <- p.adjust(pc1_df$P.value, method = "BH")
}

# Clean & order
pc1_df <- pc1_df %>%
  mutate(
    CHROM  = gsub("^Chr", "", as.character(CHROM)),
    CHROM  = gsub("^0",   "", CHROM),
    CHROM  = factor(CHROM, levels = as.character(1:10), ordered = TRUE),
    POS    = as.numeric(POS),
    logFDR = -log10(pmax(FDR, .Machine$double.xmin))
  ) %>%
  filter(CHROM %in% as.character(1:10), !is.na(POS)) %>%
  arrange(CHROM, POS)

# ------- Cumulative genome position (reuse Panel A spacing if available) -------
if (exists("chrom_offsets")) {
  pc1_df <- pc1_df %>%
    left_join(chrom_offsets, by = "CHROM") %>%
    mutate(offset = ifelse(is.na(offset), 0, offset),
           BPcum  = POS + offset)
} else {
  chrom_lengths_b <- pc1_df %>% group_by(CHROM) %>%
    summarise(chr_len = max(POS, na.rm = TRUE), .groups = "drop")
  chrom_offsets_b <- chrom_lengths_b %>%
    mutate(offset = dplyr::lag(cumsum(chr_len), default = 0)) %>%
    dplyr::select(CHROM, offset)
  pc1_df <- pc1_df %>%
    left_join(chrom_offsets_b, by = "CHROM") %>%
    mutate(BPcum = POS + offset)
}

# Axis centers: reuse Panel A if present; else compute
if (!exists("axis.set")) {
  axis.set <- pc1_df %>%
    group_by(CHROM) %>%
    summarise(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop")
}

# ------- Plot -------
thr_line <- -log10(0.05)
maxFDR_B <- max(pc1_df$logFDR, na.rm = TRUE)
maxFDR_B <- ifelse(maxFDR_B > thr_line, ceiling(maxFDR_B), ceiling(thr_line))

panelB <- ggplot(pc1_df, aes(x = BPcum, y = logFDR)) +
  geom_point(size = 1.8, alpha = 0.9, color = pc1_col) +  # same PC1 color as Panel A
  geom_hline(yintercept = thr_line, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = axis.set$center, labels = as.character(axis.set$CHROM)) +
  scale_y_continuous(limits = c(0, maxFDR_B), expand = c(0, 0.01)) +
  labs(x = "Chromosome", y = expression(-log[10](FDR)), title = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background  = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    axis.line         = element_line(color = "black"),
    axis.text.x       = element_text(color = "black"),
    axis.text.y       = element_text(color = "black"),
    legend.position   = "none"
  )

# Save Panel B
ggsave("PanelB_PC1_logFDR.png", panelB, width = 12, height = 4.5, dpi = 400, bg = "white")
#ggsave("PanelB_PC1_logFDR.svg", panelB, width = 12, height = 4.5, bg = "white")
# =======================
# Panel C — Loadings (top +/− per PC) with your class palette
# =======================
# ===== Panel C — Loadings (top ±5 per PC) =====

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(colorspace)
loadings<- read.csv("TraitLoadings_Top10PCs_Winsorized_Annotated.csv")
# If your class column is named differently, normalize it:
if (!"Class" %in% names(loadings) && "PhenotypeClass" %in% names(loadings)) {
  loadings <- dplyr::rename(loadings, Class = PhenotypeClass)
}

# PCs to show (change as needed)
pcs_to_plot <- c("PC1","PC2","PC3","PC6")
stopifnot(all(c("Trait","Class", pcs_to_plot) %in% names(loadings)))

# ---- Class colors: same mapping, darker shade ----
library(colorspace)

ordered_classes <- sort(c("Agronomic","Biochemical","Disease","Reproductive","Root","Seed","Vegetative"))

# darker Dynamic palette (tweak l/c if you want)
class_colors <- qualitative_hcl(
  length(ordered_classes),
  palette = "Dynamic",
  l = 40,   # lower = darker (try 35–45)
  c = 90    # chroma; bump a bit to keep vividness
)
names(class_colors) <- ordered_classes

# If your data has classes beyond that list, extend the palette but keep existing colors unchanged
extra_classes <- setdiff(unique(loadings$Class), ordered_classes)
if (length(extra_classes) > 0) {
  extra_cols <- qualitative_hcl(length(extra_classes), palette = "Dynamic", l = 40, c = 90)
  names(extra_cols) <- extra_classes
  class_colors <- c(class_colors, extra_cols)
}
loadings$Class <- factor(loadings$Class, levels = names(class_colors))

# ---- Your loadings logic (top ±5) → Panel C ----
top_k_each <- 5  # how many positives and negatives per PC
if (!exists("pcs_to_plot")) pcs_to_plot <- c("PC1","PC2","PC3","PC6")

load_long <- loadings %>%
  dplyr::select(Trait, Class, dplyr::all_of(pcs_to_plot)) %>%
  tidyr::pivot_longer(dplyr::all_of(pcs_to_plot),
                      names_to = "PC", values_to = "Loading") %>%
  dplyr::mutate(PC = factor(PC, levels = pcs_to_plot, ordered = TRUE))

top_extremes <- load_long %>%
  dplyr::group_by(PC) %>%
  dplyr::mutate(
    rank_pos = dplyr::if_else(Loading > 0, dplyr::dense_rank(dplyr::desc(Loading)), NA_integer_),
    rank_neg = dplyr::if_else(Loading < 0, dplyr::dense_rank(Loading),              NA_integer_)
  ) %>%
  dplyr::filter(
    (!is.na(rank_pos) & rank_pos <= top_k_each) |
      (!is.na(rank_neg) & rank_neg <= top_k_each)
  ) %>%
  dplyr::ungroup()

Lmax <- max(abs(top_extremes$Loading), na.rm = TRUE)

top_extremes <- top_extremes %>%
  dplyr::group_by(PC) %>%
  dplyr::arrange(Loading, .by_group = TRUE) %>%
  dplyr::mutate(
    trait_id = paste(PC, Trait, sep = " - "),
    trait_id = forcats::fct_inorder(trait_id)
  ) %>%
  dplyr::ungroup()

panelC <- ggplot(top_extremes,
                 aes(x = trait_id, y = Loading, fill = Class)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray35") +
  geom_col(width = 0.75, color = "black", linewidth = 0.25, alpha = 1) +  # keep alpha=1 for darker look
  scale_fill_manual(values = class_colors, name = "Class") +
  scale_y_continuous(limits = c(-Lmax, Lmax)) +
  scale_x_discrete(labels = function(x) sub("^PC\\d+\\s*-\\s*", "", x)) +
  labs(x = "Trait", y = "PC Loading", title = NULL) +
  facet_grid(. ~ PC, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background  = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    axis.text.x       = element_text(angle = 60, hjust = 1, size = 8),
    legend.position   = "bottom",
    strip.text        = element_text(face = "bold", size = 12, color = "black"),
    strip.background  = element_blank()
  )

# Save if you want
ggsave("PanelC_Loadings.png", panelC, width = 16, height = 7, dpi = 400, bg = "white")
ggsave("PanelC_Loadings.svg", panelC, width = 16, height = 7, bg = "white")


#############I WANTED TO FIX SOMETHING WITH NAMES HERE 

# If 'top_extremes' doesn't exist yet, build it quickly (same logic as your Panel C)
if (!exists("top_extremes")) {
  if (!exists("pcs_to_plot")) pcs_to_plot <- c("PC1","PC2","PC3","PC6")
  top_k_each <- 5
  
  load_long <- loadings %>%
    dplyr::select(Trait, Class, dplyr::all_of(pcs_to_plot)) %>%
    tidyr::pivot_longer(dplyr::all_of(pcs_to_plot),
                        names_to = "PC", values_to = "Loading") %>%
    dplyr::mutate(PC = factor(PC, levels = pcs_to_plot, ordered = TRUE))
  
  top_extremes <- load_long %>%
    dplyr::group_by(PC) %>%
    dplyr::mutate(
      rank_pos = dplyr::if_else(Loading > 0, dplyr::dense_rank(dplyr::desc(Loading)), NA_integer_),
      rank_neg = dplyr::if_else(Loading < 0, dplyr::dense_rank(Loading),              NA_integer_)
    ) %>%
    dplyr::filter(
      (!is.na(rank_pos) & rank_pos <= top_k_each) |
        (!is.na(rank_neg) & rank_neg <= top_k_each)
    ) %>%
    dplyr::ungroup()
}

# Create a template with unique trait names appearing in the plot
label_template <- top_extremes %>%
  dplyr::distinct(Trait) %>%
  dplyr::arrange(Trait) %>%
  dplyr::mutate(New_Label = "")

write.csv(label_template, "TraitLabelTemplate_PanelC.csv", row.names = FALSE)
# -> Edit "TraitLabelTemplate_PanelC.csv" and fill New_Label as you like (you can use \n for line breaks).
# Read edited labels
label_map_df <- read.csv("TraitLabelTemplate_PanelC.csv", stringsAsFactors = FALSE)
label_map_df$New_Label <- trimws(label_map_df$New_Label)
label_map_df$Display <- ifelse(nchar(label_map_df$New_Label) > 0,
                               label_map_df$New_Label, label_map_df$Trait)

# Named vector: Trait -> Display label
trait_label_map <- setNames(label_map_df$Display, label_map_df$Trait)

# Reuse your already-built 'top_extremes' and class_colors
# If you used 'trait_id' like "PC1 - Trait", keep the ordering but show mapped labels on the x-axis
panelC_labeled <- panelC +
  scale_x_discrete(
    labels = function(x) {
      base <- sub("^PC\\d+\\s*-\\s*", "", x)          # strip the "PCx - " part
      mapped <- trait_label_map[base]
      ifelse(is.na(mapped), base, mapped)
    }
  )

# Save updated figure
ggsave("PanelC_Loadings_Labeled.png", panelC_labeled, width = 16, height = 7, dpi = 400, bg = "white")
ggsave("PanelC_Loadings_Labeled.svg", panelC_labeled, width = 16, height = 7, bg = "white")

panelC_labeled

# =======================
# Combine Panel A, B, C
# =======================
library(patchwork)

# If you built a labeled version of Panel C, use it; otherwise fall back to panelC
pC <- if (exists("panelC_labeled")) panelC_labeled else panelC

combined <- panelA / panelB / pC +
  plot_layout(heights = c(1.2, 1.0, 1.0)) +
  plot_annotation(tag_levels = "A")  # adds A, B, C labels

# Preview
combined

# Save
ggsave("Figure_PC_GWAS_TWAS_Loadings.png",
       combined, width = 12, height = 14, dpi = 400, bg = "white")
ggsave("Figure_PC_GWAS_TWAS_Loadings.svg",
       combined, width = 12, height = 14, bg = "white")


















