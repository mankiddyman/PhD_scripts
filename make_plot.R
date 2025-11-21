library(ggplot2)
library(ggtree)
library(treeio)
library(ape)

### ----------------------------------------------------
### 1. LOAD TREE AND TRAITS
### ----------------------------------------------------

traits <- read.csv("traits.csv", stringsAsFactors = FALSE)
tree <- read.tree("angiosperm_families.tre")

# Fix spelling mismatch
tree$tip.label <- gsub("Convolvulvaceae", "Convolvulaceae", tree$tip.label)

### ----------------------------------------------------
### 2. BUILD TRAIT CATEGORY
### ----------------------------------------------------

trait_class <- ifelse(
  traits$carnivory == 1 & traits$holocentric == 1, "both",
  ifelse(traits$carnivory == 1, "carnivorous",
  ifelse(traits$holocentric == 1, "holocentric", "none"))
)

traits$trait_class <- factor(
  trait_class,
  levels = c("none", "carnivorous", "holocentric", "both")
)

### ----------------------------------------------------
### 3. ATTACH TRAITS TO TREE
### ----------------------------------------------------

tree_dat <- data.frame(
  label = tree$tip.label,
  trait_class = traits$trait_class[match(tree$tip.label, traits$clade)],
  stringsAsFactors = FALSE
)

### ----------------------------------------------------
### 4. COLORS FOR TRAIT CATEGORIES
### ----------------------------------------------------

cols <- c(
  none = "grey60",
  carnivorous = "#1f78b4",
  holocentric = "#33a02c",
  both = "#e31a1c"
)

### ----------------------------------------------------
### 5. BASE GGTREE PLOT
### ----------------------------------------------------

p <- ggtree(tree, size = 0.8) %<+% tree_dat +
  geom_tiplab(
    aes(color = trait_class),
    size = 4,
    fontface = "italic"
  ) +
  scale_color_manual(
    values = cols,
    name = "Trait category",
    labels = c(
      none = "\u2022 None",
      carnivorous = "\u2022 Carnivorous",
      holocentric = "\u2022 Holocentric",
      both = "\u2022 Both"
    )
  ) +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 40, 10, 10)
  )

# Add space so tip labels are not cropped
p <- p + xlim(0, max(p$data$x) + 1.5)

### ----------------------------------------------------
### 6. GET INTERNAL NODES BY LABEL (HANDLE DUPLICATES)
### ----------------------------------------------------
# Your tree has internal node labels:
# Asterids, Rosids, Dicots, Monocots, Angiosperms, Root
# Some (like Asterids) appear multiple times (nested clades).
# We’ll take the *last* occurrence, which corresponds to the most inclusive clade.

get_node_by_label <- function(tree, label_name) {
  idx <- which(tree$node.label == label_name)
  if (length(idx) == 0) return(NA)
  idx <- tail(idx, 1)  # use the most inclusive (outermost) one
  idx + length(tree$tip.label)  # convert to actual node number
}

n_asterids   <- get_node_by_label(tree, "Asterids")
n_rosids     <- get_node_by_label(tree, "Rosids")
n_dicots     <- get_node_by_label(tree, "Dicots")
n_monocots   <- get_node_by_label(tree, "Monocots")
n_angios     <- get_node_by_label(tree, "Angiosperms")

### ----------------------------------------------------
### 7. ADD ONLY THE CLADE LABELS YOU WANT
### ----------------------------------------------------

add_label <- function(p, node_id, label) {
  if (is.na(node_id)) return(p)
  
  # extract that node's coordinates from the plot data
  node_df <- p$data[p$data$node == node_id, ]
  if (nrow(node_df) == 0) return(p)
  
  node_df$lab <- label
  
  p + geom_nodelab(
    data = node_df,
    aes(x = x, y = y, label = lab),
    fontface = "bold",
    hjust = -0.2,
    nudge_x = 0.05,
    inherit.aes = FALSE
  )
}

p <- add_label(p, n_asterids,   "Asterids")
p <- add_label(p, n_rosids,     "Rosids")
p <- add_label(p, n_dicots,     "Dicots")
p <- add_label(p, n_monocots,   "Monocots")
p <- add_label(p, n_angios,     "Angiosperms")

### ----------------------------------------------------
### 8. EXPORT PDF
### ----------------------------------------------------

ggsave(
  "Drosera_angio.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

