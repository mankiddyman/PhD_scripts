library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggimage)
library(magick)

# ------------------------------------------------------
# Load tree + traits
# ------------------------------------------------------
tree <- read.tree("carnivory_tree.tre")

traits <- data.frame(
  label = tree$tip.label,
  richness = c(12500, 1, 1, 150, 250, 0),
  holocentric = c("mono", "mono", "mono", "mono", "holo", "mono")
)

traits <- traits %>%
  mutate(label2 = paste0(label, " ~", richness))

# ------------------------------------------------------
# Get tree coordinates
# ------------------------------------------------------
p_base <- ggtree(tree) %<+% traits
d <- p_base$data

# ------------------------------------------------------
# USER-ADJUSTABLE WEDGE SCALE
# ------------------------------------------------------
wedge_scale <- 0.003

# ------------------------------------------------------
# Triangle wedge (unchanged)
# ------------------------------------------------------
make_terminal_wedge <- function(d, tip, richness,
                                scale_width = wedge_scale) {

  tip_node <- d$node[d$label == tip & !is.na(d$label)]
  tip_row  <- d[d$node == tip_node, ]
  
  parent_node <- tip_row$parent
  parent_row  <- d[d$node == parent_node, ]
  
  x0 <- parent_row$x
  y0 <- tip_row$y
  
  x1 <- tip_row$x
  y1 <- tip_row$y
  
  W <- sqrt(richness) * scale_width
  
  data.frame(
    x = c(x0, x1, x1),
    y = c(y0, y1 - W, y1 + W),
    tip = tip
  )
}

# ------------------------------------------------------
# Build wedges
# ------------------------------------------------------
wedges <- rbind(
  make_terminal_wedge(d, "Drosera",   250),
  make_terminal_wedge(d, "Nepenthes", 150),
  make_terminal_wedge(d, "Non-carnivorous_caryophyllales", 12500)
)

# ------------------------------------------------------
# IMAGE TABLE
# ------------------------------------------------------
img_df <- data.frame(
  label = c("Dionaea", "Drosera", "Nepenthes",
            "Triphyophyllum", "Drosophyllum",
            "Non-carnivorous_caryophyllales"),
  
  image = c("Dionaea.png",
            "Drosera.jpg",
            "Nepenthes.png",
            "Triphyophyllum.png",
            "Drosophyllum.png",
            "noncarn.png")
)

# Tip y + x positions
tip_coords <- d %>%
  filter(isTip) %>%
  select(label, y, tip_x = x)

traits2 <- traits %>%
  left_join(img_df, by = "label") %>%
  left_join(tip_coords, by = "label")

# ------------------------------------------------------
# IMAGE SIZE NORMALIZATION
# ------------------------------------------------------

# User controls
TARGET_VISUAL_WIDTH  <- 0.14
MAX_SIZE             <- 0.20
HEIGHT_COMPENSATION  <- 1.00

# read pixel dims
get_img_dim <- function(path) {
  img <- magick::image_read(path)
  info <- magick::image_info(img)
  c(width = info$width, height = info$height)
}

dims <- lapply(traits2$image, get_img_dim)
dims <- do.call(rbind, dims)

traits2$pixel_width  <- dims[, 1]
traits2$pixel_height <- dims[, 2]

width_scale <- TARGET_VISUAL_WIDTH / max(traits2$pixel_width)

traits2$img_size <- traits2$pixel_width * width_scale
traits2$img_size <- traits2$img_size * HEIGHT_COMPENSATION
traits2$img_size <- pmin(traits2$img_size, MAX_SIZE)

# ------------------------------------------------------
# CUSTOM MULTIPLIERS PER TAXON (YOUR REQUEST)
# ------------------------------------------------------
traits2 <- traits2 %>%
  mutate(
    img_size = case_when(
      label == "Drosera"                          ~ img_size * 2.5,
      label == "Triphyophyllum"                  ~ img_size * 5,
      label == "Drosophyllum"                    ~ img_size * 5,
      label == "Non-carnivorous_caryophyllales"  ~ img_size * 1.2,
      TRUE ~ img_size
    )
  )

# ------------------------------------------------------
# BASE TREE PLOT
# ------------------------------------------------------
p <- ggtree(tree) %<+% traits +
  theme_tree() +
  
  geom_polygon(
    data = wedges,
    aes(x = x, y = y, group = tip),
    fill = "black",
    alpha = 1,
    inherit.aes = FALSE
  ) +
  
  geom_tippoint(aes(color = holocentric), size = 4) +
  scale_color_manual(values = c(none="grey60", mono="black", holo="red")) +
  
  geom_tiplab(aes(label = label2),
              size = 4,
              offset = 0.5,
              hjust = 0) +
  
  ggtitle("Lineage-specific species richness") +
  coord_cartesian(clip = "off") +
  xlim(NA, max(d$x) + 5)

# ------------------------------------------------------
# NO-OVERLAP IMAGE POSITIONING
# ------------------------------------------------------
TIP_LABEL_OFFSET    <- 0.5      # must match geom_tiplab(offset = ...)
CHAR_WIDTH_SCALING  <- 0.06     # tweak to match visual label length
LABEL_IMAGE_PADDING <- 2.5     # gap between label and image

traits2 <- traits2 %>%
  mutate(
    label_nchar   = nchar(label2),
    label_start_x = tip_x + TIP_LABEL_OFFSET,
    label_end_x   = label_start_x + CHAR_WIDTH_SCALING * label_nchar,
    image_x       = label_end_x + LABEL_IMAGE_PADDING
  )

# ------------------------------------------------------
# ADD IMAGES
# ------------------------------------------------------
p <- p +
  geom_image(
    data = traits2,
    aes(x = image_x, y = y, image = image, size = img_size),
    inherit.aes = FALSE
  ) +
  scale_size_identity()

# ------------------------------------------------------
# SAVE PLOT
# ------------------------------------------------------
ggsave("Drosera_caryophyllales.png", p,
       width = 10, height = 6, dpi = 300)

