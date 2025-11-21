library(ggtree)
library(treeio)
library(ggplot2)

tree <- read.tree("Drosera_genus.tree")

# Build base tree (italic tip labels)
p <- ggtree(tree, layout="rectangular") +
  geom_tiplab(aes(label = paste0("italic(", label, ")")),
              parse = TRUE,
              hjust=-0.05, size=4)

df <- p$data
max_x <- max(df$x)

# Extra room for brackets + labels
p <- p + xlim(0, max_x + 10)   # a little more room overall
df <- p$data


# ============================
# Character state (holocentric / monocentric)
# ============================

holos <- c("D.regia", "D.paradoxa", "D.scorpioides","D.roseana")
df$chromotype <- ifelse(df$label %in% holos & df$isTip,
                        "holocentric", "monocentric")
df$chromotype[!df$isTip] <- NA

p <- p +
  geom_point(data=df[df$isTip,],
             aes(x=x, y=y, color=chromotype),
             size=3) +
  scale_color_manual(values=c(
    "holocentric"="firebrick3",
    "monocentric"="gray20"
  ),
  name="Chromosome type",
  breaks=c("holocentric","monocentric"),
  labels=c("Holocentric","Monocentric"))


# ============================
# Subgenus brackets
# ============================

get_y_range <- function(df, tips){
  ys <- df$y[df$label %in% tips & df$isTip]
  c(min(ys), max(ys))
}

regiae_r    <- get_y_range(df, "D.regia")
ergaleium_r <- get_y_range(df, c("D.binata","D.paradoxa","D.roseana","D.scorpioides"))
drosera_r   <- get_y_range(df, c("D.filiformis","D.tokaiensis","D.capensis","D.aliciae"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MOVE BRACKETS TO THE RIGHT HERE
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
x_bracket <- max(df$x) + 3.5     # was 2.5
x_label   <- max(df$x) + 5.5     # ensure label is further right

square_bracket <- function(y1, y2, col) {
  p <<- p + geom_segment(aes(x=x_bracket, xend=x_bracket+1.2,
                             y=y2, yend=y2),
                         linewidth=1.3, color=col)
  p <<- p + geom_segment(aes(x=x_bracket+1.2, xend=x_bracket+1.2,
                             y=y1, yend=y2),
                         linewidth=1.3, color=col)
  p <<- p + geom_segment(aes(x=x_bracket, xend=x_bracket+1.2,
                             y=y1, yend=y1),
                         linewidth=1.3, color=col)
}

# Draw brackets
square_bracket(regiae_r[1],    regiae_r[2],    "#8B0000")
square_bracket(ergaleium_r[1], ergaleium_r[2], "#E68A00")
square_bracket(drosera_r[1],   drosera_r[2],   "#66BB00")


# Subgenus labels
p <- p +
  annotate("text", x=x_label, y=mean(regiae_r),
           label="Subg. Regiae", color="#8B0000", size=5, hjust=0) +
  annotate("text", x=x_label, y=mean(ergaleium_r),
           label="Subg. Ergaleium", color="#E68A00", size=5, hjust=0) +
  annotate("text", x=x_label, y=mean(drosera_r),
           label="Subg. Drosera", color="#66BB00", size=5, hjust=0)


# Save
ggsave("Drosera_genus.png", p, width=8, height=6)

