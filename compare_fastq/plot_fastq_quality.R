library(jsonlite)
library(tidyverse)

file1_json <- "file1.json"
file2_json <- "file2.json"

j1 <- fromJSON(file1_json)
j2 <- fromJSON(file2_json)

# ============================================================
# Extract per-base mean quality
# ============================================================

# For your JSON:
# j$read1_before_filtering$quality_curves$mean is a large numeric vector

pb1 <- tibble(
  position = seq_along(j1$read1_before_filtering$quality_curves$mean),
  mean_quality = j1$read1_before_filtering$quality_curves$mean,
  sample = "File 1"
)

pb2 <- tibble(
  position = seq_along(j2$read1_before_filtering$quality_curves$mean),
  mean_quality = j2$read1_before_filtering$quality_curves$mean,
  sample = "File 2"
)

df_pb <- bind_rows(pb1, pb2)

# ============================================================
# Plot per-base mean quality
# ============================================================

p1 <- ggplot(df_pb, aes(x = position, y = mean_quality, color = sample)) +
  geom_line(linewidth = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Per-Base Mean Quality",
    x = "Read Position",
    y = "Mean Phred Score"
  )

ggsave("quality_per_base.png", p1, width = 10, height = 6, dpi = 200)


# ============================================================
# Extract per-base quality by nucleotide (A, C, G, T)
# ============================================================

extract_nuc_df <- function(j, label) {
  qc <- j$read1_before_filtering$quality_curves
  tibble(
    position = seq_along(qc$A),
    A = qc$A,
    C = qc$C,
    G = qc$G,
    T = qc$T,
    sample = label
  ) |>
  pivot_longer(c(A, C, G, T), names_to="nucleotide", values_to="quality")
}

df_nuc <- bind_rows(
  extract_nuc_df(j1, "File 1"),
  extract_nuc_df(j2, "File 2")
)

# ============================================================
# Plot per-base quality by nucleotide
# ============================================================

p2 <- ggplot(df_nuc, aes(position, quality, color=nucleotide)) +
  geom_line(alpha=0.8) +
  facet_wrap(~ sample, ncol = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Per-base Quality by Nucleotide",
    x = "Read Position",
    y = "Quality"
  )

ggsave("quality_per_base_by_nucleotide.png", p2, width = 10, height = 10, dpi = 200)


print("Generated plots:")
print(" - quality_per_base.png")
print(" - quality_per_base_by_nucleotide.png")

