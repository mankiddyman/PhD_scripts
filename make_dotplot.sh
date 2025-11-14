#!/bin/bash
# Working directory
wd=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/self_dotplot
mkdir -p "$wd" && cd "$wd"
source ~/.bashrc
micromamba activate r_lang

# Genome path
GENOME=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/Files/genome/out_JBAT.FINAL.fa
GENOME16=genome16.fa
PAF=genome16_self.paf
SVG=genome16_self.svg

# Load modules
module load seqkit

# Step 1: Extract first 16 sequences
echo "Extracting first 16 sequences..."
seqkit head -n 16 "$GENOME" > "$GENOME16"

# Step 2: Run minimap2 self-alignment
echo "Running minimap2 self-alignment..."
minimap2 -x asm5 -DP -t 20 "$GENOME16" "$GENOME16" > "$PAF"

echo "Minimap2 alignment complete: $PAF"
echo "PAF file size: $(du -h $PAF)"

# Step 3: Run inline Python to generate SVG dotplot
echo "Generating SVG dotplot..."
python3 - << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# Input/Output
paf_file = "genome16_self.paf"
output_svg = "genome16_self.svg"

# --- Parameters ---
bin_size = 1_000  # Bin coordinates to reduce overplotting (1 kb)
line_alpha = 0.3
line_width = 0.5
fig_size = (30, 30)  # Large figure
dpi = 300            # High-res for raster export (optional)

# --- Read PAF ---
df = pd.read_csv(
    paf_file,
    sep='\t',
    header=None,
    usecols=range(12)
)
df.columns = ['qname', 'qlen', 'qstart', 'qend', 'strand',
              'tname', 'tlen', 'tstart', 'tend', 'nmatch',
              'alen', 'mapq']

print(f"Loaded {len(df)} alignments")

# --- Optional: Bin coordinates to reduce file size ---
df['qstart_bin'] = (df['qstart'] // bin_size) * bin_size
df['qend_bin']   = (df['qend']   // bin_size) * bin_size
df['tstart_bin'] = (df['tstart'] // bin_size) * bin_size
df['tend_bin']   = (df['tend']   // bin_size) * bin_size

# --- Prepare line segments for LineCollection ---
segments_plus = [((row.qstart_bin, row.tstart_bin), (row.qend_bin, row.tend_bin))
                 for _, row in df[df.strand=='+'].iterrows()]
segments_minus = [((row.qstart_bin, row.tstart_bin), (row.qend_bin, row.tend_bin))
                  for _, row in df[df.strand=='-'].iterrows()]

# --- Plot ---
fig, ax = plt.subplots(figsize=fig_size)

lc_plus = LineCollection(segments_plus, colors='blue', linewidths=line_width, alpha=line_alpha)
lc_minus = LineCollection(segments_minus, colors='red', linewidths=line_width, alpha=line_alpha)

ax.add_collection(lc_plus)
ax.add_collection(lc_minus)

# --- Axes and labels ---
ax.set_xlabel("Query (bp)", fontsize=18)
ax.set_ylabel("Target (bp)", fontsize=18)
ax.set_title("Self-dotplot (first 16 sequences)", fontsize=20)
ax.grid(True, alpha=0.2)

# Set limits
ax.set_xlim(0, df['qlen'].max())
ax.set_ylim(0, df['tlen'].max())

# --- Save ---
plt.tight_layout()
plt.savefig(output_svg, format='svg')  # Zoomable
plt.close()

print(f"Dotplot saved as {output_svg}")
EOF

echo "Dotplot generation complete: $SVG"
