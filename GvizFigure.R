# Load required libraries
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# Import GTF file
gtf_file <- "Mus_musculus.GRCm38.102.gtf"
gtf <- import(gtf_file, format = "gtf")

# Filter out the gene with the name "A930002I21Rik"
filtered_gtf <- subset(gtf, gene_name != "A930002I21Rik")
gtf <- filtered_gtf

# Save the filtered data back to a GTF file (if needed)
# export(filtered_gtf, "Mus_musculus.GRCm38.102_Mus_SEEKR.gtf")

# Input gene symbols of interest
gene_symbol <- c("Kcna2", "Kcna10", "Slc6a17", "Kcna3", "Slc16a4", "Kcnc4", "AI504432")
genometype <- "mm39"
genetrackbuffer <- 32000
meta_merge <- TRUE  # Toggle exon merging

border_color = "transparent"
line_color = "gray"
highlight_color = "blue"


#####################################################


# Get only the genes of interest
selected_genes <- subset(gtf, sapply(tolower(gtf$gene_name), function(x) any(tolower(gene_symbol) == x)))

# Ensure we only use the chromosome of these genes
target_chr <- unique(as.character(seqnames(selected_genes)))

# Now filter the GTF file to only this chromosome
gtf <- subset(gtf, seqnames(gtf) %in% target_chr)

# Find the earliest start and latest end of the selected genes
earliest_start <- min(start(selected_genes))
latest_end <- max(end(selected_genes))

# Expand the range to include all genes within this region (+ buffer)
expanded_genes <- subset(gtf, start(gtf) <= latest_end + genetrackbuffer & end(gtf) >= earliest_start - genetrackbuffer)

# Ensure correct chromosome formatting
chr <- unique(as.character(seqnames(expanded_genes)))
ichr <- paste('chr', chr, sep = '')
seqlevels(expanded_genes) <- gsub(chr, ichr, seqlevels(expanded_genes))  # Adjust chromosome labels
gdf <- as.data.frame(expanded_genes)

# Translate Ensembl column names to Gviz column names
colnames(gdf)[colnames(gdf) == "seqnames"] <- "chromosome"
colnames(gdf)[colnames(gdf) == "transcript_id"] <- "transcript"
colnames(gdf)[colnames(gdf) == "transcript_biotype"] <- "feature"
colnames(gdf)[colnames(gdf) == "gene_id"] <- "gene"
colnames(gdf)[colnames(gdf) == "exon_id"] <- "exon"
colnames(gdf)[colnames(gdf) == "gene_name"] <- "symbol"

# Apply exon merging for each gene
merged_exons <- gdf %>%
  filter(type == "exon") %>%
  group_by(symbol) %>%
  arrange(start) %>%
  mutate(transcript = paste0(symbol, "_transcript")) %>%
  ungroup()

# Remove NA transcript values to prevent Gviz errors
merged_exons <- merged_exons[!is.na(merged_exons$transcript), ]

merged_exons$symbol[is.na(merged_exons$symbol)] <- merged_exons$transcript[is.na(merged_exons$symbol)]

merged_exons$exon[is.na(merged_exons$exon)] <- "NA"

# Add color and border columns based on gene symbol
merged_exons$color <- ifelse(merged_exons$symbol %in% gene_symbol, highlight_color, "gray")

# Assign unique group names to ensure per-gene color assignment
merged_exons$group <- merged_exons$symbol  # Ensures each gene is treated separately

# Split data by strand
pos_strand_data <- subset(merged_exons, strand == "+")
neg_strand_data <- subset(merged_exons, strand == "-")

# Ensure colors and borders are correctly mapped to transcript IDs
fill_colors <- setNames(as.character(merged_exons$color), merged_exons$transcript)

# Create an IdeogramTrack
itrack <- IdeogramTrack(genome = genometype, chromosome = ichr)

# Create a GenomeAxisTrack
genomeTrack <- GenomeAxisTrack(background.panel = "#FFFEDB",
                               fontsize = 10)

# Set the plotting region
start <- earliest_start - genetrackbuffer
end <- latest_end + genetrackbuffer

# Create GeneRegionTracks with explicit color mapping
pos_track <- GeneRegionTrack(
  pos_strand_data,
  genome = genometype,
  chromosome = ichr,
  name = "+",
  transcriptAnnotation = "name",
  fill = pos_strand_data$color,
  col = border_color, # scalar
  col.line = line_color, # scalar
  group = pos_strand_data$group,
  background.panel = "#E0FFFF" # Background color for pos_track
)

neg_track <- GeneRegionTrack(
  neg_strand_data,
  genome = genometype,
  chromosome = ichr,
  name = "-",
  transcriptAnnotation = "name",
  fill = neg_strand_data$color,
  col.line = line_color, # scalar
  col = border_color, # scalar
  group = neg_strand_data$group,
  background.panel = "#FFEBE0" # Background color for neg_track
)

# Save plot to PDF
pdf(paste0(paste(gene_symbol, collapse = "_"), "_gene_plot_auto_detected.pdf"), width = 8, height = 1.8) # Adjust height to fit all tracks
plotTracks(
  list(genomeTrack, pos_track, neg_track),
  from = start,
  to = end,
  collapseTranscripts = FALSE,
  collapse = FALSE,
  stacking = "squish",
  background.title = "darkblue",
  sizes = c(0.5,0.5,0.5)
)
dev.off()

