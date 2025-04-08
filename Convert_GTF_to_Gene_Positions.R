library(data.table)
library(dplyr)

# Read GTF (adjust path)
gtf <- fread("gencode.v38.annotation.gtf.gz", skip = "chr")
genes <- gtf[gtf$V3 == "gene", ]

# Extract gene_id, chr, start, end
gene_positions <- genes %>%
  mutate(
    gene_id = gsub(".*gene_id \"([^\\\"]+)\".*", "\\1", V9),
    gene_name = gsub(".*gene_name \"([^\\\"]+)\".*", "\\1", V9)
  ) %>%
  select(chr = V1, start = V4, end = V5, gene_id, gene_name)

gene_order_clean <- gene_positions[, .(gene_name, chr, start, end)]
gene_order_clean <- gene_order_clean[!duplicated(gene_name)]

setcolorder(gene_order_clean, c("gene_name", "chr", "start", "end"))

#fwrite(gene_order_clean, "gene_order_infercnv.txt", sep = "\t", col.names = FALSE)

# Save as TXT
write.table(gene_order_clean, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/single_cell/sarcoma_before_after_radiation/test_scripts/gencode.v38.gene_positions.txt", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)

