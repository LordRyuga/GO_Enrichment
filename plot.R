setwd("D:/Books_and_asses/cse_courses/comp bio/compbioass/gene_onotology_analysis")

# 1. Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Install required Bioconductor packages (only runs if missing)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")

# 3. Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# 4. Read the target gene list generated from the terminal pipeline
# (Make sure your R working directory contains this file)
genes <- read.table("hg38_GCGC_target_genes.tsv", 
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1

# 5. Perform Gene Ontology (GO) Enrichment Analysis
ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL", # Because we extracted gene names (symbols)
                ont           = "BP",     # Biological Process
                pAdjustMethod = "BH",     # Benjamini-Hochberg correction
                qvalueCutoff  = 0.01)

# 6. Save the raw data results to a CSV
write.csv(as.data.frame(ego), "hg38_GCGC_GO_results.csv")

# 7. Generate and save the Dotplot
pdf("hg38_GCGC_GO_dotplot.pdf", height = 8, width = 8)
dotplot(ego, showCategory = 20, font.size = 6)
dev.off()

# Print a success message to the console
print("GO Enrichment Analysis Complete! Check your folder for the PDF plot.")
