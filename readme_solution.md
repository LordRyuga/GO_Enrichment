# Gene Ontology Enrichment Analysis Pipeline: From Annotations to Dotplot

## Overview
This pipeline extracts transcription start sites (TSS) from human gene annotation data, identifies promoter regions by extending these coordinates upstream, extracts the DNA sequences for these promoters, and searches for a specific GC-rich transcription factor binding motif (`GCGC..GCGC`). Finally, it extracts the genes containing this motif and performs a Gene Ontology (GO) enrichment analysis to visualize their biological functions.

## Prerequisites
*   **Python 3.x**
*   **bedtools** (for genomic coordinate manipulation and sequence extraction)
*   **EMBOSS suite (`dreg`)** (for regex-based motif searching)
*   **R & RStudio** with `clusterProfiler` and `org.Hs.eg.db` packages installed
*   **Required Reference Files:**
    *   `human_gene_annotation.tsv` (Ensembl human gene annotations)
    *   `hg38.chrom.sizes` (UCSC chromosome sizes for hg38)
    *   `hg38.fa` (Reference human genome sequence)

---

## Step 1: Extracting Transcription Start Sites (TSS)
**Objective:** Parse the raw `.tsv` annotation file to extract desired chromosomes and reformat the data into a standard `.bed` file containing TSS coordinates. We use terminal standard input (`stdin`) to process the file line-by-line, maintaining an **O(1)** space complexity to minimize memory usage.

### 1A: Filtering for Valid Chromosomes
We wrote `filter.py` to keep only standard chromosomes (1-22, X, Y) using an O(1) dictionary lookup.

```python
# filter.py
import sys

desired_chromosomes = {str(i): True for i in range(1, 23)}
desired_chromosomes['X'] = True
desired_chromosomes['Y'] = True

for line in sys.stdin:
    if line.startswith('ensembl_transcript_id'):
        continue
        
    columns = line.strip('\n').split('\t')
    if len(columns) > 4:
        chromosome_name = columns[4]
        if chromosome_name in desired_chromosomes:
            sys.stdout.write(line)
```

**Command:**
```bash
cat human_gene_annotation.tsv | python filter.py > filtered.tsv
```

### 1B: Transforming to BED Format
We wrote `transform.py` to reorder the columns into a standard 6-column BED format. The 4th column is formatted as `Chromosome@Start-End|geneName` to carry the gene symbol forward.

```python
# transform.py
import sys

for line in sys.stdin:
    columns = line.strip('\n').split('\t')
    if len(columns) < 8:
        continue
        
    chromosome = columns[4]
    strand_val = columns[5]
    gene_name = columns[6]
    tss_str = columns[7]
    
    try:
        tss = int(tss_str)
    except ValueError:
        continue
        
    col1 = chromosome
    col2 = str(tss)
    col3 = str(tss + 1)
    col4 = f"{col1}@{col2}-{col3}|{gene_name}"
    col5 = "."
    
    if strand_val == "1":
        col6 = "+"
    elif strand_val == "-1":
        col6 = "-"
    else:
        col6 = "."

    transformed_row = [col1, col2, col3, col4, col5, col6]
    sys.stdout.write("\t".join(transformed_row) + "\n")
```

**Command:**
```bash
cat filtered.tsv | python transform.py > hg38_tss.bed
```

---

## Step 2: Defining Promoter Regions
**Objective:** Extend the TSS coordinates 500 base pairs upstream in a strand-aware manner to capture the core promoter regions where transcription factors bind. 

### 2A: Harmonizing Chromosome Names (Ensembl vs. UCSC)
The `.bed` file uses Ensembl naming (e.g., `1`), while the `.genome` sizes file uses UCSC naming (`chr1`). We append the "chr" prefix using `awk`.

**Command:**
```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {$1 = "chr"$1; print $0}' hg38_tss.bed > hg38_tss_chr.bed
```

### 2B: Extending Coordinates with `bedtools slop`
We use `bedtools slop` to extend the coordinates 500bp upstream (`-l 500`), 0bp downstream (`-r 0`), while respecting the strand direction (`-s`) and preventing out-of-bounds errors using the chromosome sizes file (`-g`).

**Command:**
```bash
bedtools slop -i hg38_tss_chr.bed -g hg38.chrom.sizes -l 500 -r 0 -s > hg38_promoters_500bp.bed
```

---

## Step 3: Extracting Promoter Sequences
**Objective:** Fetch the actual DNA sequences for our newly defined 500bp promoter regions. The `-name` flag ensures our custom column 4 (`Chromosome@Start-End|geneName`) is used as the FASTA header.

**Command:**
```bash
bedtools getfasta -fi hg38.fa -bed hg38_promoters_500bp.bed -fo hg38_promoters_500bp.fa -name
```

---

## Step 4: Motif Searching
**Objective:** Search the promoter sequences for the specific transcription factor binding motif `GCGC..GCGC` using the EMBOSS `dreg` tool.

**Command:**
```bash
dreg -sequence hg38_promoters_500bp.fa -pattern "GCGC..GCGC" -outfile hg38_promoters_500bp_GCGC_hits.txt
```

---

## Step 5: Extracting the Target Gene List
**Objective:** The `dreg` tool strips away sequence names, leaving only coordinates. We must extract the coordinates of positive hits and map them back to our BED file to recover the original gene symbols.

### 5A: Extract Positive Hit Coordinates
Extract the coordinates for any sequence with a `HitCount > 0`.

**Command:**
```bash
awk '/Sequence:/ {seq=$3} /HitCount:/ {if ($3 > 0) print seq}' hg38_promoters_500bp_GCGC_hits.txt > hit_coords.txt
```

### 5B: Map Back to BED and Extract Gene Symbols
Match these coordinates against columns 2 and 3 of the BED file, extract the 4th column, parse out the gene symbol using the pipe delimiter (`|`), and remove duplicates.

**Command:**
```bash
awk 'NR==FNR {hits[$1]=1; next} {coord=$2"-"$3; if (coord in hits) print $4}' hit_coords.txt hg38_promoters_500bp.bed | cut -d '|' -f 2 | sort -u > hg38_GCGC_target_genes.tsv
```

---

## Step 6: Gene Ontology (GO) Enrichment Analysis
**Objective:** Use R to perform a biological process enrichment analysis on the final list of target genes and generate a dotplot visualizing the results.

**R Script (`go_analysis.R`):**
```R
# 1. Install required packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")

# 2. Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# 3. Read the target gene list
genes <- read.table("hg38_GCGC_target_genes.tsv", 
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1

# 4. Perform GO Enrichment Analysis (Biological Process)
ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL", 
                ont           = "BP",     
                pAdjustMethod = "BH",     
                qvalueCutoff  = 0.01)

# 5. Save results and generate plot
write.csv(as.data.frame(ego), "hg38_GCGC_GO_results.csv")

pdf("hg38_GCGC_GO_dotplot.pdf", height = 8, width = 8)
dotplot(ego, showCategory = 20, font.size = 6)
dev.off()
```

**Final Outputs:** 
*   `hg38_GCGC_GO_results.csv`: Raw tabular data of enriched biological processes.
*   `hg38_GCGC_GO_dotplot.pdf`: A visual representation of the top 20 enriched biological processes driven by genes with the `GCGC..GCGC` motif in their promoters.