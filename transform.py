import sys

for line in sys.stdin:
    columns = line.strip('\n').split('\t')
    
    # Ensure we have all 8 columns
    if len(columns) < 8:
        continue
        
    # Extract needed columns based on the provided index order:
    # 4: chromosome_name, 5: strand, 6: external_gene_name, 7: transcription_start_site
    chromosome = columns[4]
    strand_val = columns[5]
    gene_name = columns[6]
    tss_str = columns[7]
    
    # Convert TSS to integer to perform math
    try:
        tss = int(tss_str)
    except ValueError:
        continue  # Skip rows with missing or malformed TSS data
        
    # Build the required BED columns
    col1 = chromosome
    col2 = str(tss)
    col3 = str(tss + 1)
    col4 = f"{col1}@{col2}-{col3}|{gene_name}"
    col5 = "."
    
    # Determine strand (+ if 1, - if -1)
    if strand_val == "1":
        col6 = "+"
    elif strand_val == "-1":
        col6 = "-"
    else:
        col6 = "." # Fallback for unexpected strand values

    # Store transformed columns in a list as requested
    transformed_row = [col1, col2, col3, col4, col5, col6]
    
    # Separate the contents by \t and print
    sys.stdout.write("\t".join(transformed_row) + "\n")