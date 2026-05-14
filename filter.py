import sys

# Create a dictionary for O(1) lookup of desired chromosomes (1-22, X, Y)
# We assign 'True' as the value, though we only care about the keys.
desired_chromosomes = {str(i): True for i in range(1, 23)}
desired_chromosomes['X'] = True
desired_chromosomes['Y'] = True

for line in sys.stdin:
    # Skip the header row if it exists
    if line.startswith('ensembl_transcript_id'):
        continue
        
    # Split the line by tabs
    columns = line.strip('\n').split('\t')
    
    # Ensure the row has enough columns to avoid IndexError
    if len(columns) > 4:
        # chromosome_name is the 5th column (index 4)
        chromosome_name = columns[4]
        
        # Check if the chromosome is in our dictionary (avg O(1) lookup)
        if chromosome_name in desired_chromosomes:
            # sys.stdout.write is faster than print() and avoids adding extra newlines
            sys.stdout.write(line)