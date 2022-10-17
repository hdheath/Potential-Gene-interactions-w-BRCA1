## LAB 3 - BIT 150

library(tidyverse)

# View the files of the directory we are working in :
#list.files()

# Open a dataset and convert it to tabuklar format :
brca1_data <- read_tsv('BRCA1_String.tsv')

# View the data in the newly created tsv
#view(brca1_data)

# Change the name of the first column to not include # :
colnames(brca1_data)[1] <- 'node1'
colnames(brca1_data)

# Filter the table for rows in column 'node1' that contain 'BRCA1'
filter(brca1_data,node1 == 'BRCA1')

# Filter the table for rows in column "exp' for values above .5
filter(brca1_data,experimentally_determined_interaction > 0.5)

# First sort
brca1_data_temp <- arrange(brca1_data,desc(experimentally_determined_interaction))

# Then filter. Note we use the temporary file as the argument to filter now
brca1_data_temp <- filter(brca1_data_temp,node1 == 'BRCA1')

# Now view
#View(brca1_data_temp)

# Create a new table with one row
new_interaction <- tibble(node1='SEPT2', node2='MARCH1')

# View it
#new_interaction

# Bind it to the bottom of the old table
brca1_data_temp <- bind_rows(brca1_data_temp,new_interaction)

# View the bottom of the new table
tail(brca1_data_temp)

# Save new table to tsv
write_tsv(brca1_data_temp,file = 'processed_BRCA1_interactions.tsv')

#         ACTIVITY        #

# Filter gene interaction data for BRCA
brca1_data_subset <- filter(brca1_data,node1 == 'BRCA1')
#View(brca1_data_subset)
brca1_data_subset <- filter(brca1_data_subset,experimentally_determined_interaction > .5)

# Turn filtered dataset to a tsv 
write_tsv(brca1_data_subset,file = 'brca1_data_subset.tsv')
#view(brca1_data_subset)

# Load expression and ID mapping data :
tissue_expression=read_tsv('normal_tissue.tsv')
gene_to_protein=read_tsv('Gene_to_Protein.tsv')

# Make column names easier to access :
colnames(tissue_expression)[2] = 'Gene_name'
colnames(tissue_expression)[4] = 'Cell_type'
#glimpse(tissue_expression)

colnames(gene_to_protein)[1] = 'Gene_stable_ID'
colnames(gene_to_protein)[2] = 'Protein_stable_ID'
#glimpse(gene_to_protein)

# remove the data that does not include a protein ID : 
gene_to_protein = filter(gene_to_protein, !is.na(Protein_stable_ID))

# Add a new column to the brca1_data table for the ENSEMBL Gene IDs of the node2 genes :
brca1_data_subset <- separate(brca1_data_subset,col='node2_string_id',into=c("node2_ID","ENSP_2"))
#glimpse(brca1_data_subset)

# check that each of the Protein ID's from the brca1_data_subset table is present in the gene_to_protein table:
brca1_data_subset$ENSP_2 %in% gene_to_protein$Protein_stable_ID

#check is that our desired key for the new data (Protein_stable_ID in the gene_to_protein table) is a valid key, that all values are unique
anyDuplicated(gene_to_protein$Protein_stable_ID)

gene_to_protein = unique(gene_to_protein)
anyDuplicated(gene_to_protein$Protein_stable_ID)

# Join the brca1_data_subset and gene_to_protein tables based on the node1 protein
brca1_data_subset = left_join(brca1_data_subset,gene_to_protein,by=c('ENSP_2' = 'Protein_stable_ID'))

#           Analysis          #
# First, find the tissues of expression for BRCA1
# Find the ensembl GENE ID :
BRCA1_ENSP = "ENSP00000418960"
filter(gene_to_protein,Protein_stable_ID == BRCA1_ENSP)

TP53_ENSP = "ENSG00000141510"
filter(gene_to_protein,Protein_stable_ID == TP53_ENSP)

# Equate the IDs to match each other
BRCA1_ENSG = "ENSG00000012048"

# Filter the data for the correct ID to BRCA
brca1_expression = filter(tissue_expression, Gene == BRCA1_ENSG)

# Refine this for only tissues with medium or high expression:
brca1_expression = filter(brca1_expression,Level %in% c('Medium','High'))
#view(brca1_expression)

# Extract for values with Medium or High expression
brca1_expression$Tissue

# Now, pull out the ID of the first interacting gene, and find its medium or high-expressed tissues:
node2_gene = brca1_data_subset$Gene_stable_ID[1]

node2_expression = filter(tissue_expression, Gene == node2_gene & Level %in% c("High","Medium"))
node2_expression$Tissue

common_tissues = intersect(brca1_expression$Tissue,node2_expression$Tissue)
#View(common_tissues)
length(common_tissues)

summary_table = tibble() # this creates an empty table

for(row in 1:nrow(brca1_data_subset)){  
  print(row)
  
  # extract the Gene_stable_ID for the node2 gene
  node2_gene = brca1_data_subset$Gene_stable_ID[row]
  
  # use filter to get the gene expression data for each protein
  #node2_expression = filter(tissue_expression, Gene == node2_gene & Level %in% c("High","Medium"))
  node2_expression = filter(tissue_expression, Gene == node2_gene & Level == 'High')
  
  # make a list of the tissues in common
  common_tissues = intersect(brca1_expression$Tissue,node2_expression$Tissue)
  
  # count the number of tissues in common
  num_common_tissues = length(common_tissues)
  
  # form the list of tissues into a single String so we can store the result in a cell of the table
  tissue_list = toString(common_tissues)
  
  # build a new one-row table with all these results
  interaction_table <- tibble(
    node1_name = 'BRCA1', 
    node2_name = brca1_data_subset$node2[row],# this pulls out the actual gene name
    num_common_tissues = num_common_tissues,
    tissue_list = tissue_list
  ) 
  
  # append this disease's results to the full table
  summary_table <- bind_rows(summary_table, interaction_table)    
  
}


View(summary_table)
hist(summary_table$num_common_tissues)
write_tsv(summary_table,file = 'BRCA1_interactors_common_expression.tsv')
filter(summary_table,num_common_tissues == 0)
