# FULL DATA WITH NO TPM FILTERING threshold = 30 000  pval = 0.000267
upreg_ERVs_gene_close <- 55  # Number in Category A
upreg_ERVs_no_gene_close <- 47   # Number in Not Category A
downreg_ERVs_gene_close <- 6  # Number in Category A
downreg_ERVs_no_gene_close <- 28

# Create a matrix with the data
contingency_table <- matrix(c(upreg_ERVs_gene_close, downreg_ERVs_gene_close, upreg_ERVs_no_gene_close, downreg_ERVs_no_gene_close), nrow = 2, byrow = TRUE)

# Add row and column names
rownames(contingency_table) <- c("Upregulated ERV with upregulated gene nearby", "Downregulated ERV with upregulated gene nearby")
colnames(contingency_table) <- c("Upregulated ERV with no upregulated gene nearby", "Downregulated ERV with no upregulated gene nearby")

# Print the contingency table
print(contingency_table)
fisher.test(contingency_table)


#DATA WITH TPM > 10 - this doesnt work threshold  = 30 000 p-val = 1 
upreg_ERVs_gene_close <- 5  # Number in Category A
upreg_ERVs_no_gene_close <- 8   # Number in Not Category A
downreg_ERVs_gene_close <- 2  # Number in Category A
downreg_ERVs_no_gene_close <- 5 

# Create a matrix with the data
contingency_table <- matrix(c(upreg_ERVs_gene_close, downreg_ERVs_gene_close, upreg_ERVs_no_gene_close, downreg_ERVs_no_gene_close), nrow = 2, byrow = TRUE)

# Add row and column names
rownames(contingency_table) <- c("Upregulated ERV with upregulated gene nearby", "Downregulated ERV with upregulated gene nearby")
colnames(contingency_table) <- c("Upregulated ERV with no upregulated gene nearby", "Downregulated ERV with no upregulated gene nearby")

# Print the contingency table
print(contingency_table)
fisher.test(contingency_table)


#DATA OF LOW EXPRESSED ERV'S TPM < 5 10 - threshold = 30 000, pval = 0.0001442
upreg_ERVs_gene_close <- 50   # Number in Category A
upreg_ERVs_no_gene_close <- 39   # Number in Not Category A
downreg_ERVs_gene_close <- 4  # Number in Category A
downreg_ERVs_no_gene_close <- 23 


#TPM
sum_lowly_expreseed_ERVs <- 50 + 39 + 4 + 23
sum_lowly_expreseed_ERVs

# Create a matrix with the data
contingency_table <- matrix(c(upreg_ERVs_gene_close, downreg_ERVs_gene_close, upreg_ERVs_no_gene_close, downreg_ERVs_no_gene_close), nrow = 2, byrow = TRUE)

# Add row and column names
rownames(contingency_table) <- c("Upregulated ERV with upregulated gene nearby", "Downregulated ERV with upregulated gene nearby")
colnames(contingency_table) <- c("Upregulated ERV with no upregulated gene nearby", "Downregulated ERV with no upregulated gene nearby")

# Print the contingency table
print(contingency_table)
fisher.test(contingency_table)


#CHECK IF ERVGROUP IS OVERREPRESSENTED IN THE OVEREXPRESSED ERVS 
#DATA WITH TPM > 10 - this doesnt work threshold  = 30 000 p-val = 1 
diff_reg_ervs_total <- 136  # Number in Category A
diff_reg_ervs_tpm10 <- 20   # Number in Not Category A
herv_ervs_total <- 26  # Number in Category A
herv_ervs_tpm10 <- 8 

diff_reg_ervs_total <- 90  # Number in Category A
diff_reg_ervs_tpm10 <- 12   # Number in Not Category A
herv_ervs_total <- 16  # Number in Category A
herv_ervs_tpm10 <- 8 


diff_reg_ervs_total <- 136  # Number in Category A
diff_reg_ervs_tpm10 <- 9   # Number in Not Category A
herv_ervs_total <- 7# Number in Category A
herv_ervs_tpm10 <- 4

# Create a matrix with the data
contingency_table <- matrix(c(diff_reg_ervs_total, diff_reg_ervs_tpm10, herv_ervs_total, herv_ervs_tpm10), nrow = 2, byrow = TRUE)

# Add row and column names
rownames(contingency_table) <- c("Upregulated ERV with upregulated gene nearby", "Downregulated ERV with upregulated gene nearby")
colnames(contingency_table) <- c("Upregulated ERV with no upregulated gene nearby", "Downregulated ERV with no upregulated gene nearby")

# Print the contingency table
print(contingency_table)
fisher.test(contingency_table)





