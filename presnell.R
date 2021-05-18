main_run <- "results/SSLB/real/presnell/deseq_sf/raw/expressed/tensor/run_seed_8080_K_60"


sample_info <- read.csv("data/real/presnell/deseq_sf/raw/expressed/tensor/sample_info_tidied.txt",
                       sep='\t')
sample_info$age <- sample_info$Characteristics..age.

# We found that all samples with no sex given are in fact male, as can be told by
# the expression of genes on Y chromosome such as DDX3Y
sample_info$sex <- replace(sample_info$sex, sample_info$sex == "  ", "male")

sample_names <- sample_info[['sample_description']]
gene_info <- read.csv("data/real/presnell/deseq_sf/raw/expressed/tensor/gene_info.txt",
                      sep="\t",
                      stringsAsFactors = FALSE)
colnames(gene_info) <- gsub("ID..hsa\\d+.NAME..", "", colnames(gene_info))
colnames(gene_info) <- gsub("...Homo.sapiens..human..", "", colnames(gene_info))

# Gene symbols, but for the 28 genes which do not have a unique gene symbol (14 which are duplicated)
# add a suffix ___1 to the second such gene
gene_symbols <- make.unique(gene_info$GeneSymbol, sep="___")
rownames(gene_info) <- gene_symbols

Y <- read_matrix_from_folder("data/real/presnell/deseq_sf/raw/expressed/tensor",
                             "Y.txt")
colnames(Y) <- gene_symbols
rownames(Y) <- sample_names
