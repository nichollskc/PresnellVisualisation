main_run <- "/Users/kath/docs/PhD/biclustering/E-GEOD-60424/results/SSLB/real/presnell/deseq_sf/raw/expressed/tensor/run_seed_8080_K_60"

sample_info <- read.csv("/Users/kath/docs/PhD/biclustering/E-GEOD-60424/data/real/presnell/deseq_sf/raw/expressed/tensor/sample_info_tidied.txt",
                       sep='\t')

sample_names <- sample_info[['sample_description']]
gene_info <- read.csv("/Users/kath/docs/PhD/biclustering/E-GEOD-60424/data/real/presnell/deseq_sf/raw/expressed/tensor/gene_info.txt",
                      sep="\t",
                      stringsAsFactors = FALSE)
# Gene symbols, but for the 28 genes which do not have a unique gene symbol (14 which are duplicated)
# add a suffix ___1 to the second such gene
gene_symbols <- make.unique(gene_info$GeneSymbol, sep="___")
