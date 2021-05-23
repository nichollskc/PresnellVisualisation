library(limma)
library(tibble)

read_matrix_from_folder <- function(folder, file) {
  full_filename <- paste(folder, file, sep="/")
  as.matrix(read.csv(full_filename, sep="\t", header=FALSE))
}

l2_norm <- function(x) {
  norm(x, type="2")
}

scale_factors_X_B_same_norm <- function(X, B) {
  # Scale so that columns of X and columns of B have same norm
  X_norms <- apply(X, MARGIN=2, FUN=l2_norm)
  B_norms <- apply(B, MARGIN=2, FUN=l2_norm)

  # Calculate the overall norm of each factor, then adjust the matrices so that for each
  # factor k X[, k] and B[, k] have the same L2 norm
  combined_norms <- X_norms * B_norms
  X_adjustment <- diag(sqrt(combined_norms) / X_norms)
  B_adjustment <- diag(sqrt(combined_norms) / B_norms)
  X_scaled <- X %*% X_adjustment
  B_scaled <- B %*% B_adjustment
  return(list("X_scaled"=X_scaled,
              "B_scaled"=B_scaled))
}

delete_below_threshold <- function(mat, threshold) {
  # Delete the elements that would be below the threshold if the matrix was scaled to have
  # columns with norm 1
  norms <- apply(mat, MARGIN=2, FUN=l2_norm)
  ifelse(abs(mat %*% diag(1 / norms)) > threshold, mat, 0)
}

read_thresholded_factor_matrices <- function(folder, threshold=0.01) {
  # Read in matrices
  X_raw <- read_matrix_from_folder(folder, "X.txt")
  B_raw <- read_matrix_from_folder(folder, "B.txt")

  # Scale so that columns of X and columns of B have same norm
  scaled <- scale_factors_X_B_same_norm(X_raw, B_raw)

  # Remove values smaller than threshold
  X_thresh <- delete_below_threshold(scaled$X_scaled, threshold)
  B_thresh <- delete_below_threshold(scaled$B_scaled, threshold)

  # Identify which factors have at least one gene and at least one sample
  X_nonempty <- (colSums(X_thresh != 0) != 0)
  B_nonempty <- (colSums(B_thresh != 0) != 0)
  both_nonempty <- X_nonempty & B_nonempty

  # Restrict to these non-empty factors
  X = X_thresh[, both_nonempty]
  B = B_thresh[, both_nonempty]

  return(list("X"=X,
              "B"=B,
              "K"=sum(both_nonempty)))
}

calc_factor_contribution <- function(biclustering, factor_index) {
  factor_contribution <- outer(biclustering$X[, factor_index],
        biclustering$B[, factor_index])
  return(factor_contribution)
}

calc_nz_factor_contribution <- function(biclustering, factor_index) {
  nz_samples <- (biclustering$X[, factor_index] != 0)
  nz_genes <- (biclustering$B[, factor_index] != 0)
  factor_contribution <- outer(biclustering$X[nz_samples, factor_index],
                               biclustering$B[nz_genes, factor_index])
  return(factor_contribution)
}

add_row_col_names <- function(biclustering, gene_symbols, sample_names) {
  rownames(biclustering$B) <- gene_symbols
  rownames(biclustering$X) <- sample_names
  colnames(biclustering$X) <- paste("factor", 1:biclustering$K, sep="_")
  colnames(biclustering$B) <- paste("factor", 1:biclustering$K, sep="_")
  
  return(biclustering)
}

generate_gene_info_with_fac <- function(biclustering, gene_info) {
  abs_B <- abs(biclustering$B)
  colnames(abs_B) <- paste0("abs_", colnames(abs_B))

  gene_info_with_fac <- cbind(gene_info, biclustering$B)
  gene_info_with_fac <- cbind(gene_info_with_fac, abs_B)
  return(gene_info_with_fac)
}

generate_sample_info_with_fac <- function(biclustering, sample_info) {
  abs_X <- abs(biclustering$X)
  colnames(abs_X) <- paste0("abs_", colnames(abs_X))

  sample_info_with_fac <- cbind(sample_info, biclustering$X)
  sample_info_with_fac <- cbind(sample_info_with_fac, abs_X)
  return(sample_info_with_fac)
}

calculate_pathway_enrichment <- function(biclustering, pathway_info) {
  total_rows <- biclustering$K * ncol(pathway_info)
  pathway_enrichment <- data.frame(PathwayName=rep("", total_rows),
                                   FactorIndex=rep(NA, total_rows),
                                   CameraPR_pvalue=rep(NA, total_rows),
                                   GeneSetTest_pvalue=rep(NA, total_rows))

  design <- model.matrix( ~ cell, sample_info)

  index = 1
  for (factor_index in 1:biclustering$K) {
    for (pathway_name in colnames(pathway_info)) {
      genes_in_pathway <- rownames(pathway_info %>%
                                   filter(!!as.name(pathway_name) != 0))
      genesettest_pvalue <- geneSetTest(genes_in_pathway,
                                        biclustering$B[, factor_index])
      camerapr_pvalue <- cameraPR(biclustering$B[, factor_index],
                                  genes_in_pathway)$PValue
      pathway_enrichment[index, ] = list(pathway_name,
                                         factor_index,
                                         camerapr_pvalue,
                                         genesettest_pvalue)
      index = index + 1
    }
  }
  
  pathway_enrichment$CameraPR_qvalue <- p.adjust(pathway_enrichment$CameraPR_pvalue,
                                                 method="BY")
  pathway_enrichment$GeneSetTest_qvalue <- p.adjust(pathway_enrichment$GeneSetTest_pvalue,
                                                    method="BY")

  return(pathway_enrichment)
}

add_odds_ratio_and_counts <- function(pathway_enrichment, biclustering, pathway_info) {
  pathway_enrichment$GenesInPathway <- rep(colSums(pathway_info), times=biclustering$K)
  pathway_enrichment$GenesInFactor <- rep(colSums(biclustering$B != 0), each=ncol(pathway_info))
  intersections = t(as.matrix(pathway_info)) %*% as.matrix(biclustering$B != 0)
  pathway_enrichment$GenesInIntersection <- as.numeric(intersections)
  factor_not_pathway <- pathway_enrichment$GenesInFactor - pathway_enrichment$GenesInIntersection
  odds_in_factor <- pathway_enrichment$GenesInIntersection / factor_not_pathway
  
  pathway_not_factor <- pathway_enrichment$GenesInPathway - pathway_enrichment$GenesInIntersection
  neither <- nrow(pathway_info) - pathway_enrichment$GenesInIntersection - factor_not_pathway - pathway_not_factor
  odds_out_factor <- pathway_not_factor / neither
  pathway_enrichment$OddsRatio <- odds_in_factor / odds_out_factor
  return(pathway_enrichment)
}

count_samples_by_type <- function(sample_info) {
  total_sample_counts_by_cell_disease <- sample_info %>%
    group_by(cell, short_disease) %>%
    summarize(total=n()) %>%
    pivot_wider(names_from=short_disease, values_from=total) %>%
    column_to_rownames(var='cell') %>%
    select(order(colnames(.)))
  
  total_sample_counts_by_sex_disease <- sample_info %>%
    group_by(sex, short_disease) %>%
    summarize(total=n()) %>%
    pivot_wider(names_from=short_disease, values_from=total) %>%
    column_to_rownames(var='sex') %>%
    select(order(colnames(.)))
  
  return(list("by_cell_disease"=total_sample_counts_by_cell_disease,
              "by_sex_disease"=total_sample_counts_by_sex_disease))
}

counts_by_variables <- function(sample_info_with_fac, groupvar_1, groupvar_2, factor_index=NULL) {
  if (! is.null(factor_index)) {
    sample_info_with_fac$included <- factor(sample_info_with_fac[[paste0("factor_", factor_index)]] != 0)  
  } else {
    sample_info_with_fac$included <- TRUE
  }
  
  counts <- sample_info_with_fac %>%
    # Group by the two variables we eventually want as rows/columns and also by "included"
    # We do .drop=F to keep all groups that exist in the dataset, even if not present in our factor
    group_by(!!as.name(groupvar_1), !!as.name(groupvar_2), included, .drop=FALSE) %>%
    # Count by group
    summarise(count=n()) %>%
    # Restrict to the columns we want and then pivot to go from "each row is a group" to
    # essentially a 2D grid, with groupvar_1 as rows, groupvar_2 as columns
    #   (plus extra columns we will get rid of soon)
    select(!!as.name(groupvar_1), !!as.name(groupvar_2), count, included) %>%
    pivot_wider(names_from=groupvar_2, values_from=count) %>%
    # Finally, restrict to rows that should be included and then we can discard this column
    filter(included == TRUE) %>%
    select(-included) %>%
    # Use groupvar_1 as rownames
    column_to_rownames(var=groupvar_1) %>%
    # Put column names in alpahbetical order
    select(order(colnames(.)))
  
  # Replace any NA with 0
  counts[is.na(counts)] <- 0
  
  return(counts)
}

proportions_to_integral_percentages <- function(df) {
  percentages <- data.frame(apply(df * 100, MARGIN=2, as.integer), row.names = rownames(df))
  
  text <- percentages %>% apply(MARGIN=2, FUN=as.character)
  text[is.na(percentages) | percentages == "0"] <- ""
  text[text != ""] <- paste0(text[text != ""], "%")
  
  percentages[is.na(percentages)] <- -100
  return(list("percentages" = percentages,
              "text" = text))
}

convert_percentage_text <- function(x) {
  print(class(x))
  print(x)
  x[is.na(x) | x == 0] <- ""
  x[x != ""] <- paste0(x[x != ""], "%")
  
  return(x)
}

generate_hovertext_sample_groups <- function(total_count, factor_count, cell, disease) {
  paste0("Sample type: ", cell, ", ", disease,
         "<br>Samples in dataset: ", total_count,
         "<br>Samples in factor: ", factor_count)
}

sample_heatmap <- function(total_counts, factor_counts, colour_palette, show_x_labels=TRUE) {
  factor_counts_long <- factor_counts %>%
    rownames_to_column("groupvar1") %>%
    pivot_longer(-groupvar1, names_to="groupvar2", values_to="factor_count")
  total_counts_long <- total_counts %>%
    rownames_to_column("groupvar1") %>%
    pivot_longer(-groupvar1, names_to="groupvar2", values_to="total_count")
  sample_groups_info <- full_join(total_counts_long, factor_counts_long, by=c("groupvar1", "groupvar2")) %>%
    mutate(proportion=factor_count/total_count,
           percentage=as.integer(proportion * 100),
           percentage_text=convert_percentage_text(percentage),
           SampleGroup=generate_hovertext_sample_groups(total_count, factor_count, groupvar1, groupvar2))

  hm <- ggplot(sample_groups_info,
               aes(x=groupvar2, y=groupvar1, fill=percentage, label=percentage_text, text=SampleGroup)) +
    geom_tile() +
    geom_text() +
    scale_fill_distiller(palette=colour_palette, direction=1,
                         na.value="grey50",
                         limits=c(0, 100),
                         guide=FALSE) +
    theme(panel.grid=element_blank(),
          axis.title=element_blank(),
          panel.background=element_blank())
  if (show_x_labels) {
    hm <- hm + theme(axis.text.x=element_text(angle=45, hjust = "right"))
  } else {
    hm <- hm + theme(axis.text.x=element_blank())
  }
  with_hover <- ggplotly(tooltip="SampleGroup") %>%
    layout(hoverlabel=list("font"=list("family"='sans-serif', "size"=25)))
      
  return(with_hover)
}
