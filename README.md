# Visualising Presnell biclustering results

# Todo

- Hover labels on sample heatmap e.g. show "No samples of this type" for NA values
- Factor selector needs to be easier to click
- Also display plots and tables in individual tabs
- Correlation between genes in original dataset
- Correlation between genes restricted to the samples in the factor
- Add user guide box at top
- Alternative scale for factor contribution - log scale?

# Caching

Caching data so that initial load up is fast, and reviewing plots you've already looked at is fast, is probably worth implementing

https://shiny.rstudio.com/articles/caching.html

# Colour scales

scale_fill_gradient2 is quite convenient - you just provide top and bottom colour and most usefully you can specify that 0 should be the midpoint.
heatmaply(fac_12, Rowv=F, Colv=FALSE, scale_fill_gradient_fun = scale_fill_gradient2(low="blue", high="red", midpoint=0, limits=c(-100000, 100000), trans="pseudo_log"))

scale_fill_distiller allows you to use palettes, which often look better and give a better scale than simple gradient
heatmaply(fac_12, Rowv=F, Colv=FALSE, scale_fill_gradient_fun = scale_fill_distiller(palette="RdBu", limits=c(-100000, 100000), trans="pseudo_log"))
