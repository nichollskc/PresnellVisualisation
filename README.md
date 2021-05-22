# Visualising Presnell biclustering results

# Todo

- Hover labels on sample heatmap e.g. show "No samples of this type" for NA values
- Also display plots and tables in individual tabs
- Correlation between genes in original dataset
- Correlation between genes restricted to the samples in the factor
- Add user guide box at top

# Caching

Caching data so that initial load up is fast, and reviewing plots you've already looked at is fast, is probably worth implementing

https://shiny.rstudio.com/articles/caching.html

# Colour scales

scale_fill_gradient2 is quite convenient - you just provide top and bottom colour and most usefully you can specify that 0 should be the midpoint.
heatmaply(fac_12, Rowv=F, Colv=FALSE, scale_fill_gradient_fun = scale_fill_gradient2(low="blue", high="red", midpoint=0, limits=c(-100000, 100000), trans="pseudo_log"))

scale_fill_distiller allows you to use palettes, which often look better and give a better scale than simple gradient
heatmaply(fac_12, Rowv=F, Colv=FALSE, scale_fill_gradient_fun = scale_fill_distiller(palette="RdBu", limits=c(-100000, 100000), trans="pseudo_log"))

# Custom spinner

I was really annoyed with how tiny the normal numeric input spinner buttons were and decided to switch to use the jquery spinner so I had more control. I am yet to find a way to dynamically update the maximum factor, so for now I'm settling with hard coding it in the www/js/numeric_spinner.js file.

I also think the way I've done it you can only have one such spinner in the document, which obviously isn't ideal in general but is fine for now.


# Sample type heatmaps

heatmaply doesn't offer enough flexibility with hovertext

Try ggplot with geom_tile and then ggplotly(p, tooltip=hovertext)
