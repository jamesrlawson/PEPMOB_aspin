require(gridExtra)

plot_panel <- function(df, split, x, y, z, ...) {

  # split is a vector of factors for df to be split by, can have up to three elements
  # x is the dependent variable
  # y is the independent variable
  # z is an interaction term  
  
  x <- deparse(substitute(x))
  y <- deparse(substitute(y))
  z <- deparse(substitute(z))
  
#regression <- regression_output(df, split, x, y, z)

df$species_time_comb <- paste(as.character(df[[split[1]]]), df[[split[2]]], as.character(df[[split[3]]]), sep = "-")

plt <-  ggplot(phys_prot, aes_string(x=x,y=y)) +
        geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
      #  geom_smooth(method=lm,se=FALSE) +
        geom_smooth(data=df[df$species_time_comb %in% regression[regression$pval < 0.05,]$submodel,], 
                   aes_string(x=x, y=y), method='lm', se=FALSE, size = 0.5, color = 'black') +
        scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
        facet_grid(. ~ time, scales = "free_x") +
        ylab('Heat Shock Protein (mg / m2)') +
        xlab('log10(Gs) (mmol H20 / s / m2)') 

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(regression, rows=NULL, theme=tt)
# Plot chart and table into one object
grid.arrange(plt, tbl,
             nrow=2,
             as.table=TRUE,
             heights=c(3,1))
}

regression <- regression_output(phys_prot, c('time'), HSPs, log10gs)
plot_panel(phys_prot, c('time'), HSPs, log10gs)




