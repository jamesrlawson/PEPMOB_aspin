require(gridExtra)


regression_output <- function(df, split, x, y, z) { 
  
  # split is a vector of factors for df to be split by, can have up to three elements
  # x is the dependent variable
  # y is the independent variable
  # z is an interaction term
  
  x <- deparse(substitute(x))
  y <- deparse(substitute(y))
  z <- deparse(substitute(z))
  
  df$species_time_comb <- paste(as.character(df[[split[1]]]), df[[split[2]]], as.character(df[[split[3]]]), sep = "-")
  
  a <- data.frame()
  my.list <- vector("list", length(unique(df$species_time_comb)))
  
  for(i in 1:length(unique(df$species_time_comb))) {
    sub <- unique(df$species_time_comb)[i]
    df_sub <- subset(df, species_time_comb == sub)
    
    model <- lm(paste(y, '~', x, sep = ""), data = df_sub)
    
    model.stats <-  data.frame(cbind(sub,summary(model)$r.squared, summary(model)$coefficients[,4][2]))
    names(model.stats) <- c('submodel','R2','pval') 
    rownames(model.stats) <- unique(df$species_time_comb)[i]
    
    my.list[[i]] <- model.stats
    
  }
  
  a <- rbind(a, do.call(rbind, my.list))
  a$p.adj <- p.adjust(as.numeric(as.character(a$pval)), method = "BH")
  
  a$R2 <- signif(as.numeric(as.character(a$R2)), 3)
  a$pval <- signif(as.numeric(as.character(a$pval)), 5)
  a$p.adj <- signif(as.numeric(as.character(a$p.adj)), 5)
  a$submodel <- as.character(a$submodel)
  
  return(a)
  
}


plot_panel <- function(df, split, x, y, z, ...) {
  
  `+.uneval` <- function(a,b) {
    `class<-`(modifyList(a,b), "uneval")
  }
  
  # split is a vector of factors for df to be split by, can have up to three elements
  # x is the dependent variable
  # y is the independent variable
  # z is an interaction term  
  
  x <- deparse(substitute(x))
  y <- deparse(substitute(y))
  z <- deparse(substitute(z))
  
  split1 <- split[1]
  split2 <- split[2]
  split3 <- split[3]
  
  #regression <- regression_output(df, split, x, y, z)
  
  df$species_time_comb <- paste(as.character(df[[split[1]]]), df[[split[2]]], as.character(df[[split[3]]]), sep = "-")
  
  plt <-  ggplot(phys_prot, aes_string(x=x,y=y)) +
    ggtitle(paste(y, ' vs log10Gs', ' , submodels split by ', split1, split2, split3)) +
    geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
    #  geom_smooth(method=lm,se=FALSE) +
    geom_smooth(data=df[df$species_time_comb %in% regression[regression$pval < 0.05,]$submodel,], 
                aes_string(x=x, y=y), method='lm', se=FALSE, size = 0.5) +
    scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
    facet_grid(. ~ time, scales = "free_x") +
    ylab(paste(y,' (mg protein / g leaf dry mass', sep = "")) +
    xlab('log10(Gs) (mmol H20 / s / m2)') 
  
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)), base_size = 8)
  tbl <- tableGrob(regression, rows=NULL, theme=tt)
  table_size <- nrow(regression) * 0.3
  # Plot chart and table into one object
  grid.arrange(plt, tbl,
               nrow=2,
               as.table=TRUE,
               bottom=TRUE,
               heights=c(3,table_size))
}

write_panel <- function(df, split, x, y, z, ...) {
  
  `+.uneval` <- function(a,b) {
    `class<-`(modifyList(a,b), "uneval")
  }
  
  # split is a vector of factors for df to be split by, can have up to three elements
  # x is the dependent variable
  # y is the independent variable
  # z is an interaction term  
  
  x <- deparse(substitute(x))
  y <- deparse(substitute(y))
  z <- deparse(substitute(z))
  
  split1 <- split[1]
  split2 <- split[2]
  split3 <- split[3]
  
  #regression <- regression_output(df, split, x, y, z)
  
  df$species_time_comb <- paste(as.character(df[[split[1]]]), df[[split[2]]], as.character(df[[split[3]]]), sep = "-")
  
  plt <-  ggplot(phys_prot, aes_string(x=x,y=y)) +
    ggtitle(paste(y,'standardised by LMA', 'vs log10Gs', ' , submodels split by ', split1, split2, split3, sep = " ")) +
    geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
    #geom_smooth(method=lm,se=FALSE, size = 0.4, color = 'black
    geom_smooth(data=df[df$species_time_comb %in% regression[regression$pval < 0.05,]$submodel,], 
                aes_string(x=x, y=y), method='lm', se=FALSE, size = 0.7, color = 'black') +
    geom_smooth(data=df[df$species_time_comb %in% regression[regression$pval < 0.05,]$submodel,], 
                aes_string(x=x, y=y, color = paste(split3), linetype = paste(split2)), method='lm', se=FALSE, size = 0.5) +
    
    scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
    facet_grid(. ~ time, scales = "free_x") +
    ylab(paste(y,' (mg protein / g leaf dry mass)', sep = "")) +
    xlab('log10(Gs) (mmol H20 / s / m2)') 
  
  plt <- plt + theme_bw() 
  plt <- plt + theme_set(theme_bw(base_size = 8))
  plt <- plt + theme(
                 axis.title.y = element_text(hjust=0.35),
                 axis.title.x = element_text(vjust=0.35),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 axis.line = element_line(size=.5, color = "black"))
  
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)), base_size = 8)
  tbl <- tableGrob(regression, rows=NULL, theme=tt)
  table_size <- nrow(regression) * 0.2
  
  # output figure
  
  figname <- paste(y, 'standardised by LMA','vs', x, 'submodels split by', split1, split2, split3, Sys.Date(), sep = "_")
  pdf(sprintf("output/figs/%s.pdf", figname), width = 10, height = 10, pointsize=8)
  
  # Plot chart and table into one object
  grid.arrange(plt, tbl,
               nrow=2,
               as.table=TRUE,
               bottom=TRUE,
               heights=c(2,table_size))
  
  dev.off()
}



write_panel_ <- function(df, split, x, y, z, ...) {
  
  `+.uneval` <- function(a,b) {
    `class<-`(modifyList(a,b), "uneval")
  }
  
  # split is a vector of factors for df to be split by, can have up to three elements
  # x is the dependent variable
  # y is the independent variable
  # z is an interaction term  
  
  x <- deparse(substitute(x))
  y <- deparse(substitute(y))
  z <- deparse(substitute(z))
  
  split1 <- split[1]
  split2 <- split[2]
  split3 <- split[3]
  
  #regression <- regression_output(df, split, x, y, z)
  
  df$species_time_comb <- paste(as.character(df[[split[1]]]), df[[split[2]]], as.character(df[[split[3]]]), sep = "-")
  
  plt <-  ggplot(phys_prot, aes_string(x=x,y=y)) +
    ggtitle(paste(y, 'vs log10Gs', ' , submodels split by ', split1, split2, split3, sep = " ")) +
    geom_point(size = 2, aes(color=species, shape = temp_treatment)) +
    #geom_smooth(method=lm,se=FALSE, size = 0.4, color = 'black
    geom_smooth(data=df[df$species_time_comb %in% regression[regression$pval < 0.05,]$submodel,], 
                aes_string(x=x, y=y), method='lm', se=FALSE, size = 0.7, color = 'black') +
    geom_smooth(data=df[df$species_time_comb %in% regression[regression$pval < 0.05,]$submodel,], 
                aes_string(x=x, y=y, color = paste(split3), linetype = paste(split2)), method='lm', se=FALSE, size = 0.5) +
    
    scale_colour_hue(l=50) +  # Use a slightly darker palette than normal)
    facet_grid(. ~ time, scales = "free_x") +
    ylab(paste(y,' (mg / m2)', sep = "")) +
    xlab('log10(Gs) (mmol H20 / s / m2)') 
  
  plt <- plt + theme_bw() 
  plt <- plt + theme_set(theme_bw(base_size = 8))
  plt <- plt + theme(
    axis.title.y = element_text(hjust=0.35),
    axis.title.x = element_text(vjust=0.35),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(size=.5, color = "black"))
  
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)), base_size = 8)
  tbl <- tableGrob(regression, rows=NULL, theme=tt)
  table_size <- nrow(regression) * 0.2
  
  # output figure
  
  figname <- paste(y,'vs', x, 'submodels split by', split1, split2, split3, Sys.Date(), sep = "_")
  pdf(sprintf("output/figs/%s.pdf", figname), width = 10, height = 10, pointsize=8)
  
  # Plot chart and table into one object
  grid.arrange(plt, tbl,
               nrow=2,
               as.table=TRUE,
               bottom=TRUE,
               heights=c(2,table_size))
  
  dev.off()
}



