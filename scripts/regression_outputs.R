
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
  
  a$R2 <- as.numeric(as.character(a$R2))
  a$pval <- as.numeric(as.character(a$pval))
  a$p.adj <- as.numeric(as.character(a$p.adj))
  a$submodel <- as.character(a$submodel)
  
  return(a)
  
}



