library(ggplot2)

## Functions for R
# Create 2d plotting function
pca2d <- function(tab,d1,d2,perc,out){
  # Takes a table and its columns and does a plot
  png(file=paste(out,sprintf("pca_2d_d%s_d%s.png",d1,d2), sep = "/"))
  # ggplot(tab, aes(x=V1, y=V2, color=Treatment, size=20)) + geom_point() + theme_minimal()
  plot(tab[,d1], tab[,d2], type="n",
       xlab=sprintf("Dimension %s (%s %%)",d1,perc[d1]),
       ylab=sprintf("Dimension %s (%s %%)",d2,perc[d2]), main="")
  text(tab[,d1], tab[,d2], rownames(tab), cex=0.8)
  dev.off()
}

# Select the database for a list of organisms. Updateable.
select.organism <- function(organism){
  if (organism == "human") {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    database <- org.Hs.eg.db
  } else if (organism == "mouse") {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    database <- org.Mm.eg.db
  } else {
    print("Organism not recognized")
    quit()
  }
  return(database)
}

# Generate a table with the versions of the packages used.
get_versions <- function (outfile) {
  # Get packages versions
  pkglist <- sessionInfo()$otherPkgs
  
  pkgnames <- c()
  pkgvers <- c()
  
  pkgdf <- data.frame()
  for (i in 1:length(pkglist)) {
    pkgnames[i] <-names(sessionInfo()$otherPkgs[i])
    pkgvers[i] <- sessionInfo()$otherPkgs[[i]]$Version
  }
  
  pkgdf <- data.frame(list(pkgnames, pkgvers))
  colnames(pkgdf) <- c("Name","Version")
  write.table(pkgdf, file=outfile, sep="\t", row.names = FALSE)
}