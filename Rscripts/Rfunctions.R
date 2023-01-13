## Functions for R
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