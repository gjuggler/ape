# tree_io.R - a series of round-trip tests for ape's phylogenetic tree I/O.
#
# Author: Gregory Jordan (greg@ebi.ac.uk)
#
# Instructions: use the following command from within the ape base directory:
# > R --vanilla < test/tree_io.R

library(ape)

source("R/read.nhx.R")
source("R/write.tree.R")

round.trip <- function(str.in) {
  phylo <- read.nhx(str.in)

  if (!is.null(phylo$.tags)) {
    print(str(phylo))
  }

  str.out <- write.tree(phylo)
  if ((!identical(str.in, str.out))) {
    print(paste('Round trip failed!', "Expected [", str.in, "] Got [", str.out,"]", sep=''))
  }
}
round.trip("(a,(b,c));")
round.trip("(a,b,c,d,e,f,g);")
round.trip("(a,(b,c,d));")
round.trip("(a,(b,(c,d)e)f)g;")
round.trip("(a:1,(b:1,c:1):1):1;")
round.trip("(a,b);")
round.trip("(a,(b,c));")
round.trip("(a,(b));")
round.trip("((a,b));")

round.trip("(a[&&NHX:foo=bar:this=that],(b,c));")