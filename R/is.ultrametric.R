## is.ultrametric.R (2009-05-10)

##   Test if a Tree is Ultrametric

## Copyright 2003-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.ultrametric <- function(phy, tol = .Machine$double.eps^0.5)
{
    if (!inherits(phy, "phylo"))
      stop('object "phy" is not of class "phylo".')
    if (is.null(phy$edge.length))
      stop("the tree has no branch lengths.")
    phy <- reorder(phy)
    n <- length(phy$tip.label)
    n.node <- phy$Nnode

    ## xx: vecteur donnant la distance d'un
    ## noeud ou tip � partir de la racine
    xx <- numeric(n + n.node)

    for (i in 1:dim(phy$edge)[1])
      xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]

    if (identical(all.equal.numeric(var(xx[1:n]),
                                    0, tolerance = tol), TRUE)) TRUE
    else FALSE
}
