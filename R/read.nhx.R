# Reads a tree with NHX annotations and places them as a list of tags in the phylo object.
# Use the method get.phylo.tag to retrieve tags from the tree.
read.nhx <- function(str) {
  nhx.matches <- gregexpr("(\\w+)?(:?\\d*\\.?\\d*)?\\[&&NHX.*?\\]", str)
  matches <- nhx.matches[[1]]
  match.pos <- as.numeric(matches)
  match.len <- attr(matches, 'match.length')

  if (match.len[1] == -1) {
    return(read.tree(text=str))
  }

  nhx.strings <- substring(str, match.pos, match.pos+match.len-1)

  labels <- gsub("\\[&&NHX.*", "", nhx.strings)
  labels.no.bl <- gsub(":.*", "", labels)

  # Go through and splice the stripped-down labels back into the string.
  for (i in 1:length(match.pos)) {
    new.label <- gsub(labels.no.bl[i], paste('zzz',i,'zzz',sep=''), labels[i])
    str <- paste(substr(str, 0, match.pos[i]-1 ), new.label, substr(str, match.pos[i] + match.len[i], nchar(str)), sep='')
    match.pos <- match.pos - match.len[i] + nchar(new.label)
  }

  # Parse the Phylo object from the cleaned-up string.
  tree <- read.tree(text=str)

  # Create a list of NHX annotations keyed by the node ID.
  tag.values <- gsub(".*\\[&&NHX:(.*)\\]", "\\1", nhx.strings)
  tagval.list <- strsplit(tag.values, ":")
  names(tagval.list) <- labels.no.bl

  map.list <- lapply(tagval.list, function(x) {
    list.out <- list()
    cur.list <- strsplit(x, "=")
    if (length(cur.list) > 0) {
      for (i in 1:length(cur.list)) {
        vec <- cur.list[[i]]
        list.out[vec[1]] = vec[2]
      }
    }
    return(list.out)
  })

  # Replace the labels with the true labels.
  tree$.tags <- list()
  for (i in 1:(tree$Nnode+length(tree$tip.label))) {
    tree$.tags[[i]] <- list()
  }

  for (i in 1:length(match.pos)) {
    cur.node <- node.with.label(tree, paste('zzz', i, 'zzz', sep=''))

    leaf <- node.is.leaf(tree, cur.node)
    real.node.name <- names(map.list)[i]
    if (leaf) {
      tree$tip.label[cur.node] <- real.node.name
    } else {
      tree$node.label[cur.node-length(tree$tip.label)] <- real.node.name
    }
    tree$.tags[[cur.node]] <- map.list[[i]]
  }

  #tree$.tags <- map.list
  return(tree)
}

# Tries to get a tag from the given node in the tree.
get.phylo.tag <- function(phylo, node, tag) {
  tag <- phylo$.tags[[node]][[tag]]
  if(is.null(tag)) {return('')}
  return(tag)
}

# Finds the node with a given label.
node.with.label <- function(tree,label) {
  all.labels <- c(tree$tip.label,tree$node.label)
  return(which(all.labels %in% label))
}

node.is.leaf <- function(phylo,node) {
  return(node <= length(phylo$tip.label))
}
