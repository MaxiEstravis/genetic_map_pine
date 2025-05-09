# modified from onemap version 2.1.1, the last to have the "dual" plot like in the Norway spruce paper, with the diagonal splitting LOD and rf

library(ggplot2)
library(reshape2)
library(ggnewscale)

rf_graph_table_dual <- function(input.seq, main = NULL, mrk.axis = "numbers", lab.xy = NULL) {
  # Check input
  if (!any(inherits(input.seq, "sequence"))) 
    stop(deparse(substitute(input.seq)), " is not an object of class 'sequence'")
  
  if (!(mrk.axis %in% c("names", "numbers", "none"))) 
    stop("This mrk.axis argument is not defined, choose 'names', 'numbers' or 'none'")
  
  # Get recombination fraction matrix
  if (inherits(input.seq$data.name, c("outcross", "f2"))) {
    mat <- t(onemap:::get_mat_rf_out(input.seq, LOD = TRUE, max.rf = 0.501, min.LOD = -0.1))
  } else {
    mat <- t(onemap:::get_mat_rf_in(input.seq, LOD = TRUE, max.rf = 0.501, min.LOD = -0.1))
  }
  
  mat[row(mat) > col(mat) & mat > 0.5] <- 0.5
  mat[row(mat) < col(mat)][mat[row(mat) < col(mat)] < 0.1] <- 0.1
  diag(mat) <- NA
  
  colnames(mat) <- rownames(mat) <- colnames(input.seq$data.name$geno)[input.seq$seq.num]
  
  if (mrk.axis == "numbers") 
    colnames(mat) <- rownames(mat) <- input.seq$seq.num
  
  n.mrk <- length(input.seq$seq.num)
  
  if (length(input.seq$seq.rf) > 1) {
    for (i in 1:(n.mrk - 1)) {
      mat[i + 1, i] <- input.seq$seq.rf[i]
    }
  }
  
  # Get LOD matrix
  mat.LOD <- mat
  mat.rf <- mat
  mat.LOD[lower.tri(mat.LOD)] <- t(mat.LOD)[lower.tri(mat.LOD)]
  mat.rf[upper.tri(mat.rf)] <- t(mat.rf)[upper.tri(mat.LOD)]
  
  # Melt matrices
  df.rf <- reshape2::melt(mat.rf, varnames = c("x", "y"), value.name = "rf")
  df.LOD <- reshape2::melt(mat.LOD, varnames = c("x", "y"), value.name = "LOD")
  
  # Merge
  df.graph <- merge(df.rf, df.LOD, by = c("x", "y"))
  
  # For axis ordering
  if (mrk.axis == "numbers") {
    df.graph$x <- factor(df.graph$x, levels = as.character(input.seq$seq.num))
    df.graph$y <- factor(df.graph$y, levels = as.character(input.seq$seq.num))
  }
  
  # Create dual plot with correct scales
  p <- ggplot(df.graph, aes(x, y)) +
    # Upper triangle: rf (dark blue = low rf, light blue = high rf)
    geom_tile(aes(fill = rf), data = subset(df.graph, as.numeric(x) < as.numeric(y))) +
    scale_fill_gradient(low = "darkblue", high = "lightblue", na.value = "white", name = "rf") +
    ggnewscale::new_scale_fill() +
    # Lower triangle: LOD (light pink = low LOD, dark red = high LOD)
    geom_tile(aes(fill = LOD), data = subset(df.graph, as.numeric(x) > as.numeric(y))) +
    scale_fill_gradient(low = "pink", high = "darkred", na.value = "white", name = "LOD") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Labels
  if (is.null(lab.xy)) {
    p <- p + labs(x = " ", y = " ")
  } else {
    if (length(lab.xy) != 2) {
      warning("You should give a character vector with two components to axis labels")
    } else {
      p <- p + labs(x = lab.xy[1], y = lab.xy[2])
    }
  }
  
  if (mrk.axis == "none") {
    p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  }
  
  if (!is.null(main)) {
    p <- p + ggtitle(main)
  }
  
  p
}


