library(VennDiagram)
library(gridExtra)
library(grid)

pairings <- list(c("1", "1"), c("2", "2"), c("3", "4"), c("4", "5"), c("5", "6"), c("6", "7"),
                 c("7", "3"), c("8", "8"), c("9", "9"), c("10", "12"), c("11", "10"), c("12", "11"))

pdf("AC1017_loose_vs_strict.pdf", width = 11.69, height = 8.27)

grid.newpage()

pushViewport(viewport(layout = grid.layout(3, 4)))

plot_venn <- function(num1, num2, row, col) {
  file1 <- paste0("AC1017_loose_LG", num1, "_rec_map_redu.tsv")
  file2 <- paste0("AC1017_strict_LG", num2, "_rec_map_redu.tsv")
  
  markers1 <- read.table(file1, header = TRUE, sep = "\t")$id
  markers2 <- read.table(file2, header = TRUE, sep = "\t")$id
  
  venn.plot <- draw.pairwise.venn(
    area1 = length(markers1),
    area2 = length(markers2),
    cross.area = length(intersect(markers1, markers2)),
    category = c(paste0("Loose_LG", num1), paste0("Strict_LG", num2)),
    fill = c("blue", "red"),
    lty = "blank",
    cex = 1.5,
    cat.cex = 1.0,
    cat.pos = c(-20, 20),
    #cat.dist = c(0.06, 0.06),
    ind = FALSE
  )
  
  vp <- viewport(layout.pos.row = row, layout.pos.col = col)
  pushViewport(vp)
  grid.draw(venn.plot)
  popViewport()
}

for (i in seq_along(pairings)) {
  row <- ((i - 1) %/% 4) + 1  # Compute row (1 to 3)
  col <- ((i - 1) %% 4) + 1   # Compute column (1 to 4)
  plot_venn(pairings[[i]][1], pairings[[i]][2], row, col)
}

for (r in 1:3) {
  grid.lines(x = c(0, 1), y = c(r / 3, r / 3), gp = gpar(col = "black", lwd = 1))
}
for (c in 1:3) {
  grid.lines(x = c(c / 4, c / 4), y = c(0, 1), gp = gpar(col = "black", lwd = 1))
}

popViewport()

dev.off()


## continues in genetic_map_outliers.R