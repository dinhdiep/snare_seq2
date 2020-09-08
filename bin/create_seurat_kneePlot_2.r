require('httr')
require('broom')
require('igraph')
require('png')
require('hdf5r')
require('rlang')
require('crayon')
require('digest')
require('assertthat')
require('glue')
require('purrr')
require('backports')
require('Seurat')
require('dplyr')
require('tidytext')
require('DropletUtils')
require('dropestr')
library(Seurat)
library(DropletUtils)

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

merge.dtn6 <- function(rds.file, umi_cutoff=100) {
  rds <- readRDS(rds.file)
  tidy.cm <- tidy(rds$cm)  

  dt.barcodes <- read.table("~/snare-seq2_scripts/config/R1_dTN6_pairs.txt",
                            header = F, stringsAsFactors = F)[,1]
  n6.barcodes <- read.table("~/snare-seq2_scripts/config/R1_dTN6_pairs.txt",
                            header = F, stringsAsFactors = F)[,2]
  
  n6.dt.map <- dt.barcodes
  names(n6.dt.map) <- n6.barcodes
  

  dtn6.barcodes <- sapply(tidy.cm$column, function(x) substr(x, 17, 24))
  r3r2.barcodes <- sapply(tidy.cm$column, function(x) substr(x, 1, 16))

  is.dt <- sapply(dtn6.barcodes, function(x) x %in% dt.barcodes) 
  merged.barcodes <- sapply(dtn6.barcodes, 
                        function(x) {
                          if (x %in% names(n6.dt.map)) { 
                            return(n6.dt.map[[x]])
                          } else { return(x) }
                        })
  
  new.cell.barcodes <- paste0(r3r2.barcodes, merged.barcodes)
  tidy.cm$column <- new.cell.barcodes
  combined.matrix <- tidy.cm %>% group_by(row, column) %>% summarise(mergedValue = sum(value, na.rm = TRUE)) #%>% cast_sparse(row, column, mergedValue)
  total_umis <- combined.matrix %>% group_by(column) %>% summarise(total_umi = sum(mergedValue, na.rm = TRUE)) %>% filter(total_umi > umi_cutoff)
  combined.matrix <- combined.matrix %>% filter(column %in% total_umis$column) %>% cast_sparse(row, column, mergedValue)

  combined.matrix
}

merge.dtn6.aligned.reads.umis.per.cell <- function(rds.file){
  rds <- readRDS(rds.file)
  info <- data.frame(aligned_reads_per_cell = rds$aligned_reads_per_cell, aligned_umis_per_cell = rds$aligned_umis_per_cell)
  info$row <- rownames(info)
  dt.barcodes <- read.table("~/snare-seq2_scripts/config/R1_dTN6_pairs.txt",
                            header = F, stringsAsFactors = F)[,1]
  n6.barcodes <- read.table("~/snare-seq2_scripts/config/R1_dTN6_pairs.txt",
                            header = F, stringsAsFactors = F)[,2]
  
  n6.dt.map <- dt.barcodes
  names(n6.dt.map) <- n6.barcodes
  
  
  dtn6.barcodes <- sapply(rownames(info), function(x) substr(x, 17, 24))
  r3r2.barcodes <- sapply(rownames(info), function(x) substr(x, 1, 16))
  
  is.dt <- sapply(dtn6.barcodes, function(x) x %in% dt.barcodes) 
  merged.barcodes <- sapply(dtn6.barcodes, 
                            function(x) {
                              if (x %in% names(n6.dt.map)) { 
                                return(n6.dt.map[[x]])
                              } else { return(x) }
                            })
  
  new.cell.barcodes <- paste0(r3r2.barcodes, merged.barcodes)
  info$row <- new.cell.barcodes
  combined.matrix <- info %>% group_by(row) %>% summarise(sum_aligned_reads_per_cell = sum(aligned_reads_per_cell, na.rm = TRUE), sum_aligned_umis_per_cell = sum(aligned_umis_per_cell, na.rm = TRUE))
  combined.matrix
  
}

process.raw.rds <- function(data.file, id){
  # merge DT and N6 primers
  # data.file = "./rdsFiles/CPTCA_20200227.SPLiT_N701_S1.rds"
  counts = merge.dtn6((data.file))
  reads_umis_info = merge.dtn6.aligned.reads.umis.per.cell(data.file)
  reads_umis_info = reads_umis_info[reads_umis_info$row %in% colnames(counts),]

  merged_aligned_reads_per_cell = reads_umis_info$sum_aligned_reads_per_cell
  names(merged_aligned_reads_per_cell) = reads_umis_info$row

  merged_umis_per_cell = reads_umis_info$sum_aligned_umis_per_cell
  names(merged_umis_per_cell) = reads_umis_info$row

  pdf(paste0(id,".cell_quality.plot.pdf"))
  p1 = PlotCellsNumberLogLog(merged_umis_per_cell, estimate.cells.number=T)
  p2 = PlotCellsNumberLine(merged_umis_per_cell, estimate.cells.number=T)
  p3 = PlotCellsNumberHist(merged_umis_per_cell, estimate.cells.number=T)
  multiplot(p1,p2,p3, cols = 2)
  dev.off()
  
  cells_df <- PrepareLqCellsData(counts, aligned.reads.per.cell = merged_aligned_reads_per_cell)
  #for (n in names(cells_df)) {
  #  smoothScatter(cells_df[[n]], xlab = "Cell rank", ylab = n, main = n)
  #}
  
  cells_number_manual <- list(min=450, max=nrow(cells_df))
  scores <- ScoreQualityData(merged_umis_per_cell, cells_df, cells_number_manual)
  
  seurat.obj <- CreateSeuratObject(counts = counts, project = id)
  seurat.obj <- subset(seurat.obj, cells = names(scores)[which(scores > 0.9)] )
  
  #printing colnames of SeuratObject
  head(colnames(seurat.obj))
  tail(colnames(seurat.obj))
  seurat.obj

  #saves seurat.raw as an raw.rds file
  saveRDS(seurat.obj, file = paste0(id, ".filtered.seurat.rds"))
}

args = commandArgs(trailingOnly=TRUE)
data.dir = args[1]

rdsFiles <- list.files(path = data.dir, pattern = ".rds$", full = TRUE)
rdsIds <- sapply(rdsFiles, function(x) gsub(".rds$", "", x))

sapply(1:length(rdsFiles), function(i) process.raw.rds(rdsFiles[i], rdsIds[i]))
#sapply(1:1, function(i) process.raw.rds(rdsFiles[i], rdsIds[i]))

