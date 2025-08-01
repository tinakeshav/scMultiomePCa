###############################################################################
# MULTIOME_WNN_VlnDotPlots.R
# Script to generate violin, dot, bar, and UMAP plots for multiome data
#
# Usage: Rscript MULTIOME_WNN_VlnDotPlots.R --date ... [other args]
############

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dittoSeq)
})
source('scripts/sc_utilFunctions.R')

# ================================
# Argument Parsing
# ================================
parse_args <- function() {
  parser <- ArgumentParser(description = 'Generate violin, dot, bar plots, and UMAPs for multiome data')
  parser$add_argument('--date', type = 'character', required = TRUE,
                     help = 'Date script is running; used in saving object names')
  parser$add_argument('--figdir', type = 'character', required = TRUE,
                     help = 'Figures directory; must already exist')
  parser$add_argument('--datadir', type = 'character', required = TRUE,
                     help = 'Data directory; must already exist')
  parser$add_argument('--workingdir', type = 'character', required = TRUE,
                     help = 'Working directory; must already exist')
  parser$add_argument('--filteredMultiomeFile', type = 'character', nargs = '+', required = TRUE,
                     help = 'RDS file(s) to load (relative to data dir)')
  parser$add_argument('--sampleName', type = 'character', required = TRUE,
                     help = 'Shortened name to use for figures')
  parser$add_argument('--clustersofinterest', type = 'character', nargs = '+',
                     default = c('wsnn_res.0.4', 'wsnn_res.0.6', 'wsnn_res.0.8', 'wsnn_res.1'),
                     help = 'List of clusters of interest for figures')
  args <- parser$parse_args()
  return(args)
}

# ================================
# Gene Set Definitions
# ================================
get_gene_sets <- function() {
  geneset_oi <- c('AR', 'FOXA1', 'HOXB13', 'EPCAM', 'KRT18', 'TP63', 'LYZ', 'SELL', 'PTPRC', 'BTLA', 'CD2', 'CD3G', 'ACTA2', 'COL3A1')
  chengeneset <- stringr::str_split('ACTA2,PECAM1,VWF,ENG,CMA1,MS4A2,TPSAB1,TPSB2,AR,KRT19,KRT18,KRT8,TP63,KRT14,LYZ,FCGR3A,CSF1R,CD68,CD163,CD14,UCHL1,HAVCR2,PDCD1,CTLA4,CD8A,SELL,PTPRC,CD4,BTLA,IL2RA,IL7R,CCR7,CD28,CD27,SLAMF1,DPP4,CD7,CD2,CD3G,CD3E,CD3D', ',')[[1]]
  ln_set <- stringr::str_split('ACKR1,PLVAP,AQP1,CD74,VWF,RGCC,PECAM1,NTS,TFF3,PROX1,CD9,ACKR4,NUDT4,DCN,CXCL14,CCL19,COL1A2,TAGLN,ACTA2,FBLN1,COL3A1', ',')[[1]]
  karthaus_science <- stringr::str_split('CD19,MS4A1,CD4,CD8A,CD14,AIF1,CD16,CD57,CD31,VWF,PROX1,SOX10,NTM,MYH11,IGF1,ACTA2,SRL,COL5A2,FGF10,RSPO3,CXCL5,EPCAM,TP63,KRT13,KRT4,KLK3,CD26,DPP4,PLA2G2A,SCGB1A1,PSCA', ',')[[1]]
  mouse_ln <- stringr::str_split('PECAM1,PDPN,PROX1,FLT4,LYVE4,ITGA2B,MADCAM1,ACKR4', ',')[[1]]
  neuron <- stringr::str_split('NEUROD2,NEUROD4,GAD1,GAD2,GDA,SYP,PAX6,SYN1', ',')[[1]]
  all_genes <- unique(c(geneset_oi, chengeneset, ln_set, karthaus_science, mouse_ln, neuron))
  list(
    short = geneset_oi,
    cheng_prostate = chengeneset,
    LN = ln_set,
    karthaus_science = karthaus_science,
    mouse_LN = mouse_ln,
    neuron = neuron,
    all = all_genes
  )
}

# ================================
# Main Plotting Function
# ================================
run_all_plots <- function(sample, markers, clusters_of_interest, date, figures_dir, name_to_keep) {
  # Set up and normalize
  DefaultAssay(sample) <- 'SCT'
  DefaultAssay(sample) <- 'RNA'
  cat('[INFO] LogNormalizing RNA\n')
  sample <- NormalizeData(sample)

  for (cluster in clusters_of_interest) {
    for (marker in names(markers)) {
      gene_set <- markers[[marker]]
      # Violin Plots
      cat(sprintf('[INFO] Generating violin plots for cluster %s, marker set %s\n', cluster, marker))
      vlnplots_on_clusters(
        sample, gene_set, cluster = cluster, assay = 'RNA',
        middle_filename = paste0('VlnPlot_', marker, 'Genes'),
        date = date, figures_dir = figures_dir, name_to_keep = name_to_keep
      )
      # Dot Plots
      cat(sprintf('[INFO] Generating dot plots for cluster %s, marker set %s\n', cluster, marker))
      dotplots_on_clusters(
        sample, gene_set, cluster = cluster, assay = 'RNA',
        middle_filename = paste0('DotPlot_', marker, 'Genes'),
        date = date, figures_dir = figures_dir, name_to_keep = name_to_keep
      )
      # Bar Plots
      cat(sprintf('[INFO] Generating bar plots for cluster %s\n', cluster))
      barplot_on_clusters(
        sample, cluster, var = 'name', scale = 'percent',
        middle_filename = paste0(cluster, '_SampleContribution'),
        date = date, figures_dir = figures_dir, name_to_keep = name_to_keep
      )
      # WNN UMAP Feature Plots
      cat(sprintf('[INFO] Generating WNN UMAP feature plots for cluster %s, marker set %s\n', cluster, marker))
      wnn_integrated_plotname <- get_data_figure_name(
        figures_dir, name_to_keep,
        paste0('WNN_integrated-RNA_', marker, '_MARKERS_UMAP'),
        date, '.png'
      )
      wnn_umap_plot <- FeaturePlot(
        sample, features = gene_set, reduction = 'wnn.umap',
        max.cutoff = 5, ncol = 3, order = TRUE
      )
      ggsave(
        wnn_integrated_plotname, wnn_umap_plot,
        height = ceiling(3 * (length(gene_set) / 3)), width = 9, limitsize = FALSE
      )
    }
  }

  # UMAPs by name
  cat('[INFO] Generating UMAPs by sample name\n')
  plot_umaps_name <- get_data_figure_name(figures_dir, paste0('UMAPs-byName', name_to_keep), 'MULTIOME', date, '.png')
  p1 <- DimPlot(sample, reduction = 'umap.rna', group.by = 'name', label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle(paste0(name_to_keep, ' RNA'))
  p2 <- DimPlot(sample, reduction = 'umap.atac', group.by = 'name', label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle(paste0(name_to_keep, ' ATAC'))
  p3 <- DimPlot(sample, reduction = 'wnn.umap', group.by = 'name', label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle(paste0(name_to_keep, ' WNN'))
  umap_plot <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot_umaps_name, umap_plot, width = 15, height = 5)
}

# ================================
# run fxn
# ================================
main <- function() {
  args <- parse_args()
  cat('[INFO] Arguments parsed successfully\n')
  setwd(args$workingdir)
  cat(sprintf('[INFO] Working directory set to %s\n', args$workingdir))

  # Defensive: check file existence
  if (!file.exists(args$filteredMultiomeFile[1])) {
    cat(sprintf('[ERROR] File not found: %s\n', args$filteredMultiomeFile[1]))
    stop('Filtered multiome file not found.')
  }
  sample <- readRDS(args$filteredMultiomeFile[1])
  cat('[INFO] Multiome object loaded\n')

  markers <- get_gene_sets()
  run_all_plots(
    sample = sample,
    markers = markers,
    clusters_of_interest = args$clustersofinterest,
    date = args$date,
    figures_dir = args$figdir,
    name_to_keep = args$sampleName
  )
  cat('[INFO] All plots generated successfully\n')
}

if (sys.nframe() == 0) {
  main()
} 