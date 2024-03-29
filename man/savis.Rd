% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunSummary.R
\name{savis}
\alias{savis}
\title{savis}
\usage{
savis(
  expr_matrix,
  npcs = 20,
  nfeatures = 2000,
  hvg_method = NULL,
  distance_metric = "euclidean",
  cluster_method = "louvain",
  resolution = 0.1,
  resolution_sub = 0,
  adaptive = TRUE,
  max_stratification = 3,
  scale_factor_separation = 3,
  process_min_size = NULL,
  process_min_count = NULL,
  run_adaUMAP = TRUE,
  adjust_UMAP = TRUE,
  adjust_method = "all",
  adjust_rotate = TRUE,
  shrink_distance = TRUE,
  density_adjust = TRUE,
  density_adjust_via_global_umap = FALSE,
  adjust_scale_factor = 0.9,
  global_umap_embedding = NULL,
  check_differential = FALSE,
  verbose = TRUE,
  show_cluster = FALSE,
  return_cluster = FALSE,
  return_combined_PC = FALSE,
  verbose_more = FALSE,
  compressed_storage = TRUE,
  seed.use = 42L
)
}
\arguments{
\item{expr_matrix}{The expression COUNT matrix: gene(feature) as row; cell(sample) as column.}

\item{npcs}{The number of principle components will be computed as the initialization input of nonlinear low dimensional embeddings. Default is 20.}

\item{nfeatures}{The number of highly variable genes will be selected. Default is 2000.}

\item{hvg_method}{High Variable Gene Selection Method. Refer to manual of package 'mixhvg' and its function FindVariableFeaturesMix.}

\item{distance_metric}{The default is "euclidean". Recommend to use "euclidean" because we need to distinguish between global distance and local distance.}

\item{cluster_method}{The default is "louvain". User can choose from c("louvain","spectral"). But "louvain" performs much better.}

\item{resolution}{The resolution for the louvain clustering. The resolution ranges from 0 to 1. The lower resolution means conservative clustering(smaller number of clusters), while the higher resolution means aggressive clustering. The default is 0.5.}

\item{resolution_sub}{The resolution for the louvain clustering within subclusters, which means after the first step clustering and separation, we perform further clustering. The default is 0 because we are mainly interested in one step separation.}

\item{adaptive}{Whether we will run adaptive visualization. If adaptive is FALSE, we are just doing UMAP. The default is TRUE.}

\item{max_stratification}{The maximum level of stratification for subclustering. Set the maximum level to restrict too detailed exploration. The default is 3.}

\item{scale_factor_separation}{Scale factor used to distinguish global distance and local distance. The default is 3.}

\item{process_min_size}{The smallest size of cluster which we use to stop the process of subcluster evaluation. The clusters whose sizes are less than the cutoff will be not be further explored. The default is NULL.}

\item{process_min_count}{The processed clusters are determined by the ranking of the size of clusters. Only further explore the clusters whose rankings are less than the process_min_count. The default is NULL.}

\item{run_adaUMAP}{Whether we run the adaptive visualization. If the criterion is The default is TRUE.}

\item{adjust_UMAP}{The default is TRUE.}

\item{adjust_method}{The default is "all". Select from c("umap","mds").}

\item{adjust_rotate}{Adjust the rotation of each cluster. The default is TRUE.}

\item{shrink_distance}{Shrink distance of small clusters to avoid too much space in plot. The default is TRUE.}

\item{density_adjust}{Adjust density of plot. The default is TRUE.}

\item{density_adjust_via_global_umap}{Do density adjustment using the density from UMAP to make it comparable.}

\item{adjust_scale_factor}{Scale factor for adjustment. The default is 0.9. The smaller value (>0) means larger illustration of clusters.}

\item{global_umap_embedding}{The default is NULL.}

\item{check_differential}{The default is FALSE.}

\item{verbose}{The default is TRUE.}

\item{show_cluster}{The default is FALSE.}

\item{return_cluster}{The default is FALSE.}

\item{return_combined_PC}{The default is FALSE.}

\item{verbose_more}{More details are displayed. The default is FALSE.}

\item{compressed_storage}{Whether compress storage when returning. The default is FALSE.}

\item{seed.use}{The default is 42L}
}
\value{
nothing useful
}
\description{
savis: single-cell RNAseq adaptive visualiztaion
}
\details{
This function argument to the function
}
\examples{
a<-1

}
