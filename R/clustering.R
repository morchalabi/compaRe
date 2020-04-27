# _________________________________Roxygen_________________________________________________________
#' Clustering samples with controls
#'
#' It needs some samples as controls and returns maximal cliques as clusters and highly interconnected regions as
#' communities.
#'
#' @param simMat_ Similarity matrix as a matrix object
#' @param controls_ Indices of control samples in simMat_ as a string
#' @param thresh_ Minimum for two samples being considered similar as a numeric
#' @param smpl_graph If sample graph must be output (default True)
#'
#' @details TBD
#'
#' @return samples_table: samples with their assigned community (dense region)
#' @return cliques: cluster of samples
#' @return samples_graph: graph of samples colored by community
#'
#' @seealso \code{\link{compare}} for measuring spatial similarity between two samples.
#'
#' @author Morteza H. Chalabi, \url{mor.chalabi@@gmail.com}
#'
#' @examples
#' require(compaRe)
#' require(igraph)
#'
#' rm(list = ls())
#'
#' # Step 1: Reading in similarity matrix
#'
#' data(package = 'compaRe', list = c('simMat'))
#'
#' # Step 2: Clustering
#'
#' out_ = compaRe::clustering(simMat_ = as.matrix(simMat),
#'                            controls_ = "7,23,30,35,55,106,164,193,214,228,246,254,258,286,343,351,414,444,467,489,540",
#'                            thresh_ = NULL, smpl_graph = T)
#'
#' # Step 3: Plotting dispersion graph in current directory as a pdf file
#'
#' pdf(file = 'sample_graph.pdf', width = 100, height = 100)
#' par(mar = c(0,0,0,0))
#' plot(out_$samples_graph)
#' graphics.off()
#'
#' @export

clustering = function(simMat_ = NULL, controls_ = NULL, thresh_ = NULL, smpl_graph = TRUE)
{
  require(igraph)     # if igraph pacakge is already installed

  output_ = list()    # output list

  # STEP 1: Checkpoint for controling input arguments ####

  # reading in similarity matrix

  if(is.null(simMat_))
  {
    message('\nError: simMat_ cannot be left empty!')
    quit(save = 'no')
  }
  if( is.null(rownames(simMat_)) | is.null(colnames(simMat_)) )     # if similarity matrix has no row/column names
  {
    rownames(simMat_) = colnames(simMat_) = 1:nrow(simMat_)
  }

  # checking controls

  controls_ = as.integer(strsplit(controls_,'[,]')[[1]])
  if(is.null(controls_))
  {
    message('\nError: controls_ cannot be left empty!')
    quit(save = 'no')
  }

  # setting similarity cutoff

  if(is.null(thresh_))            # if thresh_ is not set by user, then it must be inferred from control samples
  {
    if(length(controls_) < 2)     # there must be at least 2 controls to infer similarity cutoff
    {
      message('\nError: for inferring similarity cutoff, there must be at least 2 control samples!')
      quit(save = 'no')
    }
    cntSimMat = simMat_[controls_, controls_]      # similarity submatirx of controls
    diag(cntSimMat) = 0
    cntSimVlas = apply(X = cntSimMat, MARGIN = 1, FUN = max)
    thresh_ = min(cntSimVlas)
  }
  message('\nSimilarity threshold set to: ', thresh_)

  # chekig if smpl_graph is requested

  simMat_org = NULL
  if(smpl_graph) { simMat_org = simMat_}

  # STEP 2: Identifying samples similar enough to controls ####

  smpls_ = rownames(simMat_)                           # sample IDs
  dt_ = data.frame(sample = smpls_,                    # dt_ is the output table
                   community = 0,                      # components/community (connected subgraphs)
                   sim_vs_control = 0,                 # median similarity of each sample with controls
                   sim_vs_all = rowMeans(simMat_),     # mean similarity of each sample with all samples including controls
                   row.names = smpls_,
                   stringsAsFactors = F)
  mat_tmp = simMat_
  diag(mat_tmp) = 0
  simsVsCntrl = apply(X = mat_tmp[controls_,], MARGIN = 2, FUN = max)     # similarity of a sample with controls
  dt_[names(simsVsCntrl), "sim_vs_control"] = simsVsCntrl                 # updating output table with sim_vs_control values

  # updating control samples

  nonCtrls = which(dt_$sim_vs_control < thresh_)
  if(length(nonCtrls) != 0)
  {
    simsVsCntrl = apply(X = simMat_[-nonCtrls, nonCtrls], MARGIN = 2, FUN = max)     # updating similarities vs new controls
    simMat_ = simMat_[nonCtrls, nonCtrls]
  }
  simMat_ = cbind(simMat_, Control = 100)
  simMat_ = rbind(simMat_, Control = 100)
  simMat_["Control", names(simsVsCntrl)] = simMat_[names(simsVsCntrl),"Control"] = simsVsCntrl

  # STEP 3: Identifying clusters (maximal cliques) ####
  message('Identifying clusters')

  # step 3.1: finding communities; controls should be removed from graph before finding communities because
  # if 2 non-control nodes form a triangle with control node, the entire triangle is reported as a cluster while the 2 non-controls are desired

  adj_mat = simMat_[-nrow(simMat_), -ncol(simMat_)]      # last row and column of dt_ is control node
  adj_mat[ adj_mat < thresh_ ] = 0
  g_ = graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected', diag = F, weighted = T)
  comms_ = components(graph = g_)     # extracting conencted subgraphs (aka components/communities)
  dt_[names(comms_$membership), "community"] = comms_$membership

  output_[['samples_table']] = dt_

  # step 3.2: xtracting mximal cliques from each identified community

  cliq_tbl = list()
  j_ = 1      # cliques counter
  for(comm in 1:comms_$no)
  {
    # step 3.2.1: extracting current community

    smpls_ = names(comms_$membership[comms_$membership %in% comm])                                              # samples in this community
    adj_mat = simMat_[smpls_, smpls_, drop = F]                                                                  # adjacency matrix of current community considering sim threshold
    adj_mat[adj_mat < thresh_] = 0
    g_comm = graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected', diag = F, weighted = T)      # current community subgraph
    cliques_ = max_cliques(graph = g_comm)                                                                      # maximal cliques (unextendible complete subgraphs)

    # step 3.2.2: replacing sample indices of cliques with sample IDs & writing to file

    for(clq_ in cliques_)
    {
      cliq_tbl[[j_]] = data.frame(Cliques = paste0(clq_$name,collapse = ','), Community = comm)
      j_ = j_ + 1
    }
  }
  output_[['cliques']] = do.call(rbind, cliq_tbl)

  # STEP 4: plotting sample graph ####

  if(smpl_graph)
  {
    message('Building samples graph')

    adj_mat = simMat_org
    adj_mat[ adj_mat < thresh_ ] = 0
    g_ = graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected', diag = F, weighted = T)

    # graph atts
    g_$layout = layout_nicely(graph = g_, dim = 3)

    # Vertex atts
    # extracting communities with at least 2 nodes and coloring them, othrewise leaving them grey

    bigComs = which(2 <= comms_$csize)      # communities with size of at least 2 nodes
    cols_ = c('grey', colorRampPalette(colors = c('red','green','blue','purple'))(length(bigComs)) )# is for community
    names(cols_) = c('0',bigComs)
    for(i_ in 1:length(V(g_)$name))
    {
      comm_ = dt_[ V(g_)$name[i_], "community" ]
      if(comm_ %in% bigComs)
      {
        V(g_)[i_]$color = adjustcolor(col = cols_[as.character(comm_)], alpha.f = .6)
      }else
      {
        V(g_)[i_]$color = adjustcolor(col = 'grey', alpha.f = .2)
      }
    }
    V(g_)$size <- 2
    V(g_)$frame.color = NA
    V(g_)$label = NA

    # edge atts
    E(g_)$width = 0.3
    E(g_)$color = adjustcolor(col = 'grey89', alpha.f = .2)

    output_[['samples_graph']] = g_
  }
  return(output_)
}
