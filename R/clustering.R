#________________________________________Roxygen________________________________________________
#
#' Clustering samples with controls
#'
#' It needs some samples as controls and returns maximal cliques as clusters and highly interconnected regions as
#' communities.
#'
#' @param simMat_ Similarity matrix as a matrix object
#' @param controls_ Indices of control samples in simMat_ as a string
#' @param thresh_ Minimum for two samples being considered similar as a numeric
#' @param smpl_graph If sample graph must be output (default True)
#' @param sim_graph If similarity graph must be output (default True)
#'
#' @details TBD
#'
#' @return samples_table: samples with their assigned community (dense region)
#' @return cliques: cluster of samples
#' @return samples_graph: graph of samples colored by community
#' @return sim_graph: graph possibly showing sample lineages
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
#' # Step 1: Reading in similarity matrix ####
#'
#' data(package = 'compaRe', list = c('simMat'))
#'
#' # Step 2: Clustering ####
#'
#' out_ = compaRe::clustering(simMat_ = as.matrix(simMat),
#'                            controls_ = "7,23,30,35,55,106,164,193,214,228,246,254,258,286,343,351,414,444,467,489,540",
#'                            thresh_ = NULL, smpl_graph = T, sim_graph = T)
#'
#' # Step 3: Plotting samples graph in current directory as a pdf file ####
#'
#' # graph atts
#' g_ = out_$samples_graph
#' g_$layout = layout_nicely(graph = g_, dim = 3)
#'
#' # vertex atts
#' x_rng = range(g_$layout[,1])
#' y_rng = range(g_$layout[,2])
#' comms_ = unique(V(g_)$comm)
#' for(comm_ in comms_)
#' {
#'   if(comm_ ==  0) { next() }
#'
#'   inds_ = which(V(g_)$comm %in% comm_)
#'   if(1 < length(inds_))
#'   {
#'     centroid_ = c(sample(seq(x_rng[1],x_rng[2],by = 0.1),1), sample(seq(y_rng[1],y_rng[2],by = 0.1),1))
#'     anch_angs = seq(0, 2*pi, length.out = length(inds_)+1)
#'     r_ = 1
#'     anch_x = r_*cos(anch_angs[-1])+centroid_[1]
#'     anch_y = r_*sin(anch_angs[-1])+centroid_[2]
#'     g_$layout[inds_,1] = anch_x
#'     g_$layout[inds_,2] = anch_y
#'   }
#' }
#' cols_ = colorRampPalette(colors = c('red','green','blue','purple','orange','pink','yellow'))(length(comms_))
#' V(g_)$color = adjustcolor(col = cols_[V(g_)$comm+1], alpha.f = .7)
#' V(g_)$color[which(V(g_)$comm %in% 0)] = adjustcolor(col = 'grey', alpha.f = .4)
#' V(g_)$size <- 4
#' V(g_)$frame.color = NA
#' V(g_)$label = out_$samples_table[V(g_)$name, 'community']
#' V(g_)$label[which(V(g_)$label %in% 0)] = NA
#' V(g_)$label.cex = 4
#' V(g_)$label.font = 2
#'
#' # edge atts
#' E(g_)$width = 0.3
#' E(g_)$color = adjustcolor(col = 'grey', alpha.f = .3)
#' E(g_)$color[E(g_)$intra_comm] = 'black'
#'
#' # plotting
#' pdf(file = 'sample_graph.pdf', width = 70, height = 70)
#' par(mai = c(0, 0, 0,0))
#' plot(g_)
#' graphics.off()
#'
#' # Step 4: Plotting similarity graph ####
#'
#' # graph atts
#' g_ = out_$similarity_graph
#' g_$layout = layout_nicely(graph = g_, dim = 2)
#'
#' # vertex atts
#' V(g_)$color = 'grey'
#' V(g_)$size <- 0.15
#' V(g_)$frame.color = NA
#' V(g_)$label.cex = 1
#' V(g_)$label.font = 2
#' V(g_)$label.dist = 0
#' V(g_)$label.degree = pi/2
#' V(g_)$label.color = adjustcolor(col = 'black', alpha.f = .6)
#'
#' # edge atts
#' cols_ = colorRampPalette(colors = c('red', 'blue'))(length(E(g_)))
#' names(cols_) = sort(E(g_)$weight)     # lower values are assigned to red shades
#' E(g_)$color = cols_[as.character(E(g_)$weight)]
#' E(g_)$width = 1
#' E(g_)$arrow.size = 0.3
#' E(g_)$label = round(E(g_)$weight,1)
#' E(g_)$label.cex = 0.7
#' E(g_)$label.font = 2
#' E(g_)$label.color = 'darkgreen'
#'
#' # plotting
#' pdf(file = 'similarity_graph.pdf', width = 70, height = 70)
#' par(mai = c(0, 0, 0,0))
#' plot(g_, add = F, mark.groups = which(V(g_)$name %in% 'Control'), mark.col = 'lightgreen', mark.expand = 2, mark.border = NA, directed = F)
#' graphics.off()
#'
#' @export

clustering = function(simMat_ = NULL, controls_ = NULL, thresh_ = NULL, smpl_graph = TRUE, sim_graph = TRUE)
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

    # finding threshold using maximum spanning tree
    # it is equivalent to a for loop starting with min similarity in ascending order and
    # stop when graph is not connected anymore using igraphh::is.connected()

    g_ = graph_from_adjacency_matrix(adjmatrix = -simMat_[controls_, controls_], mode = 'undirected',weighted = T, diag = F)      # graph with negative weights
    g_ = mst(graph = g_)                                                                                                          # maximum spanning tree
    thresh_ = min(-E(g_)$weight)     # threshold is the maximum control similarity for which control graph remains a tree
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

  # STEP 4: Sample graph ####

  if(smpl_graph)
  {
    message('Constructing samples graph')

    simMat_org[simMat_org < thresh_] = 0      # removing insignificant edges
    g_ = graph_from_adjacency_matrix(adjmatrix = simMat_org, mode = 'undirected', diag = F, weighted = T)

    # assigning community number of each node

    V(g_)$comm = dt_[ V(g_)$name, "community" ]     # adding community attribute to vertices

    # marking edges control nodes

    edges_ = as_edgelist(g_)      # all edges in similarity graph after removing insignificant ones
    E(g_)$intra_comm = FALSE      # adding intra-community edge attribute to each edge
    inds_ = which(dt_[edges_[,1],"community"] == dt_[edges_[,2],"community"] &      # edges that connect nodes of the same community
                  !dt_[edges_[,1],"community"] %in% 0)                              # except for control nodes
    E(g_)$intra_comm[inds_] = TRUE

    output_[['samples_graph']] = g_
  }

  # STEP 5: Similarity graph ####

  if(sim_graph)
  {
    message('Constructing similarity graph')

    # finding nearest (most similar and correlated) node to each node taking control node as root

    diag(simMat_) = 0                                 # to avoid self-loops
    adj_ = simMat_                                    # similarity graph adjacency matrix
    adj_['Control',] = 0                              # control node shoul not point to other nodes
    for(row_ in 1:(nrow(simMat_)-1))                  # last row/column is control
    {
      nns_ = order(simMat_[row_,], decreasing = T)[1:min(3,ncol(simMat_))]      # 3 nearest neighbors
      max_ = simMat_[row_, nns_[1]]                                             # neighbor node with the highest similarity with node_A
      nns_ = nns_[which( abs(simMat_[row_, nns_] - max_) <= 2)]                 # considers only nns_ within [max-2, max] where 0 < max-2
      cor_ = double(length = length(nns_))                                      # a vector to store correlations of nns_ with current node
      j_ = 1                                                                    # counter of cor_
      for(nn_ in nns_)                                                          # for each nns_
      {
        node_A = simMat_[row_,-c(row_,nn_)]                                     # sim values of node_A with all other nodes except for node_B
        node_B = simMat_[nn_, -c(row_,nn_)]
        cor_[j_] = cor(node_A, node_B)
        j_ = j_ + 1
      }
      nn_ = nns_[order(cor_, decreasing = T)[1]]                                # which of nns_ had highest correlation with node_A?
      adj_[row_, -nn_] = 0                                                      # removing all edges from node_A except for nn_ -> node_A
    }
    g_ = graph_from_adjacency_matrix(adjmatrix = adj_, mode = 'directed', weighted = T)

    output_[['similarity_graph']] = g_
  }

  return(output_)
}
