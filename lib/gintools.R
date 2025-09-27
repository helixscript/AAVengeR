# The functions below were originally found in the gintools package
# and have been updated to work with more recent Biconductor packages.
# https://github.com/cnobles/gintools
#-------------------------------------------------------------------------------

refine_breakpoints <- function(input.sites, counts.col = NULL, min.gap = 1L,
                               sata.gap = 3L, details = FALSE, quiet = TRUE){
  
  stopifnot(class(input.sites) == "GRanges")
  sites <- GenomicRanges::granges(input.sites)
  
  if(!quiet){
    message(paste0("Refining ", length(sites), " break points."))
    message(
      "Generating initial graph by connecting positions with 1 nt difference.")
  }
  
  # Retain original order
  sites$ori.order <- seq_along(sites)
  
  # Identify counts or abundance info, assume 1 if not given, but check for
  # a column named counts first and use, otherwise error out if column is not
  # found.
  
  if(!is.null(counts.col)){
    if(counts.col %in% names(GenomicRanges::mcols(input.sites))){
      counts_pos <- grep(counts.col, names(GenomicRanges::mcols(input.sites)))
      sites$counts <- GenomicRanges::mcols(input.sites)[,counts_pos]
    }else{
      stop("Could not identify 'counts' column.")
    }
  }else{
    if(!quiet) message("Assuming abundance of 1 for each row of sites object.")
    sites$counts <- rep(1, length(sites))
  }
  
  # Reduce the genomic locations of break points down to only unique positions,
  # and identify the abundance of the positions
  
  red_sites <- GenomicRanges::reduce(
    GenomicRanges::flank(sites, -1, start = FALSE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red_sites$breakpoint.id <- seq_along(red_sites)
  
  # Summarise count data for reduced site object
  red_counts <- GenomicRanges::mcols(sites) %>%
    as.data.frame(row.names = NULL)
  red_counts <- red_counts[unlist(red_sites$revmap),] %>%
    dplyr::mutate(grp = as.numeric(S4Vectors::Rle(
      values = seq_along(red_sites),
      lengths = S4Vectors::elementNROWS(red_sites$revmap)))) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(abund = sum(counts)) %>%
    dplyr::ungroup()
  
  red_sites$abund <- red_counts$abund
  
  # Identify which sites are adjacent to each other
  red_hits <- GenomicRanges::as.data.frame(
    GenomicRanges::findOverlaps(
      red_sites, maxgap = min.gap - 1L, drop.self = TRUE
    )
  )
  
  # Organize data and filter for constructing directed graph
  red_hits <- red_hits %>%
    dplyr::mutate(
      hits.q = queryHits,
      hits.s = subjectHits,
      pos.q = GenomicRanges::start(red_sites[hits.q]),
      pos.s = GenomicRanges::start(red_sites[hits.s]),
      abund.q = red_sites[hits.q]$abund,
      abund.s = red_sites[hits.s]$abund,
      strand = as.vector(GenomicRanges::strand(red_sites[hits.q])),
      is.downstream = ifelse(strand == "+", pos.q > pos.s, pos.q < pos.s),
      keep = abund.q > abund.s,
      keep = ifelse(abund.q == abund.s, is.downstream, keep)) %>%
    dplyr::filter(keep)
  
  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  
  g <- igraph::make_empty_graph(n = length(red_sites), directed = TRUE) %>%
    igraph::add_edges(with(red_hits, vzip(hits.q, hits.s)))
  
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet) message(paste0("Initial cluster count: ", igraph::clusters(g)$no))
  
  # Identify satalite positions that should be included in clusters up to the
  # sata.gap max. This portion of the function tries to reach out to up to the
  # sata.gap distance from the boundry of a cluster to see if there are any
  # further "satalite" positions that have been annotated. It does this by
  # iterively increasing the size from 2nt to the sata.gap by 1nt increments.
  
  if(!quiet){
    message(paste0(
      "Connecting satalite positions up to ", sata.gap, " nt apart."))
  }
  
  null <- lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red_sites, g, gap, "downstream")
    red_sites$clus.id <<- igraph::clusters(g)$membership
  })
  
  if(!quiet){
    message("Clusters after satalite connecting: ", igraph::clusters(g)$no)
  }
  
  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.
  
  g <- break_connecting_source_paths(red_sites, g, "downstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet){
    message(paste0("Final break point cluster count: ", igraph::clusters(g)$no))
  }
  
  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and correct the original
  # sites object with the adjusted breakpoints.
  
  src_nodes <- sources(g)
  
  clus_data <- data.frame(
    "clus.id" = seq_along(igraph::clusters(g)$csize),
    "chr" = GenomicRanges::seqnames(red_sites[src_nodes]),
    "strand" = GenomicRanges::strand(red_sites[src_nodes]),
    "breakpoint" = GenomicRanges::start(red_sites[src_nodes]),
    "width" = GenomicRanges::width(unlist(range(
      GenomicRanges::split(red_sites, red_sites$clus.id))))
  )
  
  sites <- sites[unlist(IRanges::as.list(red_sites$revmap))]
  sites$clus.id <- as.numeric(S4Vectors::Rle(
    igraph::clusters(g)$membership,
    S4Vectors::elementNROWS(red_sites$revmap)))
  sites$called.bp <- ifelse(
    GenomicRanges::strand(sites) == "+",
    GenomicRanges::end(sites),
    GenomicRanges::start(sites))
  sites$adj.bp <- clus_data[
    match(sites$clus.id, clus_data$clus.id), "breakpoint"]
  
  if(!quiet){
    message(paste0(
      "Cumulative displacement per range: ",
      round(
        sum( abs(sites$called.bp - sites$adj.bp) ) / length(sites),
        digits = 1)))
  }
  
  GenomicRanges::ranges(sites) <- IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(sites) == "+",
      GenomicRanges::start(sites),
      sites$adj.bp),
    end = ifelse(
      GenomicRanges::strand(sites) == "+",
      sites$adj.bp,
      GenomicRanges::end(sites))
  )
  
  # Transfer range information from sites to output_sites object
  sites <- sites[order(sites$ori.order)]
  output_sites <- input.sites
  GenomicRanges::ranges(output_sites) <- GenomicRanges::ranges(sites)
  
  if(details){
    GenomicRanges::mcols(output_sites)$called.pos <- GenomicRanges::mcols(
      sites)$called.bp
    GenomicRanges::mcols(output_sites)$adj.pos <- GenomicRanges::mcols(
      sites)$adj.bp
  }
  
  output_sites
}


standardize_sites <- function(input.sites, counts.col = NULL, min.gap = 1L,
                              sata.gap = 5L, details = FALSE, quiet = TRUE){
  stopifnot(class(input.sites) == "GRanges")
  sites <- GenomicRanges::granges(input.sites)
  
  # Retain original order
  sites$ori.order <- seq_along(sites)
  
  # Identify counts or abundance info, assume 1 if not given, but check for
  # a column named counts first and use, otherwise error out if column is not
  # found.
  
  if(!is.null(counts.col)){
    if(counts.col %in% names(GenomicRanges::mcols(input.sites))){
      counts_pos <- grep(counts.col, names(GenomicRanges::mcols(input.sites)))
      sites$counts <- GenomicRanges::mcols(input.sites)[,counts_pos]
    }else{
      stop("Could not identify 'counts' column.")
    }
  }else{
    if(!quiet) message("Assuming abundance of 1 for each row of sites object.")
    sites$counts <- rep(1, length(input.sites))
  }
  
  # Start by reducing the sites object down to only the site's positions
  # and store the revmap for going back to the original sites object.
  if(!quiet){
    message(paste0("Standardizing ", length(sites), " positions."))
    message(
      "Generating initial graph by connecting positions with 1 nt difference.")
  }
  
  # Construct reduced sites object to initialize the processing
  red_sites <- GenomicRanges::reduce(
    GenomicRanges::flank(sites, -1, start = TRUE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red_sites$site.id <- seq_along(red_sites)
  
  # Summarise count data for reduced site object
  red_counts <- GenomicRanges::mcols(sites) %>%
    as.data.frame(row.names = NULL)
  red_counts <- red_counts[unlist(red_sites$revmap),] %>%
    dplyr::mutate(grp = as.numeric(S4Vectors::Rle(
      values = seq_along(red_sites),
      lengths = S4Vectors::elementNROWS(red_sites$revmap)))) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(abund = sum(counts)) %>%
    dplyr::ungroup()
  
  red_sites$abund <- red_counts$abund
  
  # Identify which sites are adjacent to each other
  red_hits <- GenomicRanges::as.data.frame(
    GenomicRanges::findOverlaps(
      red_sites, maxgap = min.gap - 1L, drop.self = TRUE
    )
  )
  
  # Organize data and filter for constructing directed graph
  red_hits <- red_hits %>%
    dplyr::mutate(
      hits.q = queryHits,
      hits.s = subjectHits,
      pos.q = GenomicRanges::start(red_sites[hits.q]),
      pos.s = GenomicRanges::start(red_sites[hits.s]),
      abund.q = red_sites[hits.q]$abund,
      abund.s = red_sites[hits.s]$abund,
      strand = as.vector(GenomicRanges::strand(red_sites[hits.q])),
      is.upstream = ifelse(strand == "+", pos.q < pos.s, pos.q > pos.s),
      keep = abund.q > abund.s,
      keep = ifelse(abund.q == abund.s, is.upstream, keep)) %>%
    dplyr::filter(keep)
  
  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  g <- igraph::make_empty_graph(n = length(red_sites), directed = TRUE) %>%
    igraph::add_edges(with(red_hits, vzip(hits.q, hits.s)))
  
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet) message(paste0("Initial cluster count: ", igraph::clusters(g)$no))
  
  # When clusters are formed of positions with equivalent abundance, both
  # directional edges are created. This type of cluster does not have sources
  # or sinks. The following set resolves which directed edges to remove, but
  # imposes an upstream bias.
  
  #  g <- resolve_source_deficiency(red_sites, g, bias = "upstream")
  #  red_sites$clus.id <- clusters(g)$membership
  
  # Identify satalite positions that should be included in clusters up to 5nt
  # away. This portion of the function tries to reach out to up to 5nt from the
  # boundry of a cluster to see if there are any further "satalite" positions
  # that have been annotated. It does this by iterively increasing the size
  # from 2nt to 5nt by 1nt increments.
  
  if(!quiet){
    message(paste0(
      "Connecting satalite positions up to ", sata.gap, " nt apart."))
  }
  
  null <- lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red_sites, g, gap, "upstream")
    red_sites$clus.id <<- igraph::clusters(g)$membership
  })
  
  if(!quiet){
    message("Clusters after satalite connecting: ", igraph::clusters(g)$no)
  }
  
  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.
  
  g <- break_connecting_source_paths(red_sites, g, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet){
    message(paste0("Clusters after clipping: ", igraph::clusters(g)$no))
  }
  
  # In the end, sources that are within the range of satalites (5L default),
  # should be grouped together. These sources should be connected by an edge,
  # pointing toward the larger source. The chosen source will have the highest
  # abundance of widths / sonic breaks / fragment lengths.
  
  if(!quiet){
    message("Connecting clusters with source nodes within ", sata.gap, " nt.")
  }
  
  g <- connect_adjacent_clusters(red_sites, g, gap = sata.gap, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  if(!quiet) message(paste0("Final cluster count: ", igraph::clusters(g)$no))
  
  # If wide clusters are generated, these can confound the cluster source. For
  # this reason, the "true" source for the cluster will be resolved by
  # picking the source with the greatest number of fragment lengths, where ties
  # are decided by randon picking.
  
  g <- resolve_cluster_sources(red_sites, g, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership
  
  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and adjust the original
  # sites object with the standardized positions.
  
  src.nodes <- sources(g)
  
  clus.data <- data.frame(
    "clus.id" = seq_along(igraph::clusters(g)$csize),
    "chr" = GenomicRanges::seqnames(red_sites[src.nodes]),
    "strand" = GenomicRanges::strand(red_sites[src.nodes]),
    "position" = GenomicRanges::start(red_sites[src.nodes]),
    "width" = GenomicRanges::width(unlist(range(
      GenomicRanges::split(red_sites, red_sites$clus.id))))
  )
  
  # Index the sites object and add adjusted data to metadata
  sites <- sites[unlist(IRanges::as.list(red_sites$revmap))]
  sites$clus.id <- as.numeric(S4Vectors::Rle(
    igraph::clusters(g)$membership,
    S4Vectors::elementNROWS(red_sites$revmap)))
  sites$called.pos <- ifelse(
    GenomicRanges::strand(sites) == "+",
    GenomicRanges::start(sites),
    GenomicRanges::end(sites))
  sites$adj.pos <- clus.data[
    match(sites$clus.id, clus.data$clus.id), "position"]
  
  if(!quiet){
    message(paste0(
      "Cumulative displacement per range: ",
      sum(abs(sites$called.pos - sites$adj.pos))/length(sites)))
  }
  
  # Adjust range information within sites object
  GenomicRanges::ranges(sites) <- IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(sites) == "+",
      sites$adj.pos,
      GenomicRanges::start(sites)),
    end = ifelse(
      GenomicRanges::strand(sites) == "+",
      GenomicRanges::end(sites),
      sites$adj.pos)
  )
  
  # Transfer range information from sites to output_sites object
  sites <- sites[order(sites$ori.order)]
  output_sites <- input.sites
  GenomicRanges::ranges(output_sites) <- GenomicRanges::ranges(sites)
  
  if(details){
    GenomicRanges::mcols(output_sites)$called.pos <- GenomicRanges::mcols(
      sites)$called.pos
    GenomicRanges::mcols(output_sites)$adj.pos <- GenomicRanges::mcols(
      sites)$adj.pos
  }
  
  output_sites
}




vzip <- function(...){
  v <- list(...)
  if(length(unique(sapply(v, length))) > 1){
    stop("All input vectors are not of equal length.")
  }
  as.vector(t(matrix(unlist(v), ncol = length(v))))
}


connect_satalite_vertices <- function(red.sites, graph, gap, bias){
  clus_mem <- igraph::clusters(graph)$membership
  
  clus_ranges <- unlist(GenomicRanges::reduce(
    GenomicRanges::split(red.sites, clus_mem),
    min.gapwidth = (gap-1)
  ))
  
  sata_hits <- as.data.frame(
    GenomicRanges::findOverlaps(
      clus_ranges, maxgap = gap - 1L, drop.self = TRUE
    )
  )
  
  names(sata_hits) <- c("source.clus", "sata.clus")
  
  red_df <- GenomicRanges::as.data.frame(red.sites)
  
  if(nrow(sata_hits) > 0){
    clus_data <- red_df %>%
      dplyr::group_by(clus.id) %>%
      dplyr::summarize(
        clus.pos.mean = as.integer(mean(start)),
        min.abund = min(abund),
        sum.abund = sum(abund))
    
    if(bias == "upstream"){
      sata_hits <- sata_hits %>%
        dplyr::mutate(
          source.pos = clus_data[source.clus,]$clus.pos.mean,
          sata.pos = clus_data[sata.clus,]$clus.pos.mean,
          min.src.abund = clus_data[.$source.clus,]$min.abund,
          min.sat.abund = clus_data[.$sata.clus,]$min.abund,
          sum.src.abund = clus_data[.$source.clus,]$sum.abund,
          sum.sat.abund = clus_data[.$sata.clus,]$sum.abund,
          is.upstream = source.pos < sata.pos) %>%
        dplyr::filter(
          as.integer(min.src.abund) >= as.integer(min.sat.abund),
          sum.src.abund > sum.sat.abund,
          abs(source.clus - sata.clus) == 1)
    }else if(bias == "downstream"){
      sata_hits <- sata_hits %>%
        dplyr::mutate(
          source.pos = clus_data[source.clus,]$clus.pos.mean,
          sata.pos = clus_data[sata.clus,]$clus.pos.mean,
          min.src.abund = clus_data[.$source.clus,]$min.abund,
          min.sat.abund = clus_data[.$sata.clus,]$min.abund,
          sum.src.abund = clus_data[.$source.clus,]$sum.abund,
          sum.sat.abund = clus_data[.$sata.clus,]$sum.abund,
          is.downstream = source.pos > sata.pos) %>%
        dplyr::filter(
          as.integer(min.src.abund) >= as.integer(min.sat.abund),
          sum.src.abund > sum.sat.abund,
          abs(source.clus - sata.clus) == 1)
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    if(nrow(sata_hits) > 0){
      clus.map <- GenomicRanges::findOverlaps(clus_ranges, red.sites)
      clus.list <- split(
        S4Vectors::subjectHits(clus.map), S4Vectors::queryHits(clus.map))
      
      if(bias == "upstream"){
        sata_hits <- sata_hits %>%
          dplyr::mutate(
            source.node = ifelse(
              sata_hits$is.upstream,
              sapply(clus.list[sata_hits$source.clus], dplyr::last),
              sapply(clus.list[sata_hits$source.clus], dplyr::first)),
            sata.node = ifelse(
              is.upstream,
              sapply(clus.list[sata.clus], dplyr::first),
              sapply(clus.list[sata.clus], dplyr::last)))
      }else if(bias == "downstream"){
        sata_hits <- sata_hits %>%
          dplyr::mutate(
            source.node = ifelse(
              sata_hits$is.downstream,
              sapply(clus.list[sata_hits$source.clus], dplyr::first),
              sapply(clus.list[sata_hits$source.clus], dplyr::last)),
            sata.node = ifelse(
              is.downstream,
              sapply(clus.list[sata.clus], dplyr::last),
              sapply(clus.list[sata.clus], dplyr::first)))
      }
      
      sata.edges <- with(sata_hits, vzip(source.node, sata.node))
    }else{
      sata.edges <- c()
    }
  }else{
    sata.edges <- c()
  }
  igraph::add_edges(graph, sata.edges)
}


break_connecting_source_paths <- function(red.sites, graph, bias){
  src_nodes <- sources(graph)
  sources_p_clus <- IRanges::IntegerList(split(
    src_nodes, igraph::clusters(graph)$membership[src_nodes]))
  clus_w_multi_sources <- sources_p_clus[S4Vectors::elementNROWS(sources_p_clus) > 1]
  
  if(length(clus_w_multi_sources) > 0){
    adj_pairs <- do.call(c, lapply(clus_w_multi_sources, function(x){
      lapply(1:(length(x)-1), function(i) c(x[i], x[i+1]))
    }))
    
    snk_nodes <- sinks(graph)
    
    edges_to_edit <- data.frame(
      src.node.i = unlist(adj_pairs)[
        IRanges::start(IRanges::IntegerList(adj_pairs)@partitioning)],
      src.node.j = unlist(adj_pairs)[
        IRanges::end(IRanges::IntegerList(adj_pairs)@partitioning)]) %>%
      dplyr::mutate(
        src.node.i.abund = as.numeric(red.sites[src.node.i]$abund),
        src.node.j.abund = as.numeric(red.sites[src.node.j]$abund),
        sink.node = IRanges::start(
          IRanges::findOverlapPairs(
            IRanges::IRanges(src.node.i, src.node.j),
            IRanges::IRanges(snk_nodes, width = 1))@second))
    
    # Identify the nodes adjacent to sinks between connected sources
    # then filter adjacent pairs to identify which edge should be 'clipped'.
    # Filtering based first on adjacent node distance (edges with greater
    # distance get clipped), then abundance (lower abund gets clipped), then
    # biasing on upstream edges over downstream (downstream is clipped for
    # tie breaking).
    
    if(bias == "upstream"){
      target_edges <- dplyr::bind_rows(lapply(
        seq_len(nrow(edges_to_edit)), function(i){
          sink <- edges_to_edit[i, "sink.node"]
          path <- unlist(igraph::all_simple_paths(
            igraph::as.undirected(graph),
            edges_to_edit[i, "src.node.i"],
            edges_to_edit[i, "src.node.j"]))
          pos <- which(path == sink)
          data.frame(
            sink = rep(sink, 2),
            adj.node = c(path[pos-1], path[pos+1])) })) %>%
        dplyr::mutate(
          sink.pos = GenomicRanges::start(red.sites[sink]),
          adj.pos = GenomicRanges::start(red.sites[adj.node]),
          adj.abund = red.sites[adj.node]$abund,
          nt.dist = abs(sink.pos - adj.pos),
          strand = as.character(GenomicRanges::strand(red.sites[sink])),
          is.upstream = ifelse(
            strand == "+", sink.pos < adj.pos, sink.pos > adj.pos)) %>%
        dplyr::group_by(sink) %>%
        dplyr::filter(nt.dist == max(nt.dist)) %>%
        dplyr::filter(adj.abund == min(adj.abund)) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, !is.upstream)) %>%
        dplyr::filter(keep) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else if(bias == "downstream"){
      target_edges <- dplyr::bind_rows(lapply(
        seq_len(nrow(edges_to_edit)), function(i){
          sink <- edges_to_edit[i, "sink.node"]
          path <- unlist(igraph::all_simple_paths(
            igraph::as.undirected(graph),
            edges_to_edit[i, "src.node.i"],
            edges_to_edit[i, "src.node.j"]))
          pos <- which(path == sink)
          data.frame(
            sink = rep(sink, 2),
            adj.node = c(path[pos-1], path[pos+1])) })) %>%
        dplyr::mutate(
          sink.pos = GenomicRanges::start(red.sites[sink]),
          adj.pos = GenomicRanges::start(red.sites[adj.node]),
          adj.abund = red.sites[adj.node]$abund,
          nt.dist = abs(sink.pos - adj.pos),
          strand = as.character(
            GenomicRanges::strand(red.sites[sink])),
          is.downstream = ifelse(
            strand == "+", sink.pos > adj.pos, sink.pos < adj.pos)) %>%
        dplyr::group_by(sink) %>%
        dplyr::filter(nt.dist == max(nt.dist)) %>%
        dplyr::filter(adj.abund == min(adj.abund)) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, !is.downstream)) %>%
        dplyr::filter(keep) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    break_edges <- with(target_edges, vzip(sink, adj.node))
    
    edge_ids_to_break <- igraph::get.edge.ids(
      graph, break_edges, directed = FALSE)
  }else{
    edge_ids_to_break <- c()
  }
  
  igraph::delete_edges(graph, edge_ids_to_break)
}

sources <- function(graph){
  if(!igraph::is_directed(graph)){
    message("Graph provided is not a directed graph.")
    srcs <- c()
  }else{
    srcs <- which(Matrix::colSums(
      igraph::get.adjacency(graph, sparse = TRUE)) == 0)
  }
  srcs
}


connect_adjacent_clusters <- function(red.sites, graph, gap, bias){
  src_nodes <- sources(graph)
  near_sources <- GenomicRanges::findOverlaps(
    red.sites[src_nodes],
    maxgap = gap - 1L,
    drop.self = TRUE
  )
  
  if(length(near_sources) > 0){
    # Identify sources of clusters within the largets satalite gap distance
    # and identify the directionality of the edge to create based first on
    # abundance (source will likely have greater abundance) and then by
    # upstream bias (more likely the origin site is upstream of drifting mapped
    # reads).
    if(bias == "upstream"){
      near_src_df <- data.frame(
        node.i = src_nodes[S4Vectors::queryHits(near_sources)],
        node.j = src_nodes[S4Vectors::subjectHits(near_sources)]) %>%
        dplyr::mutate(
          abund.i = red.sites[node.i]$abund,
          abund.j = red.sites[node.j]$abund,
          pos.i = GenomicRanges::start(red.sites[node.i]),
          pos.j = GenomicRanges::start(red.sites[node.j]),
          strand = as.character(GenomicRanges::strand(red.sites[node.i])),
          is.upstream = ifelse(strand == "+", pos.i < pos.j, pos.i > pos.j))
    }else if(bias == "downstream"){
      near_src_df <- data.frame(
        node.i = src_nodes[S4Vectors::queryHits(near_sources)],
        node.j = src_nodes[S4Vectors::subjectHits(near_sources)]) %>%
        dplyr::mutate(
          abund.i = red.sites[node.i]$abund,
          abund.j = red.sites[node.j]$abund,
          pos.i = GenomicRanges::start(red.sites[node.i]),
          pos.j = GenomicRanges::start(red.sites[node.j]),
          strand = as.character(GenomicRanges::strand(red.sites[node.i])),
          is.downstream = ifelse(strand == "+", pos.i > pos.j, pos.i < pos.j))
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    redundant_graph <- igraph::make_graph(
      edges = with(near_src_df, vzip(
        1:nrow(near_src_df),
        match(paste(node.i, node.j), paste(node.j, node.i)))))
    redundant_groups <- igraph::clusters(redundant_graph)$membership
    
    if(bias == "upstream"){
      near_src_df <- dplyr::mutate(
        near_src_df, redundant.grp = redundant_groups) %>%
        dplyr::filter(abund.i >= abund.j) %>%
        dplyr::group_by(redundant.grp) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, is.upstream)) %>%
        dplyr::filter(keep)
    }else if(bias == "downstream"){
      near_src_df <- dplyr::mutate(
        near_src_df, redundant.grp = redundant_groups) %>%
        dplyr::filter(abund.i >= abund.j) %>%
        dplyr::group_by(redundant.grp) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, is.downstream)) %>%
        dplyr::filter(keep)
    }
    
    edges_to_connect_near_srcs <- with(near_src_df, vzip(node.i, node.j))
  }else{
    edges_to_connect_near_srcs <- c()
  }
  
  igraph::add.edges(graph, edges_to_connect_near_srcs)
}


resolve_cluster_sources <- function(red.sites, graph, bias){
  src_nodes <- sources(graph)
  sources_p_clus <- IRanges::IntegerList(split(
    src_nodes, igraph::clusters(graph)$membership[src_nodes]))
  clus_w_multi_sources <- sources_p_clus[S4Vectors::elementNROWS(sources_p_clus) > 1]
  
  if(length(clus_w_multi_sources) > 0){
    if(bias == "upstream"){
      resolve_df <- data.frame(
        node = unlist(clus_w_multi_sources),
        clus = as.numeric(S4Vectors::Rle(
          values = seq_along(clus_w_multi_sources),
          lengths = S4Vectors::elementNROWS(clus_w_multi_sources)))) %>%
        dplyr::mutate(abund = red.sites[node]$abund) %>%
        dplyr::group_by(clus) %>%
        dplyr::mutate(
          top.abund = abund == max(abund),
          strand = as.character(GenomicRanges::strand(red.sites[node])),
          pos = GenomicRanges::start(red.sites[node])) %>%
        dplyr::group_by(clus, top.abund) %>%
        dplyr::mutate(
          grp.size = n(),
          is.upstream = ifelse(strand == "+", pos == min(pos), pos == max(pos)),
          src = ifelse(
            top.abund == TRUE & grp.size > 1, is.upstream, top.abund)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else if(bias == "downstream"){
      resolve_df <- data.frame(
        node = unlist(clus_w_multi_sources),
        clus = as.numeric(S4Vectors::Rle(
          values = seq_along(clus_w_multi_sources),
          lengths = S4Vectors::elementNROWS(clus_w_multi_sources)))) %>%
        dplyr::mutate(
          abund = red.sites[node]$abund) %>%
        dplyr::group_by(clus) %>%
        dplyr::mutate(
          top.abund = abund == max(abund),
          strand = as.character(GenomicRanges::strand(red.sites[node])),
          pos = GenomicRanges::start(red.sites[node])) %>%
        dplyr::group_by(clus, top.abund) %>%
        dplyr::mutate(
          grp.size = n(),
          is.downstream = ifelse(
            strand == "+", pos == max(pos), pos == min(pos)),
          src = ifelse(
            top.abund == TRUE & grp.size > 1, is.downstream, top.abund)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }
    
    src_nodes <- resolve_df[resolve_df$src == TRUE,]$node
    snk_nodes <- lapply(seq_along(clus_w_multi_sources), function(i){
      resolve_df[resolve_df$clus == i & resolve_df$src == FALSE,]$node
    })
    
    # Accomidates multiple sinks for singular sources in a cluster.
    resolve_edges <- do.call(c, lapply(seq_along(src_nodes), function(i){
      src <- src_nodes[i]
      snk <- snk_nodes[[i]]
      vzip(rep(src, length(snk)), snk) }))
  }else{
    resolve_edges <- c()
  }
  
  igraph::add.edges(graph, resolve_edges)
}


sinks <- function(graph){
  if(!igraph::is_directed(graph)){
    message("Graph provided is not a directed graph.")
    snks <- c()
  }else{
    snks <- which(Matrix::rowSums(
      igraph::get.adjacency(graph, sparse = TRUE)) == 0)
  }
  snks
}
