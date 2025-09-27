
.sf_target_pos_vec <- function(gr, mode){
  s <- as.character(strand(gr))
  pos <- integer(length(gr))
  if (mode == "intPositions") { pos[s=="+"] <- start(gr)[s=="+"]; pos[s=="-"] <- end(gr)[s=="-"]
  } else {                       pos[s=="+"] <- end(gr)[s=="+"]  ; pos[s=="-"] <- start(gr)[s=="-"] }
  pos
}

.sf_wmedian <- function(x, w){
  o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w)/sum(w); x[which(cw >= 0.5)[1]]
}

.sf_aggregate_stats <- function(df, w_unique, w_reads){
  key <- paste(df$start, df$end, sep=":")
  n_unique  <- tapply(key,      df$pos, function(v) length(unique(v)))
  sum_reads <- tapply(df$reads, df$pos, sum)
  pos <- as.integer(names(n_unique))
  sr  <- sum_reads[match(pos, as.integer(names(sum_reads)))]
  out <- data.frame(pos=pos,
                    n_unique=as.integer(n_unique),
                    sum_reads=as.numeric(sr),
                    stringsAsFactors=FALSE)
  out$score <- w_unique*out$n_unique + w_reads*out$sum_reads
  out[order(out$pos), ]
}

.sf_pick_better <- function(i,j, n_unique, sum_reads, score, tiebreak){
  if (tiebreak == "unique_then_reads"){
    if (n_unique[i] != n_unique[j]) return(if (n_unique[i] > n_unique[j]) i else j)
    if (sum_reads[i] != sum_reads[j]) return(if (sum_reads[i] > sum_reads[j]) i else j)
  }
  if (score[i] != score[j]) return(if (score[i] > score[j]) i else j)
  return(i)  # stable: pick left if truly tied
}

.sf_find_seeds <- function(pos, score, n_unique, sum_reads,
                           seed_min_score, seed_min_sep_bp, seed_min_prom, seed_valley_frac,
                           seed_local_radius_bp = NA_real_, seed_tiebreak = c("unique_then_reads","score")){
  seed_tiebreak <- match.arg(seed_tiebreak)
  n <- length(pos)
  if (n == 0L) return(integer(0))
  
  # Windowed maxima (optional)
  if (!is.na(seed_local_radius_bp)) {
    keep <- rep(FALSE, n)
    for (i in seq_len(n)) {
      if (score[i] < seed_min_score) next
      L <- max(1L, which.max(pos >= pos[i] - seed_local_radius_bp))
      R <- max(which(pos <= pos[i] + seed_local_radius_bp))
      win <- L:R
      mx <- max(score[win])
      if (score[i] == mx) keep[i] <- TRUE
    }
    cand <- which(keep)
    if (!length(cand)) return(integer(0))
    # enforce min separation (greedy left-to-right)
    seeds <- integer(0); last <- -Inf
    for (idx in cand){
      if (pos[idx] - last >= seed_min_sep_bp) { seeds <- c(seeds, pos[idx]); last <- pos[idx] }
    }
    return(seeds)
  }
  
  # Strict local maxima
  if (n == 1L) return(if (score[1] >= seed_min_score) pos else integer(0))
  is_peak <- logical(n)
  for (i in seq_len(n)){
    if (i == 1L)      is_peak[i] <- score[i] > score[i+1]
    else if (i == n)  is_peak[i] <- score[i] > score[i-1]
    else              is_peak[i] <- (score[i] >= score[i-1]) && (score[i] > score[i+1])
  }
  pk_idx <- which(is_peak & score >= seed_min_score)
  if (length(pk_idx) <= 1L) return(pos[pk_idx])
  
  # Prune adjacent by separation/valley with tiebreak (fragments first, then reads)
  keep <- rep(FALSE, length(pk_idx))
  i <- 1L
  while (i <= length(pk_idx)){
    cur <- pk_idx[i]; keep_i <- TRUE
    if (i < length(pk_idx)){
      nxt <- pk_idx[i+1]
      valley <- if ((nxt - cur) > 1L) min(score[(cur+1):(nxt-1)]) else min(score[cur], score[nxt])
      h1 <- score[cur]; h2 <- score[nxt]; hmax <- max(h1, h2)
      sep_ok <- (pos[nxt] - pos[cur]) >= seed_min_sep_bp
      valley_ok <- ((hmax - valley) >= seed_min_prom) || (valley <= (1 - seed_valley_frac)*hmax)
      if (!sep_ok || !valley_ok){
        better <- .sf_pick_better(cur, nxt, n_unique, sum_reads, score, seed_tiebreak)
        keep[i] <- (better == cur)
        i <- i + 2L
        next
      }
    }
    keep[i] <- keep_i; i <- i + 1L
  }
  pos[pk_idx[keep]]
}

.sf_initial_assign <- function(pos, centers, reach_bp){
  if (length(centers) == 0L) return(rep(NA_integer_, length(pos)))
  d <- abs(outer(pos, centers, "-"))
  j <- max.col(-d, ties.method="first")
  mind <- d[cbind(seq_along(pos), j)]
  j[mind > reach_bp] <- NA_integer_
  j
}

.sf_cluster_centers <- function(stats, clusters, center_method, w_unique, w_reads, seeds = integer(0)){
  keep <- !is.na(clusters)
  if (!any(keep)) return(list(centers=integer(0), clus_levels=integer(0), weights_tot=numeric(0)))
  stats$w <- w_unique*stats$n_unique + w_reads*stats$sum_reads
  clus_levels <- sort(unique(clusters[keep]))
  fac <- factor(clusters[keep], levels = clus_levels)
  
  if (center_method == "seed") {
    # centers are the seed positions indexed by label
    centers_vec <- sapply(clus_levels, function(lbl) seeds[lbl])
  } else if (center_method == "topscore") {
    centers_vec <- sapply(split(seq_len(nrow(stats))[keep], fac), function(ii){
      ii[which.max(stats$score[ii])] |> (\(k) stats$pos[k])()
    })
  } else if (center_method == "wmean") {
    centers_vec <- sapply(split(seq_len(nrow(stats))[keep], fac), function(ii){
      round(stats::weighted.mean(stats$pos[ii], stats$w[ii]))
    })
  } else { # wmedian (default)
    centers_vec <- sapply(split(seq_len(nrow(stats))[keep], fac), function(ii){
      round(.sf_wmedian(stats$pos[ii], stats$w[ii]))
    })
  }
  
  weights_tot <- sapply(split(stats$w[keep], fac), sum)
  list(centers=as.integer(centers_vec),
       clus_levels=clus_levels,
       weights_tot=as.numeric(weights_tot))
}

.sf_weighted_sigma <- function(x, w){
  if (length(x) <= 1L) return(0)
  if (all(w == 0)) w <- rep(1, length(x))
  mu <- sum(w*x)/sum(w)
  v  <- sum(w*(x - mu)^2)/sum(w)
  sqrt(max(v, 0))
}

.sf_gaussian_scores <- function(x, centers, priors, sigmas){
  L <- length(centers); if (L == 0L) return(matrix(nrow=length(x), ncol=0))
  logp <- log(pmax(priors, 1e-9))
  sapply(seq_len(L), function(j){ -0.5 * ((x - centers[j]) / sigmas[j])^2 + logp[j] })
}

.sf_apply_margin <- function(scores, home_idx, x, centers, reach_bp, margin){
  if (ncol(scores) == 0L) return(rep(NA_integer_, nrow(scores)))
  jbest <- max.col(scores, ties.method="first")
  s_best <- scores[cbind(seq_along(jbest), jbest)]
  tmp <- scores; tmp[cbind(seq_along(jbest), jbest)] <- -Inf
  j2nd  <- if (ncol(scores) >= 2L) max.col(tmp, ties.method="first") else rep(NA_integer_, length(jbest))
  s_2nd <- if (ncol(scores) >= 2L) tmp[cbind(seq_along(j2nd), j2nd)] else rep(-Inf, length(jbest))
  delta <- s_best - s_2nd
  
  within_best <- abs(x - centers[jbest]) <= reach_bp
  jbest[!within_best] <- NA_integer_
  
  low <- which(delta < margin)
  if (length(low)){
    home <- home_idx[low]
    ok_home <- !is.na(home) & abs(x[low] - centers[home]) <= reach_bp
    jbest[low[ok_home]] <- home[ok_home]
  }
  
  jbest
}

# --- Main ------------------------------------------------------------------

standardizeFragments <- function(
    gr,
    mode = c("intPositions","breakPositions"),
    w_unique = 1, 
    w_reads = 0.1,
    reach_bp = 6,
    center = c("seed", "wmedian","wmean","topscore"),
    assign_rule = c("gaussian", "nearest", "power"),
    label_outside_reach = c("na","nearest"),
    seed_min_score = 1,
    seed_min_sep_bp = 3,
    seed_min_prom = 5,
    seed_valley_frac = 0.30,
    seed_local_radius_bp = NA_real_,
    seed_tiebreak = c("unique_then_reads","score"),
    power_alpha = 1.0,
    posterior_margin = 0.8,
    gaussian_iters = 1,
    gaussian_lock_mu = TRUE,
    gaussian_sigma_mode = c("seed_reads", "cluster","fixed"),
    sigma_base = 2.0,
    sigma_reads_exp = 0.25,
    sigma_cap = 5.0,
    sigma_fixed = 2.0,
    sigma_floor = 1.0
){
  stopifnot(methods::is(gr, "GRanges"))
  if (length(gr) == 0L) return(gr)
  if (is.null(gr$reads)) stop("gr must have a numeric 'reads' column.")
  if (any(is.na(gr$reads))) stop("'reads' contains NA.")
  mode <- match.arg(mode)
  center <- match.arg(center)
  assign_rule <- match.arg(assign_rule)
  label_outside_reach <- match.arg(label_outside_reach)
  seed_tiebreak <- match.arg(seed_tiebreak)
  gaussian_sigma_mode <- match.arg(gaussian_sigma_mode)
  
  # Target positions & aggregated stats
  pos_vec <- .sf_target_pos_vec(gr, mode)
  df <- data.frame(pos=pos_vec, start=start(gr), end=end(gr),
                   reads=as.numeric(gr$reads), stringsAsFactors=FALSE)
  stats <- .sf_aggregate_stats(df, w_unique, w_reads)
  
  # Seeds
  seeds <- .sf_find_seeds(stats$pos, stats$score, stats$n_unique, stats$sum_reads,
                          seed_min_score, seed_min_sep_bp, seed_min_prom, seed_valley_frac,
                          seed_local_radius_bp = seed_local_radius_bp,
                          seed_tiebreak = seed_tiebreak)
  if (length(seeds) == 0L) seeds <- stats$pos[which.max(stats$score)]
  
  # Initial nearest assignment = "home" (strict reach)
  home <- .sf_initial_assign(stats$pos, seeds, reach_bp)
  
  # Build first centers / priors
  # If center="seed", centers are the seed positions indexed by cluster label
  if (center == "seed"){
    clus_levels <- sort(unique(home[!is.na(home)]))
    centers <- if (length(clus_levels)) as.integer(seeds[clus_levels]) else integer(0)
    # approximate priors from home membership
    tmp_cc <- .sf_cluster_centers(stats, home, "wmedian", w_unique, w_reads) # weights only
    priors <- tmp_cc$weights_tot
    names(priors) <- tmp_cc$clus_levels
  } else {
    cc <- .sf_cluster_centers(stats, home, center, w_unique, w_reads, seeds = seeds)
    centers <- cc$centers; clus_levels <- cc$clus_levels; priors <- cc$weights_tot
  }
  
  clusters <- home
  
  if (assign_rule == "nearest" || length(centers) == 0L){
    clusters <- home
    
  } else if (assign_rule == "power") {
    if (length(centers) > 0L){
      logp <- log(pmax(priors, 1e-9))
      scores <- sapply(seq_along(centers), function(j){
        -abs(stats$pos - centers[j]) + power_alpha*logp[j]
      })
      clusters <- .sf_apply_margin(scores, home, stats$pos, centers, reach_bp, margin = posterior_margin)
    }
    
  } else if (assign_rule == "gaussian") {
    iters <- max(1L, as.integer(gaussian_iters))
    clusters <- home
    
    # Precompute seed reads for sigma if needed
    seed_reads_map <- NULL
    if (gaussian_sigma_mode == "seed_reads"){
      seed_reads_map <- stats$sum_reads[match(seeds, stats$pos)]
      seed_reads_map[is.na(seed_reads_map)] <- 0
    }
    
    for (t in seq_len(iters)){
      keep <- !is.na(clusters); if (!any(keep)) break
      
      # Priors from current clustering (weights by unique+reads)
      stats$w <- w_unique*stats$n_unique + w_reads*stats$sum_reads
      fac <- factor(clusters[keep], levels = sort(unique(clusters[keep])))
      priors_vec <- sapply(split(stats$w[keep], fac), sum)
      
      # Centers: either locked to seeds, or recomputed per center method
      if (gaussian_lock_mu || center == "seed"){
        clus_levels <- sort(unique(clusters[keep]))
        centers <- as.integer(seeds[clus_levels])
      } else {
        cc <- .sf_cluster_centers(stats, clusters, center, w_unique, w_reads, seeds = seeds)
        centers <- cc$centers; clus_levels <- cc$clus_levels
      }
      if (length(centers) == 0L) break
      
      # Sigmas per cluster
      if (gaussian_sigma_mode == "fixed"){
        sigmas <- rep(max(sigma_floor, sigma_fixed), length(centers))
      } else if (gaussian_sigma_mode == "seed_reads"){
        # width guided by seed reads
        sr <- if (length(centers)) seed_reads_map[clus_levels] else numeric(0)
        sigmas <- pmax(sigma_floor, pmin(sigma_cap, sigma_base * (sr + 1)^sigma_reads_exp))
      } else { # "cluster" (robust within-cluster spread)
        fac2 <- factor(clusters, levels = clus_levels)
        sigmas <- sapply(seq_along(clus_levels), function(k){
          ii <- which(fac2 == clus_levels[k])
          if (!length(ii)) return(sigma_floor)
          x <- stats$pos[ii]; w <- stats$w[ii]
          s <- .sf_weighted_sigma(x, w)
          max(s, sigma_floor)
        })
      }
      
      # Align priors to clus_levels
      priors <- as.numeric(priors_vec[as.character(clus_levels)])
      
      # Gaussian scores + margin
      scores <- .sf_gaussian_scores(stats$pos, centers, priors, sigmas)
      clusters <- .sf_apply_margin(scores, home, stats$pos, centers, reach_bp, margin = posterior_margin)
    }
  }
  
  # ---- Robust, position-indexed mapping back to GRanges ----
  pos_to_cluster <- clusters
  center_by_label <- setNames(centers, as.character(clus_levels))
  pos_to_center <- rep(NA_integer_, length(stats$pos))
  good <- !is.na(pos_to_cluster)
  if (any(good)) pos_to_center[good] <- as.integer(center_by_label[ as.character(pos_to_cluster[good]) ])
  
  idx_map    <- match(pos_vec, stats$pos)
  cl_idx     <- ifelse(is.na(idx_map), NA_integer_, pos_to_cluster[idx_map])
  pos2center <- ifelse(is.na(idx_map), NA_integer_, pos_to_center[idx_map])
  
  # Optionally fill labels for beyond-reach points (coords unchanged)
  if (label_outside_reach == "nearest" && any(is.na(cl_idx)) && length(centers) > 0L) {
    dists <- abs(outer(stats$pos, centers, "-"))
    nearest_all <- max.col(-dists, ties.method="first")
    fill <- which(is.na(cl_idx) & !is.na(idx_map))
    if (length(fill)) cl_idx[fill] <- nearest_all[idx_map[fill]]
  }
  
  # Keep original coordinate if unassigned
  if (anyNA(pos2center)) {
    na_m <- is.na(pos2center)
    pos2center[na_m] <- pos_vec[na_m]
  }
  
  # ---- Write back one end (no NAs in start/end) ----
  out <- gr; s <- as.character(strand(out))
  if (mode == "intPositions") {
    i <- which(s == "+"); if (length(i)) start(out)[i] <- pos2center[i]
    i <- which(s == "-"); if (length(i)) end(out)[i]   <- pos2center[i]
  } else {
    i <- which(s == "+"); if (length(i)) end(out)[i]   <- pos2center[i]
    i <- which(s == "-"); if (length(i)) start(out)[i] <- pos2center[i]
  }
  if (any(start(out) > end(out))) warning("Some standardized starts exceed ends; please inspect those fragments.")
  
  out$std_pos     <- as.integer(pos2center)
  out$std_cluster <- as.integer(cl_idx)
  out$std_mode    <- mode
  out
}

inspect_standardization <- function(
    gr,
    mode = c("intPositions","breakPositions"),
    w_unique = 1, w_reads = 0.1,
    seed_min_score = 3,
    seed_min_sep_bp = 3,
    seed_min_prom = 5,
    seed_valley_frac = 0.30,
    reach_bp = 6
){
  mode <- match.arg(mode)
  pos_vec <- .sf_target_pos_vec(gr, mode)
  df <- data.frame(pos = pos_vec, start = start(gr), end = end(gr), reads = as.numeric(gr$reads))
  stats <- .sf_aggregate_stats(df, w_unique, w_reads)
  seeds <- .sf_find_seeds(stats$pos, stats$score, stats$n_unique, stats$sum_reads,
                          seed_min_score, seed_min_sep_bp, seed_min_prom, seed_valley_frac)
  if (length(seeds) == 0L) seeds <- stats$pos[which.max(stats$score)]
  clus <- .sf_initial_assign(stats$pos, seeds, reach_bp)
  data.frame(stats, is_seed = stats$pos %in% seeds, cluster_home = clus, stringsAsFactors = FALSE)
}

