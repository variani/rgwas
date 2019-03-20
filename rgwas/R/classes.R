#----------------------
# Class Declatations
#----------------------

#' S3 class MFMR.
#' @name MFMR
#' @rdname MFMRClass
#' @exportClass MFMR

#----------------------
# Method Declatations
#----------------------

#' @export plotYb
plotYb <- function(x, ...) UseMethod("plotYb")

#' @export plotYq
plotYq <- function(x, ...) UseMethod("plotYq")

#' @export plotK
plotK <- function(x, ...) UseMethod("plotK")

#------------------------------
# Wrapper to mfmr (class MFMR)
#------------------------------

call_mfmr <- function(Yb, Yq, G, X = NULL, K, nrun = 10, 
  thr = 0.95, ...)
{
  # catch the call & run mfmr
  mc <- match.call()
  
  env <- parent.frame(1)
  env$mfmr <- mfmr
  mc[[1]] <- quote(mfmr)
  out <- eval(mc, env)
  
  # add more output slots
  clust <- apply(out$pmat, 1, function(x) {
    k <- which.max(x)
    p <- x[k]
    `if`(p > thr, k, 0)
  })
  
  # return
  out <- c(out,
    list(call = mc,
      # input
      N = nrow(Yb), B = ncol(Y), P = ncol(Y), S = ncol(G), Q = ncol(X),
      K = K, nrun = nrun,
      # output
      clust = clust))
    
  oldClass(out) <- "MFMR"
  out
}

#----------------------
# Plot Methods
#----------------------

plotYb.MFMR <- function(x, data = c("Yb", "Gb"), ...)
{
  data <- match.arg(data)
  
  Yb <- switch(data,
    "Yb" = {
      var_Yb <- as.list(x$call)[["Yb"]] 
      get(as.character(var_Yb))
    },
    "Gb" = {
      stop("not implemented yet")
      var_Yb <- as.list(x$call)[["G"]] 
      Yb <- get(as.character(var_Yb))
    },
    stop("switch on `data`"))
    
  tab <- lapply(colnames(Yb), function(t) {
    y <- Yb[, t]
    
    clusters <- sort(unique(x$clust))
    prev <- lapply(clusters, function(k) {
      yk <- `if`(k == 0, y, y[x$clust == k])
      tibble(clust = k, prev = sum(yk) / length(yk), prev0 = sum(y) / length(y)) # prevalence
    }) %>% bind_rows
      
    mutate(prev, trait_b = t)
  }) %>% bind_rows
  
  # prepare data for plotting
  tab <- mutate(tab, 
    clust = factor(clust),
    trait_b = factor(trait_b, levels = colnames(Yb)))

  # plot
  p <- ggplot(tab, aes(clust, prev, fill = clust)) + geom_bar(stat = "identity", position = "identity") + geom_hline(aes(yintercept = prev0), linetype = 3) + facet_wrap(~ trait_b, scales = "free", nrow = 1)
  
  p <- p + theme(legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  
  # return
  p
}

plotYq.MFMR <- function(x, data = c("Yq", "G"), ...)
{
  data <- match.arg(data)
  
  Y <- switch(data,
    "Yq" = {
      var_Y <- as.list(x$call)[["Yq"]]
      get(as.character(var_Y))
    },
    "G" = {
      var_Y <- as.list(x$call)[["G"]] 
      get(as.character(var_Y))
    },
    stop("switch on `data`"))
  
  tab <- lapply(colnames(Y), function(t) {
    y <- Y[, t]
    # the case of Intercept
    scale <- (sd(y) > 0)
    y <- scale(y, center = TRUE, scale = scale)
    
    clusters <- sort(unique(x$clust))
    stats <- lapply(clusters, function(k) {
      yk <- `if`(k == 0, y, y[x$clust == k])
      tibble(clust = k, 
        mean = mean(yk), q95 = quantile(yk, 0.95), q5 = quantile(yk, 0.05),
        mean0 = mean(y))
    }) %>% bind_rows
      
    mutate(stats, trait_q = t)
  }) %>% bind_rows
  
  # prepare data for plotting
  tab <- mutate(tab, 
    clust = factor(clust),
    trait_q = factor(trait_q, levels = colnames(Y)))

  # plot
  p <- ggplot(tab, aes(clust, mean, color = clust)) + geom_point() + geom_linerange(aes(ymin = q5, ymax = q95)) + geom_hline(aes(yintercept = mean0), linetype = 3) + facet_wrap(~ trait_q, scales = "free", nrow = 1)
  
  p <- p + theme(legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  
  # return
  p
}

# https://stackoverflow.com/questions/48184645/how-can-i-put-the-labels-outside-of-piechart?noredirect=1&lq=1
# https://stackoverflow.com/questions/47752037/pie-chart-with-ggplot2-with-specific-order-and-percentage-annotations
plotK.MFMR <- function(x, ...)
{
  K <- x$K
  clust <- x$clust
 
  clusters <- sort(unique(x$clust))
  
  cnames <- clusters
  names(cnames) <- as.character(clusters)
  
  tab <- lapply(clusters, function(k) 
      tibble(k = k, kname = cnames[[as.character(k)]], 
        num = sum(clust == k), prop = sum(clust == k)/length(clust))) %>% 
    bind_rows
  
  props <- tab$prop
  nums <- tab$num
  
  # prepare data for plotting
  tab <- mutate(tab, 
    lab = factor(k, levels = clusters, 
      labels = paste0(cnames, ": ", format(num, big.mark = ","), " (", round(100*props, 0), "%)")))

  # plot
  p <- ggplot(tab, aes("", prop, fill = lab)) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y")
  
  #p <- p + geom_text(aes(label = paste0(round(100*prop, 0), "%")), position = position_stack(vjust = 0.5)) 
  
  p <- p + theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank())
  
  # return
  p             
}
