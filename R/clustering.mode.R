

#' (Internal function) Select a clustering method to be used for the IPCAPS2
#' process.
#'
#' @param node An integer representing the current node number which is being
#' analyzed.
#' @param work.dir A working directory.
#' @param method A clustering method selected from the \code{\link{ipcaps2}}
#' function. See \code{\link{ipcaps2}} for available methods.
#' @param silence.mode To enable or disable silence mode. If silence mode is
#' enabled, the fuction is processed without printing any message on the screen,
#' and it is slightly faster. Default = TRUE.
#' @param seed To specify a seed number for random generator. Default = NA,
#' which means that a seed number is automatically chose.
#' @param no.rep To specify a  number of time to replicate the internal clustering. Default = 5,
#' @param cutoff.confident To specify a cutoff for confident values. The confident values
#' are calculated based the clustering results of all replicates. The confident values
#' represent the best of average values (agreed clusters) among the member of clusters.
#' Default = 0.5.
#'
#' @return A vector of cluster assignment, for which cluster each individual
#' belongs.
#'
#' @import stats
#' @import fpc
#' @import LPCM
#' @importFrom apcluster apcluster negDistMat
#' @importFrom Rmixmod mixmodCluster mixmodPredict
#' @import KRIS
#' @importFrom methods new
#'
#' @seealso \code{\link{ipcaps2}}

clustering.mode <-
  function(node,
           work.dir,
           method,
           silence.mode = FALSE,
           seed = NULL,
           no.rep = 5,
           cutoff.confident = 0.5)
  {
    start.time <- Sys.time()
    if (!silence.mode)
      cat(paste0("Node ", node, ": Start clustering\n"))

    PCs <- NULL
    no.significant.PC <- NULL

    # load experiment condition
    file.name <-
      file.path(work.dir, "RData", paste0("node", node, ".RData"))
    if (!silence.mode)
      cat(paste0("Node ", node, ": Loading ", file.name, "\n"))
    load(file = file.name)

    all_cluster = NULL
    for (t in seq(1, no.rep)) {
      if (method == "clara")
      {
        cl <- pamk(PCs, krange = 2:2, usepam = FALSE)
        cluster <- cl$pamobject$clustering
      } else if (method == "pam")
      {
        cl <- pamk(PCs, krange = 2:2, usepam = TRUE)
        cluster <- cl$pamobject$clustering
      } else if (method == "mixmod")
      {
        subPCs <- as.data.frame(PCs[, seq(1, no.significant.PC)])
        newplan <- new("Strategy", seed = seed)
        mmc <- mixmodCluster(
          data = subPCs,
          nbCluster = seq(1, 3),
          criterion = c("BIC", "ICL", "NEC"),
          models = NULL,
          strategy = newplan
        )
        pmm <-
          mixmodPredict(data = subPCs,
                        classificationRule = mmc["bestResult"])
        cluster <- pmm["partition"]
      } else if (method == "meanshift")
      {
        subPCs <- as.data.frame(PCs[, seq(1, no.significant.PC)])
        fit <- ms(subPCs, h = 0.095, plot = 0)
        cluster <- fit$cluster.label
      } else if (method == "apcluster")
      {
        subPCs <- as.data.frame(PCs[, seq(1, no.significant.PC)])
        ap1 <-
          apcluster(negDistMat(r = 5), subPCs, q = 0.001)
        # ap2 <- aggExCluster(x=ap1)
        cluster <- rep(0, dim(subPCs)[1])
        no.cluster <- length(ap1)
        for (i in seq(1, no.cluster))
        {
          cluster[ap1@clusters[[i]]] <- i
        }
      } else if (method == "hclust")
      {
        dis <- dist(PCs)
        hc <- hclust(dis, method = "complete")
        cluster <- cutree(hc, k = 2)
      } else if (method == "rubikclust")
      {
        cluster <- rubikclust(PCs[, seq(1, 3)], min.space = 0.15)
      } else
      {
        # default Mixed clustering methods
        cluster <- rubikclust(PCs[, seq(1, 3)])
        if (length(unique(cluster)) == 1)
        {
          # swith to mixmod
          subPCs <-
            as.data.frame(PCs[, seq(1, no.significant.PC)])
          newplan = new("Strategy", seed = seed)
          mmc <- mixmodCluster(
            data = subPCs,
            nbCluster = seq(1, 3),
            criterion = c("BIC", "ICL", "NEC"),
            models = NULL,
            strategy = newplan
          )
          pmm <- mixmodPredict(data = subPCs,
                               classificationRule = mmc["bestResult"])
          cluster <- pmm["partition"]
        }
      }
      all_cluster = cbind(all_cluster, cluster)
    }

    n = dim(all_cluster)[2]
    m = dim(all_cluster)[1]
    con_mat = matrix(rep(0, m * m), nrow = m, ncol = m)
    for (i in seq(1, n)) {
      groups = unique(all_cluster[, i])
      for (j in groups) {
        idx = which(all_cluster[, i] == j)
        for (k in idx) {
          for (h in idx) {
            con_mat[k, h] = con_mat[k, h] + 1
          }
        }
      }
    }
    ave_con_mat = con_mat / n
    filtered_mat = ave_con_mat
    filtered_mat[which(filtered_mat < cutoff.confident)] =  0

    group = rep(0, m)
    cluster.confident = rep(0, m)
    groupidx = 0

    for (i in seq(1, m)) {
      if (group[i] == 0) {
        groupidx = groupidx + 1
        group[i] = groupidx
        cluster.confident[i] = max(filtered_mat[, i])
        for (j in seq(1, m)) {
          if (filtered_mat[j, i] > 0) {
            if (group[j] == 0) {
              group[j] = groupidx
              cluster.confident[j] = filtered_mat[j, i]
            } else{
              if (cluster.confident[j] < filtered_mat[j, i]) {
                group[j] = groupidx
                cluster.confident[j] = filtered_mat[j, i]
              }
            }
          }
        }
      } else{
        if (max(filtered_mat[, i]) == 0) {
          groupidx = groupidx + 1
          group[i] = groupidx
          cluster.confident[i] = 0
        }
      }
    }

    cluster <- group
    res <-
      list('cluster' = cluster, 'confident' = cluster.confident)

    end.time <- Sys.time()
    if (!silence.mode)
      cat(paste0("Node ", node, ": done for clustering\n"))
    if (!silence.mode)
      cat(paste0("Node ", node, ": processed in ",
                 format(end.time - start.time)))
    return(res)

  }
