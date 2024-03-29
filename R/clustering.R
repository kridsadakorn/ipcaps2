


#' (Internal function) Perform the clustering process of IPCAPS2
#'
#' @param dataframe A data frame containing \code{raw.data} (matrix or data
#' frame), \code{label} (vector), and \code{index} (vector). \code{raw.data}
#' represents a matrix of subset of input data. \code{label} represents a vector
#' of labels for all rows of \code{raw.data}. \code{index} represents a vector
#' of indexes that are selected from the original input data.
#' @param node An integer representing the current node number which is being
#' analyzed.
#' @param result.dir An output directory
#' @param threshold A threshold or a cutoff to stop the IPCAPS2 process. Also see
#' \code{\link{ipcaps2}} (the parameter \code{threshold}).
#' @param min.fst A number represents a cutoff for mininum Fst value.
#' @param method A clustering method selected from the \code{\link{ipcaps2}}
#' function. See \code{\link{ipcaps2}} for available methods.
#' @param min.in.group A integer represents a minimum number of group members.
#' @param datatype To specify whether the input data are 'snp' or other type.
#' Defalut = 'snp'.
#' @param nonlinear (Unimplemented) To specify whether linear or non-linear
#' method is used for IPCAPS2 analysis. If TRUE, non-linear method is used,
#' otherwise linear method is used. Default = FALSE.
#' @param no.rep To specify a  number of time to replicate the internal clustering. Default = 5,
#' @param cutoff.confident To specify a cutoff for confident values. Default = 0.5.
#'
#' @return A list containing \code{status}, \code{node}, and \code{new.index}
#' as explained below:
#' \itemize{
#' \item \code{$status} is either \code{0} representing that the criterion is
#' not met, or \code{1} representing that the criterion is met.
#' \item \code{$node} is an integer representing the current node number which
#' is being analyzed.
#' \item \code{$new.index} is a list of vectors containing a subset of indexes
#' split from \code{dataframe$index} according to a clustering result.
#' }
#'
#' @include check.stopping.R
#' @include clustering.mode.R
#'

clustering <-
  function(dataframe,
           node = 1,
           result.dir,
           threshold,
           min.fst,
           method = "mix",
           min.in.group = 20,
           datatype = "snp",
           nonlinear = FALSE,
           no.rep = 5,
           cutoff.confident = 0.5)
  {
    # the function for check which nodes are linked a purpose to use with
    # Fst matrix mymatrix is the Fst matrix of all pairs of population X is
    # the index of column to check, 1 to maximum number of column of
    # mymatrix root needs to be set as c(0)
    recursion.linked.node <- function(mymatrix, X, root, min.fst)
    {
      # fine-level structure, define Fst > 0.0008
      idx <- which(mymatrix[, X] <= min.fst)
      for (i in root)
      {
        del <- which(idx == i)
        if (length(del) > 0)
        {
          idx <- idx[-del]
        }
      }
      if (length(idx) <= 0)
      {
        return(X)
      } else if (length(idx) == 1)
      {
        return(c(X,
                 recursion.linked.node(mymatrix, idx[1], c(X), min.fst)))
      } else
      {
        tmp <- c(X)
        for (i in idx)
        {
          del <- which(idx == i)
          mark <- idx[-del]
          tmp <- c(tmp,
                   recursion.linked.node(mymatrix, i, c(X, mark),
                                         min.fst))
        }
        return(tmp)
      }
    }

    silence.mode <- NULL

    load(file.path(result.dir, "RData", "condition.RData"))

    index <- dataframe$index
    label <- dataframe$label
    confident.list <- dataframe$confident.list

    if (length(index) < min.in.group)
    {
      if (!silence.mode)
        cat(
          paste0(
            "Node ",
            node,
            ": A number of node is lower than the",
            " minimum number (",
            min.in.group,
            "), therefore split was not performed\n"
          )
        )
      file.name <-
        file.path(result.dir, "RData", paste0("node", node,
                                              ".RData"))
      PCs <- NA
      eigen.value <- NA
      eigen.fit <- NA
      save(
        PCs,
        eigen.value,
        eigen.fit,
        threshold,
        label,
        index,
        confident.list,
        file = file.name,
        compress = "bzip2"
      )
      # case of status = 1, no split, stopping criteria are met
      ret <- list(status = 1)
    } else
    {
      X <- as.matrix(dataframe$raw.data)

      if (!silence.mode)
        cat(paste0("Node ", node, ": Reducing matrix\n"))

      res <- NA

      if (nonlinear == FALSE)
      {
        PCobj <-
          cal.pc.linear(X,
                        no.pc = min.in.group,
                        data.type = datatype,
                        PCscore = FALSE)
        PCs <- PCobj$PC
        eigen.value <- PCobj$evalue

        res <- check.stopping(eigen.value, threshold)

        eigen.fit <- res$eigen.fit
        no.significant.PC <- res$no.significant.PC
        threshold <- res$threshold

      } else
      {
        # for Non-linear PCs
        # PC.non = cal.PC.non.linear(X, no.pc =
        # min.in.group, method='kpca')
        # PC.non = cal.PC.non.linear(X, no.pc =
        # min.in.group, method='other')
        # PCs=PC.non$PC eigen.value =
        # PC.non$evalue res = list('status' = 0)
        # eigen.fit = NA
        # no.significant.PC = 3
        # threshold = threshold
      }


      #### Save result files
      file.name <- file.path(result.dir, "RData",
                             paste0("node", node,
                                    ".RData"))
      save(
        PCs,
        eigen.value,
        eigen.fit,
        threshold,
        label,
        index,
        confident.list,
        no.significant.PC,
        file = file.name,
        compress = "bzip2"
      )
      if (!silence.mode)
        cat(
          paste0(
            "Node ",
            node,
            ": EigenFit = ",
            eigen.fit,
            ", Threshold = ",
            threshold,
            ", no. significant PCs = ",
            no.significant.PC,
            "\n"
          )
        )
      #### End of saving result files returned status from check stopping
      #### criteria, status = 1 means stop
      if (res$status == 0)
      {
        # multiple testing for split

        start.time <- Sys.time()

        running.method <- method
        res.cluster <- clustering.mode(
          node = node,
          work.dir = result.dir,
          method = running.method,
          silence.mode = silence.mode,
          no.rep = no.rep,
          cutoff.confident = cutoff.confident
        )

        #new.index <- list()
        #confident.list <- list()
        new.local.index <- list()
        length.index <- c()
        for (i in unique(res.cluster$cluster))
        {
          new.local.index[[i]] <- which(res.cluster$cluster == i)
          length.index <-
            c(length.index, length(new.local.index[[i]]))

          #new.index[[i]] <- index[which(res.cluster$cluster == i)]
          #confident.list[[i]] <- res.cluster$cluster.confident[new.local.index[[i]]]
          #length.index <- c(length.index, length(new.index[[i]]))
        }


        # Check for Fst at least 0.0008, only for SNP
        #no.cluster <- length(new.index)
        no.cluster <- length(new.local.index)

        if ((no.cluster > 1) && (datatype == "snp"))
        {
          X <- as.matrix(dataframe$raw.data)
          fst.mat <-
            matrix(rep(1, no.cluster * no.cluster),
                   byrow = TRUE,
                   nrow = no.cluster)
          for (i in seq(1, (no.cluster - 1)))
          {
            for (j in (i + 1):no.cluster)
            {
              if (!silence.mode)
                cat(paste0(
                  "Node ",
                  node,
                  ": Checking Fst of group ",
                  i,
                  " and group ",
                  j,
                  "\n"
                ))
              fst <-
                fst.hudson(X,
                           new.local.index[[i]],
                           new.local.index[[j]])
              fst.mat[i, j] <- fst
              fst.mat[j, i] <- fst
            }
          }

          current.group.no <- 0
          tmp.index <- list()
          length.tmp <- c()
          marked.as.group <- rep(0, no.cluster)
          for (i in seq(1, no.cluster))
          {
            if (marked.as.group[i] == 0)
            {
              current.group.no <- current.group.no + 1
              link.node <- recursion.linked.node(fst.mat, i, c(0),
                                                 min.fst)
              tmp.vec <- c()
              for (j in link.node)
              {
                #tmp.vec <- c(tmp.vec, new.index[[j]])
                tmp.vec <- c(tmp.vec, new.local.index[[j]])
                marked.as.group[j] <- current.group.no
              }

              tmp.index[[current.group.no]] <- sort(tmp.vec)
              length.tmp <- c(length.tmp, length(tmp.vec))
            }
          }

          #new.index <- tmp.index
          new.local.index <- tmp.index
        }

        save(res.cluster, new.local.index, i, dataframe, file = "/Users/kris/Desktop/Debug.RData")

        new.index <- list()
        confident.list <- list()
        no.cluster <- length(new.local.index)
        for (i in seq(1, no.cluster))
        {
          #new.index[[i]] <- index[new.local.index[[i]]] #index is modified somewhere else
          new.index[[i]] <- dataframe$index[new.local.index[[i]]]
          confident.list[[i]] <-
            res.cluster$confident[new.local.index[[i]]]
        }


        end.time <- Sys.time()
        time.taken <- end.time - start.time
        if (!silence.mode)
          cat(paste0("Node ", node, " times took ",
                     time.taken, " secs\n"))

        if (length(new.index) > 1)
        {
          if (!silence.mode)
            cat(
              paste0(
                "Node ",
                node,
                ": ",
                length(index),
                " individuals were splitted into: ",
                length(new.index),
                " groups\n"
              )
            )
          ret <- list(
            status = 0,
            node = node,
            new.index = new.index,
            confident.list = confident.list
          )

        } else
        {
          if (!silence.mode)
            cat(
              paste0(
                "Node ",
                node,
                ": No split was performed because only",
                " one cluster detected\n"
              )
            )
          ret <- list(status = 1)
        }
      } else
      {
        # case of status = 1, no split, stopping criteria are met
        if (!silence.mode)
          cat(paste0(
            "Node ",
            node,
            paste0(
              ": No split was performed because EigenFit ",
              "is lower than ",
              "threshold\n"
            )
          ))
        ret <- list(status = 1)
      }
    }

    return(ret)
  }
