#' (Internal function) Perform the iterative process for each node
#'
#' @param node An integer representing the current node number which is being
#' analyzed.
#' @param work.dir A working directory.
#' @param no.rep To specify a  number of time to replicate the internal clustering. Default = 5,
#' @param cutoff.confident To specify a cutoff for confident values. Default = 0.5.
#'
#' @return \code{NULL}
#'
#' @include clustering.R
#'
process.each.node <- function(node,
                              work.dir,
                              no.rep = 5,
                              cutoff.confident = 0.5)
{
  start.time <- Sys.time()

  label <- NULL
  raw.data <- NULL
  threshold <- NULL
  min.fst <- NULL
  method <- NULL
  min.in.group <- NULL
  datatype <- NULL
  nonlinear <- NULL
  result.dir <- NULL
  silence.mode <- NULL

  # load tree list, no need to lock file, just reading info from the file
  file.name <- file.path(work.dir, "RData", "condition.RData")
  load(file = file.name)
  file.name <- file.path(work.dir, "RData", "tree.RData")
  if (!silence.mode)
    cat(paste0("Node ", node, ": Loading ", file.name, "\n"))

  e <- globalenv()$myenv

  if (isTRUE(e$global.lock))
  {
    Sys.sleep(runif(1, min = 0, max = 2))
  } else
  {
    e$global.lock <- TRUE
    load(file = file.name)
    e$global.lock <- FALSE
  }

  if (!silence.mode)
    cat(paste0("Node ", node, ": Start the process\n"))

  which.row <- which(tree$node == node)
  # Status -1 = deleted node due to error
  if (tree$status[which.row[1]] != -1)
  {
    # load experiment condition
    file.name <- file.path(work.dir, "RData", "condition.RData")
    if (!silence.mode)
      cat(paste0("Node ", node, ": Loading ", file.name, "\n"))
    load(file = file.name)


    # load raw data
    file.name <- file.path(work.dir, "RData", "rawdata.RData")
    if (!silence.mode)
      cat(paste0("Node ", node, ": Loading ", file.name, "\n"))
    load(file = file.name)
    ref.label <- label

    # load node data
    file.name <-
      file.path(work.dir, "RData", paste0("node", node,
                                          ".RData"))
    if (!silence.mode)
      cat(paste0("Node ", node, ": Loading ", file.name, "\n"))
    load(file = file.name)
    ref.index <- index
    ref.confident.list <- confident.list

    dataframe <- list(
      raw.data = raw.data[ref.index, ],
      label = ref.label[ref.index],
      index = ref.index,
      confident.list = ref.confident.list
    )

    if (!silence.mode)
      cat(paste0("Node ", node, ": Check for splitting\n"))

    res <- clustering(
      dataframe,
      node,
      work.dir,
      threshold,
      min.fst,
      method,
      min.in.group,
      datatype,
      nonlinear,
      no.rep = no.rep,
      cutoff.confident = cutoff.confident
    )
    # list.new.node = data.frame()
    if (!silence.mode)
      cat(paste0("Node ", node, ": Return status ", res$status, "\n"))
    if (res$status == 0)
    {
      if (!silence.mode)
        cat(paste0("Node ", node, ": Split to sub-nodes\n"))

      # serial loop
      for (i in seq(1, length(res$new.index)))
      {
        file.name <- file.path(work.dir, "RData", "tree.RData")
        load(file = file.name)
        new.node <- max(tree$node) + 1

        tree <- rbind(tree,
                      c(
                        node = new.node,
                        parent.node = node,
                        status = 0
                      ))

        if (isTRUE(e$global.lock))
        {
          Sys.sleep(runif(1, min = 0, max = 2))
        } else
        {
          e$global.lock <- TRUE
          save(tree, res,
               file = file.name,
               compress = "bzip2")
          e$global.lock <- FALSE
        }


        index <- res$new.index[[i]]
        confident.list <- res$confident.list[[i]]
        file.name <- file.path(result.dir,
                               "RData",
                               paste0("node",
                                      new.node, ".RData"))
        if (!silence.mode)
          cat(paste0("Node ", node, ": Saving ", file.name, "\n"))
        save(index,
             confident.list,
             file = file.name,
             compress = "bzip2")
      }


      # load tree list
      file.name <- file.path(work.dir, "RData", "tree.RData")
      if (!silence.mode)
        cat(paste0("Node ", node, ": Loading ", file.name, "\n"))
      load(file = file.name)

      which.row <- which(tree$node == node)
      # Status -1 = deleted node due to error
      if (tree$status[which.row[1]] != -1)
      {
        tree$status[which.row[1]] <- 2

        # colnames(list.new.node) = c('node','parent.node','status')
        # tree =
        # rbind(tree,list.new.node)

        file.name <-
          file.path(result.dir, "RData", "tree.RData")
        if (!silence.mode)
          cat(paste0("Node ", node, ": Updating ", file.name, "\n"))

        if (isTRUE(e$global.lock))
        {
          Sys.sleep(runif(1, min = 0, max = 2))
        } else
        {
          e$global.lock <- TRUE
          save(tree,
               file = file.name,
               compress = "bzip2")
          e$global.lock <- FALSE
        }

      }


    } else if (res$status == 1)
    {
      # case of status = 5, no split, stopping criteria are met
      if (!silence.mode)
        cat(paste0("Node ", node,
                   ": No split was performed, Status=1\n"))

      file.name <-
        file.path(result.dir, "RData", "leafnode.RData")
      load(file = file.name)
      leaf.node <- c(leaf.node, node)
      if (!silence.mode)
        cat(paste0("Node ", node, ": Updating ", file.name, "\n"))
      save(leaf.node, file = file.name, compress = "bzip2")


      # load tree
      file.name <- file.path(work.dir, "RData", "tree.RData")
      if (!silence.mode)
        cat(paste0("Node ", node, ": Loading ", file.name, "\n"))
      load(file = file.name)

      which.row <- which(tree$node == node)
      tree$status[which.row[1]] <- 2

      file.name <-
        file.path(result.dir, "RData", "tree.RData")
      if (!silence.mode)
        cat(paste0("Node ", node, ": Updating ", file.name, "\n"))

      if (isTRUE(e$global.lock))
      {
        Sys.sleep(runif(1, min = 0, max = 2))
      } else
      {
        e$global.lock <- TRUE
        save(tree, file = file.name, compress = "bzip2")
        e$global.lock <- FALSE
      }
    }
  }

  end.time <- Sys.time()
  if (!silence.mode)
    cat(paste0("Node ", node, ": Done!\n"))
  if (!silence.mode)
    cat(paste0("Node ", node, ": processed in ", format(end.time -
                                                          start.time)))

  invisible(NULL)
}
