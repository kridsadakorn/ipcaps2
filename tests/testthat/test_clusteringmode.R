context("test clustering.mode")

test_that("test clustering.mode", {
  PCs <- matrix(runif(300), byrow = TRUE, ncol = 3)
  node <- 0
  no.significant.PC <- 3
  workdir <- tempdir()

  #Interrupt the regular workflow in order to test the function
  targetdir <- file.path(workdir, "RData")
  if (!dir.exists(targetdir)) {
    dir.create(targetdir)
  }

  file1 <- file.path(workdir, "RData", "node0.RData")
  save(PCs,
       node,
       no.significant.PC,
       file = file1,
       compress = "bzip2")

  res <-
    clustering.mode(node, workdir, "clara", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  res <- clustering.mode(node, workdir, "pam", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  res <-
    clustering.mode(node, workdir, "mixmod", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  res <-
    clustering.mode(node, workdir, "meanshift", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  res <-
    clustering.mode(node, workdir, "apcluster", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  res <-
    clustering.mode(node, workdir, "hclust", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  res <-
    clustering.mode(node, workdir, "rubikclust", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  res <- clustering.mode(node, workdir, "mix", silence.mode = TRUE)
  expect_length(res, 2)
  expect_type(res, "list")
  expect_length(res$cluster, 100)
  expect_type(res$cluster, "double")
  expect_length(res$confident, 100)
  expect_type(res$confident, "double")

  file.remove(file1)
  file.remove(targetdir)
})
