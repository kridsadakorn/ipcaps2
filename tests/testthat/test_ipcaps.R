context("test clustering.mode")

test_that("test clustering.mode",{

  BED.file <- system.file("extdata","ipcaps_example.bed",package="IPCAPS2")
  LABEL.file <- system.file("extdata","ipcaps_example_individuals.txt.gz",package="IPCAPS2")

  res <- ipcaps2(bed=BED.file,label.file=LABEL.file,lab.col=2,out=tempdir())

  expect_type(res, "list")
  expect_length(res, 2)

  unlink(res$output.dir, recursive = TRUE)
})
