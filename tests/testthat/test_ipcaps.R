context("test ipcaps2")

test_that("test ipcaps2",{

  BED.file <-
    system.file("extdata", "ipcaps_example.bed", package = "IPCAPS2")
  LABEL.file <-
    system.file("extdata", "ipcaps_example_individuals.txt.gz",
                package = "IPCAPS2")
  # Test data.type = 'snp'
  res <-
    ipcaps2(
      bed = BED.file,
      label.file = LABEL.file,
      lab.col = 2,
      out = tempdir(),
      max.thread = 1,
      silence = TRUE,
      no.rep = 1
    )

  expect_type(res, "list")
  expect_length(res, 2)

  unlink(res$output.dir, recursive = TRUE)

  # Test data.type = 'linear'
  res <-
    ipcaps2(
      bed = BED.file,
      label.file = LABEL.file,
      lab.col = 2,
      out = tempdir(),
      data.type = 'linear',
      max.thread = 2,
      silence = TRUE,
      no.rep = 1
    )

  expect_type(res, "list")
  expect_length(res, 2)

  unlink(res$output.dir, recursive = TRUE)
})
