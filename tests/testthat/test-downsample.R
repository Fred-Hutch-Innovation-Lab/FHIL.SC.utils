test_that("downsampling", {
  data <- pbmc3k
  downsampled_data <- downsample(data, nreads=5000)
  a <- colSums(downsampled_data@assays$RNA@layers$counts)[1]
  expect_equal(a, 5000)
})
