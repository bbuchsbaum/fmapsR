test_that("fm_package_info returns expected scaffold fields", {
  info <- fm_package_info()

  expect_type(info, "list")
  expect_identical(info$package, "fmapsR")
  expect_identical(info$status, "scaffold")
})
