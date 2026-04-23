test_that("fm_package_info returns expected development-state fields", {
  info <- fm_package_info()

  expect_type(info, "list")
  expect_identical(info$package, "fmapsR")
  expect_identical(info$status, "active-development")
  expect_identical(info$milestones, c("M1", "M2", "M3"))
})
