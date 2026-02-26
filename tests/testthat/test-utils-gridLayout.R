test_that("gridLayout builds requested row/column structure", {
  ui <- gridLayout("A", "B", "C", .ncol=2, .nrow=2, .colwidths=6)
  html <- as.character(ui)

  expect_s3_class(ui, "shiny.tag")
  expect_equal(length(gregexpr("class=\"row\"", html, fixed=TRUE)[[1]]), 2)
  expect_equal(length(gregexpr("col-sm-6", html, fixed=TRUE)[[1]]), 4)
})

test_that("gridLayout respects .byrow filling order", {
  ui <- gridLayout("AA", "BB", "CC", "DD", .ncol=2, .nrow=2, .byrow=FALSE, .colwidths=6)
  html <- as.character(ui)

  pos_aa <- regexpr("AA", html, fixed=TRUE)[1]
  pos_bb <- regexpr("BB", html, fixed=TRUE)[1]
  pos_cc <- regexpr("CC", html, fixed=TRUE)[1]
  pos_dd <- regexpr("DD", html, fixed=TRUE)[1]

  expect_true(pos_aa < pos_cc)
  expect_true(pos_cc < pos_bb)
  expect_true(pos_bb < pos_dd)
})

test_that("gridLayout validates capacity and column widths", {
  expect_error(
    gridLayout("A", "B", "C", .ncol=2, .nrow=1),
    "cannot exceed"
  )

  expect_error(
    gridLayout("A", .ncol=2, .nrow=1, .colwidths=c(4, 4, 4)),
    "length 1 or `.ncol`"
  )

  expect_error(
    gridLayout("A", .ncol=2, .nrow=1, .colwidths=c(7, 6)),
    "must not be larger than 12"
  )
})
