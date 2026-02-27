test_that("infer_palette_variables classifies columns and builds module-compatible specs", {
  df <- data.frame(
    sex=factor(c("F", "M", "F"), levels=c("F", "M", "U")),
    stage=as.integer(c(1, 2, 1)),
    expr=c(10.1, 11.4, 9.8),
    cohort=c("A", "B", "A"),
    stringsAsFactors=FALSE
  )

  specs <- infer_palette_variables(df)

  expect_named(specs, c("sex", "stage", "expr", "cohort"))

  expect_equal(specs$sex$type, "categorical")
  expect_equal(specs$sex$levels, c("F", "M", "U"))

  expect_equal(specs$stage$type, "ordinal")
  expect_equal(specs$stage$levels, c("1", "2"))

  expect_equal(specs$expr$type, "continuous")
  expect_equal(specs$expr$breaks, range(df$expr))

  expect_equal(specs$cohort$type, "categorical")
  expect_equal(specs$cohort$levels, c("A", "B"))
})

test_that("infer_palette_variables treats integer columns with many levels as continuous", {
  df <- data.frame(
    idx=as.integer(seq_len(12)),
    stringsAsFactors=FALSE
  )

  specs <- infer_palette_variables(df)
  expect_equal(specs$idx$type, "continuous")
  expect_equal(specs$idx$breaks, c(1, 12))
})

test_that("infer_palette_variables builds safe breaks for constant/empty numeric columns", {
  df_const <- data.frame(x=c(5, 5, 5), stringsAsFactors=FALSE)
  sp_const <- infer_palette_variables(df_const)
  expect_equal(sp_const$x$type, "continuous")
  expect_true(length(sp_const$x$breaks) == 2L)
  expect_true(diff(sp_const$x$breaks) > 0)

  df_na <- data.frame(x=c(NA_real_, NA_real_), stringsAsFactors=FALSE)
  sp_na <- infer_palette_variables(df_na)
  expect_equal(sp_na$x$type, "continuous")
  expect_equal(sp_na$x$breaks, c(0, 1))
})

test_that("infer_palette_variables validates inputs", {
  expect_error(infer_palette_variables(letters), "data frame")
  expect_error(infer_palette_variables(data.frame(x=1), ordinal_n_levels=0), "ordinal_n_levels")
})
