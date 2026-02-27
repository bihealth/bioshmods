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
  expect_error(infer_palette_variables(data.frame(x=1), cleanup=NA), "cleanup")
  expect_error(infer_palette_variables(data.frame(x=1), cleanup="yes"), "cleanup")
})

test_that("infer_palette_variables cleanup drops noisy categorical and constant variables", {
  df <- data.frame(
    id=as.character(seq_len(12)),                # categorical and all unique
    cat_many=c(paste0("L", 1:11), "L1"),        # categorical with >10 unique
    constant_num=rep(5, 12),                     # one unique value (any type)
    constant_chr=rep("A", 12),                   # one unique value (any type)
    stage=as.integer(rep(1:3, 4)),               # ordinal, should stay
    expr=seq_len(12) / 10,                       # continuous, should stay
    cohort=rep(c("A", "B"), 6),                  # categorical, should stay
    stringsAsFactors=FALSE
  )

  specs_clean <- infer_palette_variables(df, cleanup=TRUE)
  expect_named(specs_clean, c("stage", "expr", "cohort"))
  expect_equal(specs_clean$stage$type, "ordinal")
  expect_equal(specs_clean$expr$type, "continuous")
  expect_equal(specs_clean$cohort$type, "categorical")

  specs_raw <- infer_palette_variables(df, cleanup=FALSE)
  expect_true(all(c("id", "cat_many", "constant_num", "constant_chr") %in% names(specs_raw)))
})
