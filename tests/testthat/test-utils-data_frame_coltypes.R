test_that(".data_frame_coltypes classifies common column types", {
  df <- data.frame(
    int_col = 1:3,
    num_col = c(1.5, 2.5, 3.5),
    log_col = c(TRUE, FALSE, TRUE),
    chr_col = c("a", "b", "c"),
    fac_col = factor(c("x", "y", "x")),
    stringsAsFactors = FALSE
  )
  df$date_col <- as.Date("2026-01-01") + 0:2
  df$dt_col <- as.POSIXct("2026-01-01 12:00:00", tz = "UTC") + c(0, 60, 120)
  df$list_col <- I(list(1, 2, 3))

  expect_equal(
    .data_frame_coltypes(df),
    c(
      int_col = "integer",
      num_col = "numeric",
      log_col = "logical",
      chr_col = "character",
      fac_col = "factor",
      date_col = "date",
      dt_col = "datetime",
      list_col = "list"
    )
  )
})
