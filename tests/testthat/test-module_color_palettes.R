library(shiny)

test_that("colorPalettesUI provides dynamic rows container", {
  ui <- colorPalettesUI("pal")
  html <- as.character(ui)

  expect_s3_class(ui, "shiny.tag")
  expect_true(grepl("pal-rows", html, fixed = TRUE))
})

test_that("colorPalettesServer initializes palettes for all variable types", {
  variables <- list(
    sex=list(type="categorical", levels=c("female", "male")),
    stage=list(type="ordinal", levels=c("I", "II", "III")),
    expression=list(type="continuous", breaks=c(-2, -1, 0, 1, 2))
  )

  palettes <- reactiveVal(list())

  testServer(
    colorPalettesServer,
    args=list(
      variables=variables,
      palettes=palettes
    ),
    {
      session$flushReact()
    }
  )

  p <- isolate(palettes())

  expect_named(p, c("sex", "stage", "expression"))

  expect_equal(p$sex$type, "categorical")
  expect_equal(p$sex$levels, c("female", "male"))
  expect_type(p$sex$pal, "character")
  expect_equal(names(p$sex$pal), c("female", "male"))

  expect_equal(p$stage$type, "ordinal")
  expect_equal(p$stage$levels, c("I", "II", "III"))
  expect_type(p$stage$pal, "character")
  expect_equal(names(p$stage$pal), c("I", "II", "III"))

  expect_equal(p$expression$type, "continuous")
  expect_equal(p$expression$breaks, c(-2, -1, 0, 1, 2))
  expect_true(is.function(p$expression$pal))
})

test_that("colorPalettesServer updates selected palette choices", {
  variables <- list(
    sex=list(type="categorical", levels=c("female", "male"))
  )

  palettes <- reactiveVal(list())

  testServer(
    colorPalettesServer,
    args=list(
      variables=variables,
      palettes=palettes
    ),
    {
      session$setInputs(palette_var_1="viridis:turbo")
      session$flushReact()
    }
  )

  p <- isolate(palettes())
  expect_named(p, "sex")
  expect_equal(p$sex$type, "categorical")
  expect_equal(names(p$sex$pal), c("female", "male"))
  expect_match(p$sex$pal[1], "^#")
})

test_that("colorPalettesServer accepts reactive variables input", {
  variables_rv <- reactiveVal(list(
    sex=list(type="categorical", levels=c("female", "male"))
  ))
  palettes <- reactiveVal(list())

  testServer(
    colorPalettesServer,
    args=list(
      variables=variables_rv,
      palettes=palettes
    ),
    {
      session$flushReact()
      expect_named(isolate(palettes()), "sex")

      variables_rv(list(
        expression=list(type="continuous", breaks=c(-2, 0, 2))
      ))
      session$flushReact()
      expect_named(isolate(palettes()), "expression")
      expect_equal(isolate(palettes())$expression$type, "continuous")
    }
  )
})

test_that("colorPalettesServer supports reversing continuous palettes", {
  variables <- list(
    expression=list(type="continuous", breaks=c(-2, -1, 0, 1, 2))
  )

  palettes <- reactiveVal(list())

  testServer(
    colorPalettesServer,
    args=list(
      variables=variables,
      palettes=palettes
    ),
    {
      session$setInputs(reverse_var_1=TRUE)
      session$flushReact()
    }
  )

  p <- isolate(palettes())
  cols <- p$expression$pal(c(-2, 2))
  ref <- viridisLite::viridis(5, option="viridis")

  expect_equal(cols[1], ref[5])
  expect_equal(cols[2], ref[1])
})

test_that("colorPalettesServer supports cycling categorical palettes", {
  variables <- list(
    sex=list(type="categorical", levels=c("female", "male", "unknown"))
  )

  palettes <- reactiveVal(list())

  base_cols <- {
    testServer(
      colorPalettesServer,
      args=list(
        variables=variables,
        palettes=palettes
      ),
      {
        session$flushReact()
      }
    )
    isolate(palettes()$sex$pal)
  }

  cycled_cols <- {
    testServer(
      colorPalettesServer,
      args=list(
        variables=variables,
        palettes=palettes
      ),
      {
        session$flushReact()
        session$setInputs(cycle_plus_var_1=1)
        session$flushReact()
      }
    )
    isolate(palettes()$sex$pal)
  }

  expect_equal(names(base_cols), names(cycled_cols))
  expect_equal(unname(cycled_cols), unname(c(base_cols[2:3], base_cols[1])))
})

test_that("colorPalettesServer validates malformed variables specification", {
  palettes <- reactiveVal(list())

  expect_warning(
    testServer(
      colorPalettesServer,
      args=list(
        variables=list(
          sex=list(type="categorical")
        ),
        palettes=palettes
      ),
      {
        session$flushReact()
      }
    ),
    "levels"
  )

  expect_warning(
    testServer(
      colorPalettesServer,
      args=list(
        variables=list(
          expression=list(type="continuous", breaks=0)
        ),
        palettes=palettes
      ),
      {
        session$flushReact()
      }
    ),
    "breaks"
  )
})
