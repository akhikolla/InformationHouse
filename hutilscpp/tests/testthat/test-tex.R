context("test-tex")

test_that("extractMandatory works", {
  library(hutils)
  x <- c('a', '', 'b', 'd', '{', 'e', '}', '.', 'b', 'e')
  res <- extractMandatory(x, c('b', 'd'), nCommands = 1L)[[1]]
  expect_equal(res,
               if_else(nzchar(res), x, ''))

  xop <- strsplit("qxy[ab]{jys}", split = "")[[1]]
  res_op <- extractMandatory(xop, c("q", "x", "y"), nCommands = 1L)[[1]]
  expect_true("y" %in% res_op)
})

test_that("Works with space before brace", {
  x <- c('a', '', 'b', 'd', ' ', '{', 'e', '}', '.')
  res <- extractMandatory(x, c('b', 'd'), nCommands = 1L)[[1]]
  expect_true(any(nzchar(x)))
})

test_that("Multiple optionals", {
  x <- strsplit("a \\Def[a [b] c]{df} x", split = "")[[1]]
  res <- extractMandatory(x, c("D", "e", "f"), 1L)
  expect_false("[" %in% res$support)
  expect_true("d" %in% res$support)

  x <- strsplit("a \\Def[a [b{q}] c]{df} x", split = "")[[1]]
  res <- extractMandatory(x, c("D", "e", "f"), 1L)
  expect_false("[" %in% res$support)
  expect_true("d" %in% res$support)
#
  x <- strsplit("a \\Def[a [b{q{}}] c]{df} xDe", split = "")[[1]]
  res <- extractMandatory(x, c("D", "e", "f"), 1L)
  expect_false("[" %in% res$support)
  expect_true("d" %in% res$support)

  x <- strsplit("a \\Defg[a [b{q{}}] c]{df} \\Def{a} b", split = "")[[1]]
  res <- extractMandatory(x, c("D", "e", "f"), 1L)
})

test_that("Bad document", {
  x <- strsplit("a{b", split = "")[[1]]
  res <- extractMandatory(x, c("foo"), 1L)
  expect_false(any(nzchar(res$support)))
  x <- strsplit("a[{b", split = "")[[1]]
  res <- extractMandatory(x, strsplit(c("foo"), split = "")[[1]], 1L)
  expect_false(any(nzchar(res$support)))
  x <- strsplit("b \\ad[s]", split = "")[[1]]
  res <- extractMandatory(x, c("a", "d"), 1L)
  expect_false(any(nzchar(res$support)))
})

test_that("Braces may appear within command", {
  x <- strsplit(c("the \\XYZ{cliff \\za{ba} wood.} flighty."), split = "")[[1L]]
  res <- extractMandatory(x, strsplit("XYZ", split = "")[[1]], 1L)
  expect_true(all(strsplit("cliff \\za{ba} wood.", split = "")[[1L]] %in% res$support))
})

test_that("Big popper", {
  skip_if_not_installed("TeXCheckR")
  skip_if_not_installed("data.table")
  skip_if_not(file.exists("~/AP-2018-retirement/report.tex"))
  library(TeXCheckR)
  library(data.table)
  Housing <- read_tex_document("~/AP-2018-retirement/report.tex")
  Housing_split = unlist(strsplit(Housing, split = ""))
  nFootnotes <- (length(grep("\\\\footnote(?![A-Za-z])[^\\{]*\\{", Housing, perl = TRUE)))
  footnote <- strsplit("footnote", split = "")[[1L]]

  res <- extractMandatory(Housing_split, footnote, nCommands = nFootnotes)

  Seq_All <- function(froms, tos) {
    unlist(lapply(seq_along(froms), function(i) {
      seq.int(from = froms[i], to = tos[i], by = 1L)
    }))
  }

  DT <-
    data.table(Ope = res$Openers,
               Clo = res$Closers)[, I := .I][, .(Text = paste0(res$support[.BY[[1]]:.BY[[2]]],
                                                               collapse = "")),
                                             keyby = c("Ope", "Clo", "I")]

  expect_true(any(grepl("For example, see: \\textcites[][]{IndustrySuperAustralia2015inquiryintoeconom",
                        DT$Text[1:10],
                        fixed = TRUE)))


})
