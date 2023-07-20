test_that("ipg_singlesample provides errors", {
  df1 <- data.frame(time = seq(1,10, 1),
                    value = runif(10))
  expect_error(ipg_singlesample(df1, "time", "y.wrong.name"), regexp = "Y is not found in the data")
  expect_error(ipg_singlesample(df1, "time.wrong.name", "value"), regexp = "Time is not found in the data")
  expect_error(ipg_singlesample(df1, "time", "value", 1.5), regexp = "epsilon is outside the recommended range \\(0-1\\)")
  expect_error(ipg_singlesample(df1, "time", "value", -1.5), regexp = "epsilon is outside the recommended range \\(0-1\\)")

  df2 <- df3 <- df1
  df2$value <- as.character(df2$value)
  expect_error(ipg_singlesample(df2, "time", "value"), regexp = "Y must be a numeric variable")

  df3$time <- as.character(df3$time)
  expect_error(ipg_singlesample(df3, "time", "value"), regexp = "Time must be a numeric variable")
})

test_that("ipg_singlesample returns expected values", {
  df.test <- data.frame(time = seq(1,20,1),
                        value = c(rep(1,5), seq(2,9,1), rep(10, (20-5-(9-2+1)))))
  test.result <- ipg_singlesample(df.test, "time", "value")

  expect_equal(1, round(test.result$estimates$`peak growth rate`))
  expect_equal((13-5)/2+5, round(test.result$estimates$`peak growth time`))
  expect_equal(5, round(test.result$estimates$`lag time`))
  expect_equal(10, round(test.result$estimates$`max y`))
  expect_equal(round(log(2)/1, 1), round(test.result$estimates$`doubling time`, 1))
})

