test_that("define_index lists and returns index definitions", {
  all_defs <- define_index(return = "data.frame")
  pianka <- define_index("pianka", return = "list")
  weis <- define_index("weis_sync", return = "data.frame")

  expect_true(is.data.frame(all_defs))
  expect_true("pianka_overlap" %in% all_defs$function_name)
  expect_equal(pianka$function_name, "pianka_overlap")
  expect_true("citation" %in% weis$field)
})
