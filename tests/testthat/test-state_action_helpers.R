data(Maze)

# actions
res <- structure(1:2, levels = A(Maze), class = "factor")

expect_identical(normalize_action(1:2, Maze), res)
expect_identical(normalize_action(c("up", "right"), Maze), res)

expect_identical(normalize_action_id(1:2, Maze), as.integer(res))
expect_identical(normalize_action_id(c("up", "right"), Maze), as.integer(res))

expect_identical(normalize_action_label(1:2, Maze), as.character(res))
expect_identical(normalize_action_label(c("up", "right"), Maze), as.character(res))

# states
res <- structure(1:2, levels = S(Maze), class = "factor")
id <- 1:2
label <- c("s(1,1)", "s(2,1)")
factored <- rbind(c(1, 1), c(2, 1))

expect_identical(normalize_state(id, Maze), res)
expect_identical(normalize_state(label, Maze), res)
expect_identical(normalize_state(factored, Maze), res)

expect_identical(normalize_state_id(id, Maze), as.integer(res))
expect_identical(normalize_state_id(label, Maze), as.integer(res))
expect_identical(normalize_state_id(factored, Maze), as.integer(res))

expect_identical(normalize_state_label(id, Maze), as.character(res))
expect_identical(normalize_state_label(label, Maze), as.character(res))
expect_identical(normalize_state_label(factored, Maze), as.character(res))

# factored states

expect_identical(state2features("s(1,2)"), structure(
  c(1, 2),
  dim = c(1L, 2L),
  dimnames = list(c("s(1,2)"), c("x1", "x2"))))

expect_identical(state2features(c("s(1,1)", "s(2,1)")), structure(
  c(1, 2, 1, 1),
  dim = c(2L, 2L),
  dimnames = list(c("s(1,1)", "s(2,1)"), c("x1", "x2"))))

expect_error(state2features(c("s(1,1)", "s(2,1,3)")))
expect_error(state2features(c("s(1,)", "s(2,1)")))
expect_error(suppressWarnings(state2features(c("s(1,a)", "s(2,1)"))))

expect_identical(features2state(t(c(1, 1))), "s(1,1)")
expect_identical(features2state(rbind(c(1, 1), c(2, 1))), c("s(1,1)", "s(2,1)"))

