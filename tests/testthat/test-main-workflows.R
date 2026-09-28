two_state_mdp <- function() {
  MDP(
    states = c("start", "goal"),
    actions = c("go", "wait"),
    transition_model = list(
      go = rbind(c(0, 1), c(0, 1)),
      wait = "identity"
    ),
    reward = rbind(
      R_(value = 0),
      R_(action = "go", start.state = "start", end.state = "goal", value = 1)
    ),
    discount = 0.9,
    start = "start"
  )
}

test_that("an MDP specifies transitions, rewards, and a start state", {
  model <- two_state_mdp()

  expect_s3_class(model, "MDP")
  expect_identical(S(model), c("start", "goal"))
  expect_identical(A(model), c("go", "wait"))
  expect_equal(start_vector(model, sparse = FALSE), c(start = 1, goal = 0))
  expect_equal(transition_matrix(model, "go", "start", sparse = FALSE),
               c(start = 0, goal = 1))
  expect_equal(reward_matrix(model, "go", "start", "goal"), 1)
  expect_equal(reward_matrix(model, "wait", "start", "start"), 0)
  expect_identical(absorbing_states(model, sparse = "states"), "goal")

  step <- act(model, "start", "go")
  expect_equal(step$reward, 1)
  expect_identical(as.character(step$state_prime), "goal")
})

test_that("solving an MDP produces a usable policy and value function", {
  model <- two_state_mdp()
  expect_false(is_solved_MDP(model))
  expect_error(policy(model), "not solved")

  solved <- solve_MDP(model, progress = FALSE)
  expect_true(is_solved_MDP(solved))
  expect_true(is_converged_MDP(solved))
  expect_identical(policy(solved)$state, S(model))
  expect_identical(as.character(policy(solved)$action[1]), "go")
  expect_equal(value_function(solved), c(start = 1, goal = 0))
  expect_equal(expected_return(solved), 1)
  expect_equal(unname(Q_values(solved)["start", ]), c(1, 0.9))
  expect_identical(action(solved, "start", as = "label"), "go")
  expect_equal(act(solved, "start")$reward, 1)

  finite <- solve_MDP(model, horizon = 2, progress = FALSE)
  expect_length(policy(finite, drop = FALSE), 2)
  expect_equal(unname(value_function(finite, drop = FALSE)["start", ]), c(1, 1))
  expect_identical(as.character(policy(finite, epoch = 1)$action[1]), "go")
})

test_that("manual policies change actions and expected returns", {
  model <- two_state_mdp()
  waiting <- manual_policy(model, c("wait", "go"), estimate_V = TRUE)
  expect_identical(waiting$state, S(model))
  expect_equal(waiting$V, c(0, 0))

  with_policy <- add_policy(model, waiting)
  expect_identical(action(with_policy, "start", as = "label"), "wait")
  expect_equal(policy_evaluation(model, waiting, progress = FALSE),
               c(start = 0, goal = 0))
  expect_equal(expected_return(with_policy, method = "policy_evaluation", progress = FALSE), 0)
  expect_equal(unname(induced_transition_matrix(with_policy)), diag(2))
})

test_that("both MDP sampling engines follow a deterministic policy", {
  solved <- solve_MDP(two_state_mdp(), progress = FALSE)

  for (engine in c("r", "cpp")) {
    simulation <- sample_MDP(solved, n = 3, horizon = 2,
                             engine = engine, trajectories = TRUE, progress = FALSE)
    expect_equal(simulation$reward, rep(1, 3))
    expect_equal(simulation$avg_return, 1)
    expect_identical(as.character(simulation$trajectories$s), rep("start", 3))
    expect_identical(as.character(simulation$trajectories$s_prime), rep("goal", 3))
    expect_identical(as.character(simulation$trajectories$a), rep("go", 3))
  }
})

test_that("sample-access MDPs support actions and trajectories without stored states", {
  transition <- function(model, action, state) {
    if (all(state == s(1)))
      list(reward = 0, state_prime = state)
    else
      list(reward = 1, state_prime = s(1))
  }
  model <- MDPSample(actions = "go", transition_model = transition,
                     start = s(0), absorbing_states = s(1))

  expect_s3_class(model, "MDPSample")
  expect_null(S(model))
  step <- act(model, s(0), "go")
  expect_equal(step$reward, 1)
  expect_equal(step$state_prime, s(1))

  simulation <- sample_MDP(model, n = 3, horizon = 2,
                           trajectories = TRUE, progress = FALSE)
  expect_equal(simulation$reward, rep(1, 3))
  first_steps <- simulation$trajectories[simulation$trajectories$time == 0, ]
  expect_identical(first_steps$s, rep("s(0)", 3))
  expect_identical(first_steps$s_prime, rep("s(1)", 3))
})

test_that("gridworld helpers create the expected state layout and movement", {
  model <- gw_maze_MDP(c(2, 2), start = "s(1,1)", goal = "s(2,2)")

  expect_equal(gw_matrix(model, what = "states"),
               matrix(c("s(1,1)", "s(2,1)", "s(1,2)", "s(2,2)"), nrow = 2))
  expect_equal(unname(transition_matrix(model, "right", "s(1,1)", sparse = FALSE)),
               c(0, 0, 1, 0))
  expect_equal(unname(transition_matrix(model, "up", "s(1,1)", sparse = FALSE)),
               c(1, 0, 0, 0))
  expect_identical(absorbing_states(model, sparse = "states"), "s(2,2)")
})
