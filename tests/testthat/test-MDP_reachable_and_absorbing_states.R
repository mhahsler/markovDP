# absorbing states

data(Maze)
models_solve <- list(Maze)

s_abs <- c("s(1,4)", "s(2,4)")

verbose <- interactive()
#verbose <- TRUE

models_solve_no_chaching <- models_solve
for (m in models_solve_no_chaching) {
  # force recomputation
  m$absorbing_states <- NULL
  if (verbose)
    cat(m$name, "\n")
  tr <- absorbing_states(m, use_precomputed = FALSE, sparse = "states")
  if (verbose)
    print(tr)
  expect_equal(tr, s_abs)
}
  
# unreachable states
(s_unreach <- character(0))


for (m in models_solve_no_chaching) {
  # force recomputation
  m$unreachable_states <- NULL
  if (verbose)
    cat(m$name, "\n")
  tr <- unreachable_states(m, sparse = FALSE)
  if (verbose)
    print(tr)
  expect_equal(names(which(tr)), s_unreach)
}
  

# try with a simple two state model with an absorbing start data and and an
# unreachable state

m <- MDP(states = c("s1", "s2"), 
    actions = c("a1", "a2"), 
    transition_prob = rbind(P_(end.state = "s1", probability = 1)),
    reward = rbind(R_(value = -1)),
    start = "s1",
    name = "simple 2-state problem with absorbing start state"
    )

#str(m)

ms <- list(m = m,
        m_sparse = normalize_MDP(m, sparse = TRUE),
        m_dense = normalize_MDP(m, sparse = FALSE)
)

s_abs <- "s1"

for (m in ms) {
  if (verbose)
    cat(m$name, "\n")
  tr <- absorbing_states(m, use_precomputed = FALSE, sparse = "states")
  if (verbose)
    print(tr)
  expect_equal(tr, s_abs)
}


# unreachable states
s_unreach <- "s2"

for (m in ms) {
  if (verbose)
    cat(m$name, "\n")
  tr <- unreachable_states(m, sparse = FALSE)
  if (verbose)
    print(tr)
  expect_equal(names(which(tr)), s_unreach)
}

