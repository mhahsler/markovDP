# Linear Function Approximation

Approximate a Q-function a value function or a policy using linear
function approximation.

## Usage

``` r
q_approx_linear(model, transformation = transformation_linear_basis, ...)

v_approx_linear(model, transformation = transformation_linear_basis, ...)

pi_approx_linear(model, transformation = transformation_linear_basis, ...)

approx_value(f, state, action = NULL, w = NULL, model = NULL)
```

## Arguments

- model:

  a [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) model
  with defined state features.

- transformation:

  a transformation function. See
  [transformation](http://michael.hahsler.net/markovDP/reference/transformation.md).

- ...:

  further parameters are passed on to the
  [transformation](http://michael.hahsler.net/markovDP/reference/transformation.md)
  function.

- f:

  a linear approximation object.

- state:

  a state or state features.

- action:

  an action. If `NULL`, then the value for all available actions will be
  calculated.

- w:

  a weight vector to be used instead of `f$w`.

## Details

### Linear Approximation

#### Approximate Q Values

The state-action value function is approximated by \$\$\hat{q}(s,a) =
\boldsymbol{w}^\top\phi(s,a),\$\$

where \\\boldsymbol{w} \in \mathbb{R}^n\\ is a weight vector and \\\phi:
S \times A \rightarrow \mathbb{R}^n\\ is a feature function that maps
each state-action pair to a feature vector. Linear approximation has a
single optimum and can be optimized using a simple update rule following
the gradient of the state-action function \$\$\nabla
\hat{q}(s,a,\boldsymbol{w}) = \phi(s,a).\$\$

#### Approximate Value Function

Value function approximation works with state features \\\phi(s)\\:
\$\$\hat{v}(s) = \boldsymbol{w}^\top\phi(s)\$\$

The gradient is \$\$\nabla \hat{v}(s,\boldsymbol{w}) = \phi(s).\$\$

#### Approximate Policy

The most popular method implemented here uses a linear preference
function \$\$h(s,a,\boldsymbol{w}) = \boldsymbol{w}^\top\phi(s,a).\$\$

The action with the highest preference is the greedy action. To
represents a stochastic policy, the preference scores can be converted
into probabilities using the softmax function
\$\$\hat{\pi}(a\|s,\boldsymbol{w}) =
\frac{e^{h(s,a,\boldsymbol{w})}}{\sum_b e^{h(s,b,\boldsymbol{w})}}\$\$

The gradient is \$\$\nabla \hat{\pi}(a\|s,\boldsymbol{w}) =
\hat{\pi}(a\|s,\boldsymbol{w}) \[\phi(s,a) - \sum_b
\hat{\pi}(b\|s,\boldsymbol{w}) \phi(s,b)\].\$\$

The last term represents the feature vector reduced by the expected
feature vector across all actions under the current policy. This pushes
the approximation to make the chosen action \\a\\ more likely.

We use here the gradient of the log-policy which avoids multiplying by
\\\hat{\pi}(a\|s,\boldsymbol{w})\\ and is numerically more stable:
\$\$\nabla log\\ \hat{\pi}(a\|s,\boldsymbol{w}) = \phi(s,a)-\sum_b
\hat{\pi}(b\|s,\boldsymbol{w}) \phi(s,b).\$\$

The gradient of the log-policy points in the same direction as the
original gradient.

### State-action Feature Vector Construction

For a small number of actions, we can follow the construction described
by Geramifard et al (2013) which uses a state feature function \\\phi: S
\rightarrow \mathbb{R}^{m}\\ to construct the complete state-action
feature vector. Here, we also add an intercept term. The state-action
feature vector has length \\1 + \|A\| \times m\\. It has the intercept
and then one component for each action. All these components are set to
zero and only the active action component is set to \\\phi(s)\\, where
\\s\\ is the current state. For example, for the state feature vector
\\\phi(s) = (3,4)\\ and action \\a=2\\ out of three possible actions \\A
= \\1, 2, 3\\\\, the complete state-action feature vector is \\\phi(s,a)
= (0,0,0,1,3,4,0,0,0)\\. Each action component has three entries and the
1 represent the intercept for the state feature vector. The zeros
represent the components for the two not chosen actions.

The construction of the state-action values is implemented in
`add_linear_approx_Q_function()`.

The state feature function \\\phi(s)\\ starts with raw state feature
vector \\\mathbf{x} = (x_1,x_2, ..., x_m)\\ that are either
user-specified or constructed by parsing the state labels of form
`s(feature list)`. Then an optional nonlinear transformation can be
performed (see
[transformation](http://michael.hahsler.net/markovDP/reference/transformation.md)).

### Internal Representation

All approximations are lists with the elements:

- `x(s, a)` ... function to construct from state features, action-state
  features.

- `f(s, a, w)` ... approx. function

- `gradient(s, a, w)` ... gradient of f at w

- `w` ... the weight vector (initially all 0s)

- `transformation` ... a transformation kernel function that is applied
  to state features in x.

### Fitting

The weight vector are fitted using gradient descent by the MDP solver
(e.g.,
[solve_MDP_APPROX](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)).

### Prediction

`approx_value()` calculates approximate value given the weights in the
model or a specified weight vector.

## References

Alborz Geramifard, Thomas J. Walsh, Stefanie Tellex, Girish Chowdhary,
Nicholas Roy, and Jonathan P. How. 2013. A Tutorial on Linear Function
Approximators for Dynamic Programming and Reinforcement Learning.
Foundations and Trends in Machine Learning 6(4), December 2013, pp.
375-451. [doi:10.1561/2200000042](https://doi.org/10.1561/2200000042)

Konidaris, G., Osentoski, S., & Thomas, P. 2011. Value Function
Approximation in Reinforcement Learning Using the Fourier Basis.
Proceedings of the AAAI Conference on Artificial Intelligence, 25(1),
380-385.
[doi:10.1609/aaai.v25i1.7903](https://doi.org/10.1609/aaai.v25i1.7903)

## See also

Other approximation:
[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md),
[`transformation`](http://michael.hahsler.net/markovDP/reference/transformation.md)

## Examples

``` r
data(Maze)

# Approx Q function
f_q <- q_approx_linear(Maze)
f_q
#> q_approx_linear, approx_linear
#> 
#> transformation:
#> function (x) 
#> {
#>     x <- (x - min)/(max - min)
#>     if (intercept) 
#>         x <- c(x0 = 1, x)
#>     x
#> }
#> <bytecode: 0x55fb6ff495b0>
#> <environment: 0x55fb69ff7b00>
#> 
#> weights:
#>    up.x0    up.x1    up.x2 right.x0 right.x1 right.x2  down.x0  down.x1 
#>        0        0        0        0        0        0        0        0 
#>  down.x2  left.x0  left.x1  left.x2 
#>        0        0        0        0 

approx_value(f_q, state = "s(3,1)", action = "up", model = Maze)
#> [1] 0

# update the weights using a learning rate of .1 and a delta of 100
# ([] is used to preserve the names)
f_q$w[] <- .1 * 100 * f_q$gradient(s = s(3,1), a = "up", f_q$w)  
f_q
#> q_approx_linear, approx_linear
#> 
#> transformation:
#> function (x) 
#> {
#>     x <- (x - min)/(max - min)
#>     if (intercept) 
#>         x <- c(x0 = 1, x)
#>     x
#> }
#> <bytecode: 0x55fb6ff495b0>
#> <environment: 0x55fb69ff7b00>
#> 
#> weights:
#>    up.x0    up.x1    up.x2 right.x0 right.x1 right.x2  down.x0  down.x1 
#>       10       10        0        0        0        0        0        0 
#>  down.x2  left.x0  left.x1  left.x2 
#>        0        0        0        0 

approx_value(f_q, state = "s(3,1)", action = "up", model = Maze)
#> [1] 20

# find an approximate Q function using a solver
sol <- solve_MDP_APPROX(Maze, horizon = 1000, n = 100)
sol$solution$q_approx_linear
#> q_approx_linear, approx_linear
#> 
#> transformation:
#> function (x) 
#> {
#>     x <- (x - min)/(max - min)
#>     if (intercept) 
#>         x <- c(x0 = 1, x)
#>     x
#> }
#> <bytecode: 0x55fb6ff495b0>
#> <environment: 0x55fb69d83698>
#> 
#> weights:
#>        up.x0        up.x1        up.x2     right.x0     right.x1     right.x2 
#> -0.274100881 -0.165001024  0.044786426 -0.078516323 -0.554452173  0.334560775 
#>      down.x0      down.x1      down.x2      left.x0      left.x1      left.x2 
#> -0.468657328 -0.254496145 -0.005858767 -0.463828478 -0.254573472 -0.039261815 
approx_value(sol$solution$q_approx_linear, 
             state = "s(3,1)", action = "up", model = Maze)
#> [1] -0.4391019

# Approx V function
f_v <- v_approx_linear(Maze)
f_v
#> v_approx_linear, approx_linear
#> 
#> transformation:
#> function (x) 
#> {
#>     x <- (x - min)/(max - min)
#>     if (intercept) 
#>         x <- c(x0 = 1, x)
#>     x
#> }
#> <bytecode: 0x55fb6ff495b0>
#> <environment: 0x55fb6d6d2ab0>
#> 
#> weights:
#> x0 x1 x2 
#>  0  0  0 
approx_value(f_v, state = "s(3,1)", model = Maze)
#> [1] 0

# Approx Policy
f_pi <- pi_approx_linear(Maze)
f_pi
#> pi_approx_linear, approx_linear
#> 
#> transformation:
#> function (x) 
#> {
#>     x <- (x - min)/(max - min)
#>     if (intercept) 
#>         x <- c(x0 = 1, x)
#>     x
#> }
#> <bytecode: 0x55fb6ff495b0>
#> <environment: 0x55fb6e945020>
#> 
#> weights:
#>    up.x0    up.x1    up.x2 right.x0 right.x1 right.x2  down.x0  down.x1 
#>        0        0        0        0        0        0        0        0 
#>  down.x2  left.x0  left.x1  left.x2 
#>        0        0        0        0 

approx_value(f_pi, state = "s(3,1)", model = Maze)
#>    up right  down  left 
#>  0.25  0.25  0.25  0.25 

# update the weights using a learning rate of 0.1 and the gradient of
# the log-policy.
f_pi$w <- .1 * f_pi$gradient(s = s(3,1), a = "up", f_pi$w)  
approx_value(f_pi, state = "s(3,1)", model = Maze)
#>        up     right      down      left 
#> 0.2893358 0.2368881 0.2368881 0.2368881 
```
