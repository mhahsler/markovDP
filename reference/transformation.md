# Transformation Functions for Linear Function Approximation

Several popular transformation functions applied to state features used
in linear function approximation for
[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md).

## Usage

``` r
transformation_linear_basis(model, min = NULL, max = NULL, intercept = TRUE)

transformation_polynomial_basis(
  model,
  min = NULL,
  max = NULL,
  order,
  coefs = NULL
)

transformation_RBF_basis(model, min = NULL, max = NULL, centers, var = NULL)

transformation_fourier_basis(
  model,
  min = NULL,
  max = NULL,
  order,
  coefs = NULL
)

create_basis_coefs(dim, order)
```

## Arguments

- model:

  the [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) model.

- min, max:

  vectors with the minimum and maximum values for each feature. This is
  used to scale the feature to the \\\[0,1\]\\ interval for the Fourier
  basis.

- intercept:

  logical; add an intercept term to the linear basis?

- order:

  order for the Fourier basis.

- coefs:

  an optional matrix or data frame to specify the set of coefficient
  values for the Fourier basis (overrides `order`).

- centers:

  a scalar with the number of centers to create a a regular grid with
  that many steps per feature dimension. Alternatively, a matrix with
  the centers for the RBF can be supplied.

- var:

  a scalar with the variance used for the RBF.

- dim:

  number of features to describe a state.

## Value

A transformation function

## Details

The state feature function \\\phi()\\ uses the raw state feature vectors
\\\mathbf{x} = (x_1,x_2, ..., x_m)\\ which is either user-specified or
constructed by parsing the state labels of form `s(feature list)` and
then applies a transformation functions called basis functions.
Implemented basis functions are:

- Linear: no additional transformation is applied giving \\\phi_0(s) =
  1\\ for the intercept and \\\phi_i(s) = x_i\\ for \\i = \\1, 2, ...,
  m\\\\.

- Polynomial basis: \$\$\phi_i(s) = \prod\_{j=1}^m x_j^{c\_{i,j}},\$\$
  where \\c\_{1,j}\\ is an integer between 0 and \\n\\ for and order
  \\n\\ polynomial basis.

- Radial Basis: RBF.

- Fourier basis: \$\$\phi_i(s) = \text{cos}(\pi\mathbf{c}^i \cdot
  \mathbf{x}),\$\$ where \\\mathbf{c}^i = \[c_1, c_2, ..., c_m\]\\ with
  \\c_j = \[0, ..., n\]\\, where \\n\\ is the order of the basis. The
  components of the feature vector \\x\\ are assumed to be scaled to the
  interval \\\[0,1\]\\. The fourier basis transformation is implemented
  in `transformation_fourier_basis()`. `min` and `max` are the minimums
  and maximums for each feature vector component used to resale them to
  \\\[0,1\]\\ using \\\frac{x_i - min_i}{max_i - min_i}\\

  Details of this transformation are described in Konidaris et al
  (2011).

## References

Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning:
An Introduction. Second. The MIT Press.
[http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.md).

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

[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)

Other approximation:
[`linear_function_approximation`](http://michael.hahsler.net/markovDP/reference/linear_function_approximation.md),
[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
