# Minimization for ill-conditioned problems
## Regularized quasi-Newton optimisation
Currently the only function, `rnewt` implements general-purpose regularized quasi-Newton optimisation routines as presented in [Kanzow and Steck (2023)](https://link.springer.com/article/10.1007/s12532-023-00238-4). The C++ code is written from scratch, and the More-Thuente linesearch script is an R-port specifically written for this implementation, but translated from [the python implementation associated to the article](https://github.com/dmsteck/paper-regularized-qn-benchmark/blob/d6777fa872bebcc38ebe2d7aa9dc21862d3b7ffd/utility/morethuente.py#L4).

## References
[Kanzow, C., & Steck, D. (2023). Regularization of limited memory quasi-Newton methods for large-scale nonconvex minimization. Mathematical Programming Computation, 15(3), 417-444.](https://link.springer.com/article/10.1007/s12532-023-00238-4)

[Sugimoto, S., & Yamashita, N. (2014). A regularized limited-memory BFGS method for unconstrained minimization problems. inf. téc.](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=abe09e57ca2985b7387f5875dfec307e22dacc4b)
