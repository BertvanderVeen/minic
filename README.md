# Regularized quasi-Newton optimisation
This package implements general-purpose regularized quasi-Newton optimisation routines as presented in [Kanzow and Steck (2023)](https://link.springer.com/article/10.1007/s12532-023-00238-4). The C++ code is written from scratch, and the More-Thuente linesearch script is an R-port specifically written for this implementation, but translated from [the python implementation associated to the article](https://github.com/dmsteck/paper-regularized-qn-benchmark/blob/d6777fa872bebcc38ebe2d7aa9dc21862d3b7ffd/utility/morethuente.py#L4).

## References
Kanzow, C., & Steck, D. (2023). Regularization of limited memory quasi-Newton methods for large-scale nonconvex minimization. Mathematical Programming Computation, 15(3), 417-444.