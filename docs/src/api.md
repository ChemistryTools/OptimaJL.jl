```@meta
CurrentModule = Optima
```

# API Reference

## Problem definition

```@docs
OptimaProblem
OptimaOptions
OptimaState
OptimaResult
```

## Canonicalizer

```@docs
Canonicalizer
```

## Solver

`solve` and `solve!` are intentionally **not exported** to avoid naming conflicts
with `SciMLBase.solve`. Use the qualified names `Optima.solve(...)` /
`Optima.solve!(...)`, or add `import Optima: solve` at the top of your file.

```@docs
Optima.solve
Optima.solve!
```

## Sensitivity

```@docs
SensitivityResult
sensitivity
```

## SciML interface

```@docs
OptimaOptimizer
reset_cache!
```

## Internal components

The following symbols are exported for testing and extension purposes.
They are not needed for typical usage.

### KKT residual and Hessian

```@docs
KKTResidual
kkt_residual
hessian_diagonal
gibbs_hessian_diag
```

### Newton step

```@docs
NewtonStep
compute_step!
clamp_step
```

### Line search

```@docs
LineSearchFilter
line_search
```

### Variable stability

```@docs
classify_variables
reduced_step_for_unstable!
stability_measure
```
```
