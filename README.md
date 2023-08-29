# NodeDynamicsDemo.jl

This module and Jupyter notebook aim to illustrate the differences between simple
explicit, implicit, and analytical solutions to a system of linear differential
equations depicting a chain of energy nodes connected via diffusion.


## Status

This module isn't actively developed further, as its purpose was simply
to demonstrate the differences between different numerical methods for
solving systems of linear differential equations.


## Key contents

1. `demo.ipynb` is the main Jupyter notebook containing the derivations and examples on the different solutions.
2. `src/NodeDynamicsDemo.jl` contains module and functions to solve the dynamics for those who are interested in the implementation details.


## Installation

This module isn't available in any Julia registry, so the only way to obtain
it is to clone it on your own machine.
Afterwards, you can install the dependencies of this module by navigating
into its root folder *(the one containing the `Project.toml`)*,
and running
```julia
julia> using Pkg
julia> Pkg.instantiate()
```
or by running the `instantiate` command from julia [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/)
directly.

Since the `demo.ipynb` is a Jupyter notebook, you might need to install some
[Jupyter](https://jupyter.org/) stuff to get it working.
Personally, I just run it through  the appropriate vscode extensions.


## Usage

The `demo.ipynb` contains the math and examples to mess around with.
Go nuts!


## License

This code is licensed under the [MIT License](https://opensource.org/license/mit/).
See `LICENSE` for more information.


## Acknowledgements

All me, on my free time no less. And just to prove a point.
Why do I do this to myself?
