# NodeDynamicsDemo.jl

This module and Jupyter notebook aim to illustrate the differences between
explicit, implicit, and analytical solutions to a system of linear differential
equations depicting a chain of energy nodes connected via diffusion.
Before we get started, it's important to understand the concept of an
*energy node* for the purposes of this demonstration:

> Note that there are a plethora of methods for discretising linear differential equations. However, they are not discussed as they lose in simplicity to the Euler methods, and in accuracy to the analytical solution.


## Energy node

For the purposes of this demonstration, an *energy node* is defined entirely
by 5 properties:

1. `energyContent`: The current energy content on the node.
2. `selfDischarge`: The rate of self-discharge from the node relative to its `energy_content`.
3. `diffCoeffPrev`: The rate of energy diffusion from this node to the previous node in the chain.
4. `diffCoeffNext`: The rate of energy diffusion from this node to the next node in the chain.
5. `influx`: An optional influx of energy into the node at given times.

The exact units for the above properties don't matter for the purposes of this demonstration.
The first node in the chain is connected to a set `boundary_first`, and the last
node is connected to a set `boundary_last`. The boundaries can also be left empty,
in which case the first and last nodes are only connected to one node each.

> Note that in this demonstration, we're modelling a 1D chain of interlinked nodes. There's nothing preventing from connecting nodes arbitrarily, the math just becomes a tad more complicated.


## Discretising and solving the dynamics

For solving the dynamics of the energy node chain, we start with the differential
equation describing the energy balance on each node:
$$
\frac{dE_n(t)}{dt} = - \sigma_{n} E_n(t) + \rho_{n,n-1} \big( E_{n-1}(t) - E_n(t) \big) + \rho_{n,n+1} \big( E_n(t) - E_{n+1}(t) \big) + P_n(t)
$$
- $E_n(t)$ is energy content of node `n` at time `t`.
- $\sigma_{n}$ is the `selfDischarge` rate from node `n`.
- $\rho_{n,n-1}$ is the `diffCoeffPrev` of node `n`.
- $\rho_{n,n+1}$ is the `diffCoeffNext` of node `n`.
- $P_n(t)$ is the optional energy `influx` on node `n`.

The above differential equation serves as the basis for all the different methods detailed in the following sections.

> Note that the self-discharge and diffusion coefficients are assumed to be time-invariant for simplicity. For the simpler methods they don't really need to be, but can require a lot of matrix inversions if they cannot be assumed constant.


## Euler method

The simplest method to discretise and solve the root differential equation is the
[Euler Method](https://en.wikipedia.org/wiki/Euler_method), sometimes called
*foward Euler method* or *explicit Euler method*.
The method works simply by assuming that over a short enough period $\Delta t$,
the different time-varying parameters can be assumed to stay approximately constant.
Thus, the differential equation can be discretised and the energy content on the
next time step can be solved simply via:
$$
\frac{E_{n,t+\Delta t} - E_{n,t}}{\Delta t} = - \sigma_{n} E_{n,t} + \rho_{n,n-1} \big( E_{n-1,t} - E_{n,t} \big) + \rho_{n,n+1} \big( E_{n+1,t} - E_{n,t} \big) + P_{n,t} \\
\Leftrightarrow E_{n,t+\Delta t} = E_{n,t} + \big( - \sigma_{n} E_{n,t} + \rho_{n,n-1} \big( E_{n-1,t} - E_{n,t} \big) + \rho_{n,n+1} \big( E_{n+1,t} - E_{n,t} \big) + P_{n,t} \big) \Delta t
$$
> Note that solving the energy content on the next time step isn't really a concern for MILP models, since both of the above equations can easily be cast as linear constraints.

However, the key assumption here is *a short enough period $\Delta t$*. With dynamics too fast
for the used time step, this method is prone to numerical instabilities and oscillations.


## Backward Euler method

The next simplest method is the so-called
[Backward Euler method](https://en.wikipedia.org/wiki/Backward_Euler_method) or
*implicit Euler method*. Similar to its explicit counterpart, the method
assumes that time-varying stuff stays approximately constant over the time step.
However, instead of projecting the energy content of the node on the next time step,
the implicit method calculates the current energy content of the node accounting
for current changes, using the previous energy content only as a starting point:
$$
\frac{E_{n,t} - E_{n,t-\Delta t}}{\Delta t} = - \sigma_{n} E_{n,t} + \rho_{n,n-1} \big( E_{n-1,t} - E_{n,t} \big) + \rho_{n,n+1} \big( E_{n+1,t} - E_{n,t} \big) + P_{n,t} \\
$$
> Note that solving the energy content on the current time step isn't a concern for MILP models, as the above equation can be cast as a constraint as is

For an isolated node, the current energy content can be solved trivially.
However, for a system of nodes such as our chain, the whole system needs to be
solved in unison, as the current energy contents of the nodes depend on each other.
The problem can be cast as a system of linear equations:
$$
\underbrace{
    \begin{bmatrix}
    1 + \sigma_0 + \rho_{0,1} & - \rho_{0,1} & 0 & 0 & \dots & 0 & 0 & 0 \\
    - \rho_{1,0} & 1 + \sigma_1 + \rho_{1,2} + \rho_{1,0} & - \rho_{1,2} & 0 & \dots & 0 & 0 & 0 \\
    \dots & \dots & \dots & \dots & \dots & \dots & \dots & \dots \\
    0 & 0 & 0 & 0 & \dots & - \rho_{N-1,N-2} & 1 + \sigma_{N-1} + \rho_{N-1,N-2} + \rho_{N-1,N} & - \rho_{N-1,N} \\
    0 & 0 & 0 & 0 & \dots & 0 & - \rho_{N,N-1} & 1 + \sigma_N + \rho_{N,N-1}
    \end{bmatrix}
}_{\underline{A}}
\begin{bmatrix}
E_{0,t} \\ E_{1,t} \\ E_{2,t} \\ \dots \\ E_{N-2,t} \\ E_{N-1,t} \\ E_{N,t}
\end{bmatrix}
=
\begin{bmatrix}
E_{0,t-\Delta t} + P_{0,t} + \rho_{0,-1} E_{-1,t} \\
E_{1,t-\Delta t} + P_{1,t} \\
\dots \\
E_{N-1,t-\Delta t} + P_{N-1,t} \\
E_{N,t-\Delta t} + P_{N,t} + \rho_{N,N+1} E_{N+1,t}
\end{bmatrix}
$$
> Note that $E_{-1,t}$ and $E_{N+1,t}$ represent the set boundary conditions, if any. If the boundaries are left free, the terms are simply omitted.

Assuming that we know or can set the initial state of our system,
the right-hand side on the initial time step $t-\Delta t$ is known.
Thus, we can solve the vector of current $t$ energy contents of the system
by inverting matrix $\underline{A}$ and multiplying both sides with its inverse from the left.


## Analytical solution

With simplifying assumptions similar to those required by the Euler methods,
the root linear differential equation system can actually be solved analytically.
Casting the linear differential equation system in matrix form:
$$
\begin{bmatrix}
\dot{E}_0(t) \\
\dot{E}_1(t) \\
\dots \\
\dot{E}_{N-1}(t) \\
\dot{E}_{N}(t) \\
\end{bmatrix}
=
\underbrace{
    \begin{bmatrix}
    -\sigma_{0} -\rho_{0,1} & \rho_{0,1} & 0 & \dots & 0 & 0 \\
    \rho_{1,0} & -\sigma_{1} -\rho_{1,0} - \rho_{1,2} & \rho_{1,2} & 0 & \dots & 0 \\
    \dots & \dots & \dots & \dots & \dots & \dots \\
    0 & \dots & 0 & \rho_{N-1,N-2} & -\sigma_{N-1} -\rho_{N-1,N-2} - \rho_{N-1,N} & \rho_{N-1,N} \\
    0 & 0 & \dots & 0 & \rho_{N,N-1} & -\sigma_{N} -\rho_{N,N-1}
    \end{bmatrix}
}_{B}
\begin{bmatrix}
E_0(t) \\
E_1(t) \\
\dots \\
E_{N-1}(t) \\
E_{N}(t) \\
\end{bmatrix}
+
\begin{bmatrix}
P_0(t) \\
P_1(t) \\
\dots \\
P_{N-1}(t) \\
P_{N}(t) \\
\end{bmatrix}
\\
\Leftrightarrow \underline{\dot{E}}(t) = \underline{B} \underline{E}(t) + \underline{P}(t)
$$
> Here, the dot indicates time derivative and the underlines indicate matrices and/or vectors.

It can be shown, that for an invertible matrix $\underline{B}$, all matrices
$\underline{B}$, $\underline{B}^{-1}$, and $e^{\underline{B}t}$ are mutually
commutative. Furthermore, the derivative of the matrix exponential function
can be shown to be $\frac{d}{dt}e^{\underline{B}t} = \underline{B}e^{\underline{B}t}$.
**As long as $\underline{B}$ is time-invariant**
*(or at the very least can be assumed to be locally time-invariant)*,
the linear differential equation
system much like one would can be solved similar to
a regular linear differential equation:
$$
\underline{\dot{E}}(t) = \underline{B} \underline{E}(t) + \underline{P}(t) \\
\Leftrightarrow \underline{\dot{E}}(t) - \underline{B} \underline{E}(t) = \underline{P}(t) \qquad \big|\, e^{-\underline{B}t} \times \\
\Leftrightarrow e^{-\underline{B}t} \underline{\dot{E}}(t) - \underline{B} e^{-\underline{B}t} \underline{E}(t) = e^{-\underline{B}t} \underline{P}(t) \\
\Leftrightarrow \frac{d}{dt}\big( e^{-\underline{B}t} \underline{E}(t) \big) = e^{-\underline{B}t} \underline{P}(t) \qquad \big|\, \int_{t}^{t+\Delta t}dt \\
\Leftrightarrow \big( e^{-\underline{B}t} \underline{E}(t) \big)\big|_{t}^{t+\Delta t} = \int_{t}^{t+\Delta t} e^{-\underline{B}t} \underline{P}(t) dt
$$
> In order to proceed, we'll need some information about the `influx` and its time dependency. E.g. assuming it is constant over the desired time step $\Delta t$.

$$
\Leftrightarrow \big( e^{-\underline{B}t} \underline{E}(t) \big)\big|_{t}^{t+\Delta t} = -\underline{B}^{-1} \big( e^{-\underline{B}t} \big)\big|_{t}^{t+\Delta t} \underline{P}_t \\
\Leftrightarrow e^{-\underline{B}(t + \Delta t)} \underline{E}(t + \Delta t) - e^{-\underline{B}t} \underline{E}(t) = -\underline{B}^{-1} \big( e^{-\underline{B}(t + \Delta t)} - e^{-\underline{B}(t)} \big) \underline{P}_t \qquad \big|\, e^{\underline{B}(t + \Delta t)} \times \\
\Leftrightarrow \underline{E}(t + \Delta t) = e^{\underline{B}\Delta t} \underline{E}(t) - (\underline{I} - e^{\underline{B}\Delta t})\underline{B}^{-1} \underline{P}_t
$$
allowing us to solve the energy content of the nodes at time $t + \Delta t$
based on the energy content at time $t$ and some fancy matrix math on the
coefficient matrix $\underline{B}$.

> Note that here we're talking about *time*, not *time steps*. The only "discretization" we need is some assumption on how $\underline{P}(t)$ behaves. Assuming it constant is not necessary, but greatly simplifies tha calculations.