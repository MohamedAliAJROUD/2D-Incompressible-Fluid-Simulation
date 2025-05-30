# 2D-Incompressible-Fluid-Simulation
A MATLAB implementation of incompressible fluid flow simulation using Jacobi iterative methods for pressure projection and Runge-Kutta 4th order integration for particle advection.
# 2D Incompressible Flow Simulation with Lagrangian Particle Tracking

## Overview

This MATLAB project simulates two-dimensional incompressible fluid flow using a fractional step method and visualizes the dynamics via Lagrangian tracer particles. It demonstrates essential computational fluid dynamics (CFD) techniques, including velocity-pressure coupling, semi-Lagrangian advection, and Runge-Kutta integration.

Particles are injected through a horizontal jet, enabling the visualization of vortical structures and complex flow behavior.

## Features

- Incompressible flow solver using a velocity-pressure fractional step method
- Pressure projection via Jacobi iterative solver
- Particle advection using 4th-order Runge-Kutta integration
- Stable semi-Lagrangian advection for velocity transport
- Real-time visualization of the flow field and particle dynamics
- Optional video export to MP4 format for analysis and presentation

## Mathematical Methods

### Pressure Projection

To satisfy the incompressibility condition `div(v) = 0`, the following fractional step method is used:

1. Compute the provisional velocity field `v*`
2. Solve the Poisson equation:  
   `∇²p = -div(v*)`
3. Apply pressure correction:  
   `v = v* - grad(p)`

### Iterative Solver

The Poisson equation is solved using Jacobi iterations:  
`p(k+1) = J * p(k) + rhs / 2`  
where `J` is the Jacobi stencil operator and `rhs` is the divergence of the provisional velocity.

### Particle Advection (RK4)

Lagrangian particle trajectories are computed using 4th-order Runge-Kutta integration:

- Evaluate the velocity field at four intermediate positions
- Combine these to compute the final position update
- Temporal accuracy is maintained at order `O(Δt⁴)`

## Example Results
The simulation reproduces essential flow features including:

Jet injection: horizontal velocity source initiating momentum

Vortex generation: shear-induced rotational flow structures

Lagrangian dispersion: particle trajectories revealing streamlines

Flow evolution: time-dependent mixing and advection dynamics

## Applications
This simulation framework is relevant for educational and research purposes in several engineering domains:

HVAC and ventilation system optimization

Environmental modeling of pollutant transport

Microfluidics and lab-on-chip design

Wind farm wake flow analysis

Process engineering involving mixing and separation
