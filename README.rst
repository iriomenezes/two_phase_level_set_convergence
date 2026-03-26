This code is used to produce the convergence test results in 

  "Numerical simulations and universal saturation profiles for viscous fingering patterns in Hele-Shaw flow" 
  by I.M. Coutinho, L.C. Morrow and S.W. McCue. 

The solver implements two‑phase Hele‑Shaw flow in radial geometry using a level‑set method to evolve the fluid–fluid interface. The formulation follows the standard approach: the pressure field is first computed, the interface is evolved by solving the level‑set equation, and reinitialisation is periodically performed to maintain the level‑set function as a signed‑distance function.

.. code-block:: bash

  N     : Number of grid points in each direction
  sigma : Non‑dimensional surface‑tension parameter
  eta1  : Viscosity of invading fluid
  eta2  : Viscosity of defending fluid
  tmax  : Final simulation time

The initial condition, lines 22–24, sets the initial radius of the bubble to be unity plus a 6‑fold sinusoidal perturbation of amplitude 0.1. Simulations are performed by running main.m, please note that the simulation time will increase with N.

For information on citing this work, please see https://zenodo.org/records/19210908

This work is licensed under a **Creative Commons Attribution 4.0 International License (CC BY 4.0)**.
You are free to share and adapt the material for any purpose, provided that appropriate credit is given.
