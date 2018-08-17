EXT (Exact)
===========

The EXT code solves the many-electron time-independent Schrödinger equation to calculate the fully correlated, ground-state wavefunction for a one-dimensional finite system of 2 or 3 spinless electrons interacting via the softened Coulomb repulsion :math:`(|x-x'|+1)^{-1}`. A perturbing potential is then applied to the ground-state system and its evolution is calculated exactly through solving the many-electron time-dependent Schrödinger equation.

Calculating the ground-state
----------------------------

An arbitrary wavefunction :math:`\Psi` is constructed (preferably close to the ground-state of the system) and then propagated through imaginary time using the Crank-Nicolson method. Using the eigenstates of the system :math:`\{\psi_{m}\}` as a basis:

.. math:: \Psi (x_{1}, x_{2}, \dots, x_{N}, \tau) = \sum\limits_{m} c_{m} e^{-E_{m}\tau}\psi_{m}.

Providing the wavefunction remains normalised, the limiting value is the ground-state of the system:

.. math:: E_{m+1} > E_{m} \ \ \forall \ m \in \mathbb{N}^{0} \implies \lim_{\tau \to \infty} \Psi (x_{1}, x_{2}, \dots, x_{N}, \tau) = \psi_{0}.

Time-dependence
---------------

A perturbing potential is applied to the Hamiltonian, :math:`\hat{H} = \hat{H}_{0} + \hat{\delta V}_{\mathrm{ext}}`. The system is initially in its ground-state and its evolution is calculated by propagating the ground-state wavefunction through real time using the Crank-Nicolson method.

One-electron systems
--------------------

EXT also works for systems of 1 electron. Unlike for 2 or 3 electron systems, the ground-state is calculated using an eigensolver. When a perturbation is applied to the system, its evolution is calculated by propagating the ground-state wavefunction through real time using the Crank-Nicolson method (like in the 2 or 3 electron systems).
