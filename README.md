# OpenSMOKE++ KineticOptimizer
Optimization of kinetic parameters of reactions included in a detailed kinetic mechanism

The OpenSMOKE++ Kinetic Optimizer is a tool for optimizing the kinetic parameters 
(frequency factor, temperature exponent, and activation energy) of reactions included in a detailed
kinetic mechanism to better fit experimental data and/or data coming from a more detailed/accurate
kinetic mechanism.

In the current version the optimization is carried out by minimizing the distance of ignition delay times
calculated in the following systems:
- batch reactors (adiabatic or with heat exchange)
- plug flow reactors (adiabatic or with heat exchange)
