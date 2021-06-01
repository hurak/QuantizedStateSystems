# QuantizedStateSystems.jl

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hurak.github.io/QuantizedStateSystems.jl/stable)-->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hurak.github.io/QuantizedStateSystems.jl/dev)
[![Build Status](https://github.com/hurak/QuantizedStateSystems.jl/workflows/CI/badge.svg)](https://github.com/hurak/QuantizedStateSystems.jl/actions)
[![Coverage](https://codecov.io/gh/hurak/QuantizedStateSystems.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hurak/QuantizedStateSystems.jl)

A Julia package for solving initial value problems (IVP) for ordinary differential equations (ODE) based on discretization (quantization) of the state space (in contrast to the discretization/quantization of the time axis used in common numerical methods for solving ODEs). An approach (planned to be) implemented in this package is called [Quantized State System (QSS)](https://en.wikipedia.org/wiki/Quantized_state_systems_method) and has been developed and promoted in a series of papers by [Ernesto Kofman](https://scholar.google.com/citations?user=WdXDZEkAAAAJ&hl=en) and his colleagues, starting with the paper

- E. Kofman and S. Junco, “Quantized-state systems: a DEVS Approach for continuous system simulation,” Trans. Soc. Comput. Simul. Int., vol. 18, no. 3, pp. 123–132, Sep. 2001.


In addition, chapters dedicated to QSS are also in the following monographs

- F. E. Cellier and E. Kofman, Continuous System Simulation, Springer, 2010.

- B. P. Zeigler, A. Muzy, and E. Kofman, Theory of Modeling and Simulation: Discrete Event & Iterative System Computational Foundations, 3rd ed. Academic Press, 2018.


An implementation of the method has already appeared in [PowerDEVS](https://sourceforge.net/projects/powerdevs/)

- F. Bergero and E. Kofman, “PowerDEVS: a tool for hybrid system modeling and real-time simulation,” SIMULATION, vol. 87, no. 1–2, pp. 113–132, Jan. 2011, doi: 10.1177/0037549710368029.

The primary original motivation for the development of this package was just to get familiar with the algorithm(s) based on discretization/quantization of state variables, as they are not quite frequently described/discussed elsewhere. It remains to be seen if the package evolves to any reasonable maturity and if it can be of any other use for any other user.

## Related packages

Only after setting the name of the package I learnt that there is an identically named Julia package [QuantizedStateSystems.jl](https://github.com/BenLauwens/QuantizedStateSystems.jl). It seems to rest in the [WIP](https://en.wikipedia.org/wiki/Work_in_process) mode for quite some time. It is well possible that the presented package will follow the same trajectory, and therefore I am not going to resolve the unintended clash of the names of the two packages for the time being.  
