"""
    t,x = qss1(f,x₀,tspan,Δq)

Solve the initial value problem (IVP) for a given first-order ODE using the method of hysteretically quantized state system (QSS) or order 1.

Consider a single first-order explicit ODEs (aka state equation) `ẋ = f(x)`, with the initial value `x₀` specified at an initial time `t₀`. Solution is to be found on a time span `tspan=(t₀,t₁)`. The only parameter is the quantum `Δq`. The hysteresis band has the same width as the quantum.

The outputs are two arrays of arrays: the 1D array `t` of 1D arrays (vectors) of times and the correspoding 1D array `x` of the 1D arrays corresponding to the evolution of individual quantized state variables.

# References

1. E. Kofman and S. Junco, "Quantized-state systems: a DEVS Approach for continuous system simulation", Trans. Soc. Comput. Simul. Int., vol. 18, no. 3, pp. 123–132, Sep. 2001.
2. F. E. Cellier and E. Kofman, Continuous System Simulation, Springer, 2010.
3. B. P. Zeigler, A. Muzy, and E. Kofman, Theory of Modeling and Simulation: Discrete Event and Iterative System Computational Foundations, 3rd ed. Academic Press, 2018.

# Example

```julia
julia> f(x) = -1.0*x;
julia> x₀ = 10.0;
julia> tspan = (0.0,1.0);
julia> Δq = 0.1;
julia> tarray,xarray,qarray = qss1(f,x₀,tspan,Δq);
```
"""
function qss1(f,x₀,tspan,Δq)
    (t₀,t₁) = tspan             # Initial and final times.
    tarray = [t₀,]              # Initializing an array of times. The length not known apriori.
    xarray = [x₀,]              # Initializing an array of state variables.
    x = x₀                      # Initializing the state variable.
    t = t₀                      # Initializing the current time.
    q = Δq*floor(x/Δq)          # Quantized initial state.
    qarray = [q,]               # Initializing an array of quantized state variables.
    while t <= t₁
        k = f(q)                # Derivative corresponding to the quantized state.
        Δt = abs(Δq/k)          # Time increment.
        t += Δt                 # Next time.
        q += sign(k)*Δq         # Next quantized state.
        push!(tarray,t)         # Append the time to the array of times.
        push!(xarray,q)         # Append the quantized state to the array of (samples of) states.
                                # In the scalar case with no external events, the quantized states
                                # coincide with the true states at the transition times.
        push!(qarray,q)         # The same for the quantized states up to the first value.
                                # In presence of external events (functionality possibly added
                                # later) they could differ.
    end
    return (tarray,xarray,qarray)
end
