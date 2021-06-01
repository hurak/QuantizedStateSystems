"""
    txarray,xarray,tqarray,qarray = qss1(f,x₀,tspan,Δq)

Solve the initial value problem (IVP) for a given first-order ODE using the method of hysteretically quantized state system (QSS) or order 1.

Consider a single first-order explicit ODEs (aka state equation) `ẋ = f(x)`, with the initial value `x₀` specified at an initial time `t₀`. Solution is to be found on a time span `tspan=(t₀,t₁)`. The only parameter for the method is the quantum `Δq`. The hysteresis band has the same width as the quantum.

The outputs are four arrays:
- a 1D array (aka vector) `t` of times when the value of `x` is given,
- a vector of values of `x` corresponding to the times,
- 1D array `x` of the 1D arrays corresponding to the evolution of the state variable,
- 1D array `q` of the 1D arrays corresponding to the evolution of the quantized state variable.

Note that for a first-order system with no external events the `x` and `q` arrays are identical up to the initial values (`x₀` before and after quantization), but we decided to keep track of both separately just in case we extend the functionality for handling external events (such as step inputs).

# Example

```julia
julia> f(x) = -1.0*x;
julia> x₀ = 10.0;
julia> tspan = (0.0,1.0);
julia> Δq = 0.1;
julia> txarray,xarray,tqarray,qarray = qss1(f,x₀,tspan,Δq);
```
"""
function qss1(f,x₀,tspan,Δq)
    (t₀,t₁) = tspan             # Initial and final times.
    tarray = [t₀,]              # Initializing an array of times. The length not known apriori.
    xarray = [x₀,]              # Initializing an array of state variables.
    t = t₀                      # Initializing the current time.
    x = x₀
    q = Δq*floor(x₀/Δq)         # Quantized initial state.
    qarray = [q,]               # Initializing an array of quantized state variables.
    while t <= t₁
        k = f(q)                # Derivative corresponding to the quantized state.
        Δt = abs(Δq/k)          # Time increment.
        t += Δt                 # Next time.
        q += sign(k)*Δq         # Next quantized state.
        x = q                   # Somewhat redundant in the no-input case but still...
        push!(tarray,t)         # Append the time to the array of times.
        push!(xarray,x)         # Append the quantized state to the array of (samples of) states.
                                # In the scalar case with no external events, the quantized states
                                # coincide with the true states at the transition times.
        push!(qarray,q)         # The same for the quantized states up to the first value.
                                # In presence of external events (functionality possibly added
                                # later) they could differ.
    end
    return (tarray,xarray,qarray)
end

"""
txarray,xarray,tqarray,qarray = qss1(f,x₀,tspan,Δq,tuarray,uarray)

Solve the initial value problem (IVP) for a given first-order ODE with inputs using the method of hysteretically quantized state system (QSS) or order 1.

Consider a single first-order explicit ODEs (aka state equation) `ẋ = f(x,u)`, with the initial value `x₀` specified at an initial time `t₀` and a piecewise (control) input `u` given by a pair of vectors `tuarray` and `uarray`of times and values. Solution is to be found on a time span `tspan=(t₀,t₁)`. The only parameter for the method is the quantum `Δq`. The hysteresis band has the same width as the quantum.

The outputs are four arrays:
- a 1D array (aka vector) `t` of times when the value of `x` is given,
- a vector of values of `x` corresponding to the times,
- 1D array `x` of the 1D arrays corresponding to the evolution of the state variable,
- 1D array `q` of the 1D arrays corresponding to the evolution of the quantized state variable.

Note that for a first-order system with no external events the `x` and `q` arrays are identical up to the initial values (`x₀` before and after quantization), but we decided to keep track of both separately just in case we extend the functionality for handling external events (such as step inputs).
"""
function qss1(f,x₀,tspan,Δq,tuarray,uarray)
    (t₀,t₁) = tspan             # Initial and final times for the simulation.
    if tuarray[1]==t₀           # Setting the initial value of the input:
        u = uarray[1]           # - either to the first value of the (piecewise) input,
        popfirst!(uarray)
        popfirst!(tuarray)
    else
        u = 0.0                 # - or to zero.
    end
    txarray = [t₀,]             # Initializing an array of times. The length not known apriori.
    xarray = [x₀,]              # Initializing an array of state variables.
    t = t₀                      # Initializing the current time.
    x = x₀                      # Initializing the current state.
    q = Δq*floor(x₀/Δq)         # Setting the quantized initial state.
    tqarray = [t,]              # Initializing an array of times for the quantized state variable.
    qarray = [q,]               # Initializing an array of (samples of) the quant. state variable.
    while t <= t₁
        k = f(q,u)              # Derivative corresponding to the quantized state.
        Δt = abs(Δq/k)          # Time increment based only on f and not u.
        if (~isempty(tuarray) && t+Δt <= tuarray[1]) || isempty(tuarray)   # If the next crossing of quant. level is before the change of u.
            t += Δt             # Set the next time.
            q += sign(k)*Δq     # Set the next quantized state.
            x = q               # State variable x is following q at the moment q changes.
            push!(txarray,t)    # Append the time to the array of times.
            push!(xarray,q)     # Append the quantized state to the array of (samples of) states.
            push!(tqarray,t)
            push!(qarray,q)     # The same for the quantized states (up to the first value).
        elseif ~isempty(tuarray)# If the change in u arrives before x crosses the quant. level.
            Δt = tuarray[1]-t   # The time from now till the moment u changes.
            t = tuarray[1]      # The moment when u changes.
            u = uarray[1]       # The new value of u.
            popfirst!(tuarray)  # Discard the time so that it is not used in the next iteration.
            popfirst!(uarray)   # Discard also the value of control for the same reason.
            x += k*Δt           # The state x gets updated but q does not.
            push!(txarray,t)    # Store the time when x experiences a change in derivative and
            push!(xarray,x)     # the actual value of x when that happens. No change to q and its t.
        end
    end
    return (txarray,xarray,tqarray,qarray)
end
