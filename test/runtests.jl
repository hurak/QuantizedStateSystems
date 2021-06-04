using QuantizedStateSystems
using LinearAlgebra
using Test

@testset "QSS1 for scalar systems with no input" begin
    @test begin
        a = -1.0
        f = (x,u) -> a*x
        x₀ = 10.0
        tspan = (0.0,1.0)
        Δq = 0.01
        txarray,xarray,tqarray,qarray = qss1(f,x₀,tspan,Δq,[],[])
        xtrue(t) = x₀*exp(a*t)
        xtruearray = xtrue.(txarray)
        relerror = (xarray-xtruearray)./xtruearray
        norm(relerror)<0.01
    end
end

@testset "QSS1 for scalar systems with (control) inputs" begin
    @test begin
        f = (x,u) -> -1.0*x + u
        x₀ = 10.0
        tspan = (0.0,5.0)
        tuarray = [1.76]
        uarray = [10.0]
        Δq = 0.01
        heaviside(t,τ) = t>=τ ? 1.0 : 0.0
        xtrue(t) = x₀*exp(-t) + 10.0*(heaviside(t,1.76)*(1.0-exp(-t+1.76)))
        (txarray,xarray,tqarray,qarray) = qss1(f,x₀,tspan,Δq,tuarray,uarray)
        xtruearray = xtrue.(txarray)
        relerror = (xarray-xtruearray)./xtruearray
        norm(relerror)<0.05
    end
end
