include("PhasePortrait.jl")

# Function for a dual repression system without cooperativity

# x1: range of x1 values (i.e. A values)

# x2: range of x2 values (i.e. R values)

# We use `@.` to apply the calculations across all rows.

# Note that model parameters are specified within the function

# Returns computed (dx1/dt, dx2/dt) over the range of (x1, x2)

function toggleMono(x1, x2)

    u = @. 1/(1+10*(((x2^2)/(0.1+x2^2))^2))-x1


    v = @. 1/(1+10*(((x1^2)/(0.1+x1^2))^2))-x2   #eqn 12

    return (u,v)

end

#Range of x1, x2 values

x1range = (0,1,40)          #Has the form (min, max, num points)

x2range = (0,1,40)          #Has the form (min, max, num points)

x₀ = ([0.0,0.0],)  #initial state vectors; a common must be included after the first array

tspan=(0.0,300.0)             #time span

#Call the phaseplot functon to construct the phase portrait

phaseplot(toggleMono, x1range, x2range, xinit=x₀, t=tspan, clines=true,

        norm=true, scale=0.5)
