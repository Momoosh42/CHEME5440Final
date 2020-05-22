README


PhasePortrai.jl: This program is similar to the one provided earlier in the course
    except that the initial condition plot has been disabled.

ToggleMono.jl: This program was used to create the plot in 1b. This program
    creates a phase portrait of the ODEs related to delta expression in cells.
    The plot shows three fixed points with two being stable and the central
    point being unstable. The two stable points are expressions of the two
    different cell fates, while the unstable point is when there is exactly
    an equal amount of delta in both cells. Any perturbation at this point will lead
    to a definitive cell fate determination.

2d.jl: This script produced the plot in 2d. Using the given values an expression
    of km(z) was found for any value of z. This value was used to solve for value
    of Rs_star. Rs_star and Ri_star are related so using this relationship R_star
    could be easily calculated. These values were then normalized and plotted against
    z, the distance from the beginning of the monolayer

3b.jl: This script produced the plot for 3b and 3c. Using the given values and some values
    from Bionumbers, along with provided for equations in both the exam and previous
    lectures, a steady state protein concentration was found. This concentration
    was then utilized to find the maximum flux of the system and therefore the
    predicted protein concentration. The polysome amplification constant was manipulated
    to show its influence on protein production.

q4.jl: This script produced the plot for 4c. Using assumptions of when the bound function
    were equal to zero or nearly equal to 1, the constants W1 and W2 were found
    using the Roots package. This formulation was inspired by Moon et al. The actual
    data was plotted and using visual comparison values of K and n were found that
    seemed to best-fit the predicted values to the actual values. A plot was then produced
    that showed the 95% confidence interval of the known data and where the predicted data fell
    within the confidence interval.
    W1 = 0.045
    W2 = 98.95
    K = 0.5 mM
    n = 3.25
