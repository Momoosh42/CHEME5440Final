import Roots
import Plots
using Plots
using Roots
ke = 10^-4
ke_star = 5e-3
kf = 5.14e-21
kr = 2.5e-2
kdeg = 8e-4
V = 18
q = 10^3
n = 3e8


Kss = (ke_star*kf)/(ke*(kr+ke_star))

gam = 10^2
D = 10^-10


R_total = zeros(1000) #array of steady state values
z_val = zeros(1000)
j = 1

for z in 0.001:0.001:1
g(km) = (gam*z^2/D)^(1/3)- (km/(D*z))
km = fzero(g,1)

f(Rs_star) = Rs_star - Kss*V/ke_star*(n*kr*Rs_star+n*q)/(km + n*kf*(km/(n*kr*Rs_star+n*q))*((((((kr+ke)/kf)*Rs_star)^-1)-((n*kf)/((n*kr*Rs_star)+n*q)))^-1))
Rs_star = fzero(f, 1)
Ri_star = (ke_star/kdeg)*Rs_star
R_star = Rs_star+Ri_star
R_total[j] = R_star
z_val[j] = z
global j = j+1

end
norm_R = R_total./maximum(R_total)
plot(z_val,norm_R, xlabel = "z distance", ylabel = "predicted mitotic activity")
println(Kss)
