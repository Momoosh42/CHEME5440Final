import Plots
using Plots

Cell_doubling = 40 #minutes
Cell_volume = 1 #um^3
Cell_volume_L = Cell_volume*((1/1e6)^3)*((100)^3)*(1/1000)        # units: L
Cell_weight = 4.3e-13  #g
Water = 0.70      #proportion of cell that is water
char_enzyme_con = 0.3*(Cell_volume_L/((1-Water)*Cell_weight))*(1e9/1e6) #nmol/gDW
Prot_half_life = 24 #hr
Tl_initiation = 1.5/3600 #hr
Tl_sat_coefficient = 200*(Cell_volume_L/((1-Water)*Cell_weight))*(1e9/1e6)  #nmol/gDW
k_I_tl = 1/Tl_initiation #1/hr
Prot_deg = 1/(Prot_half_life/60)*log(2) #1/hr
dilution = log(2)/(Cell_doubling/60) #1/hr
Char_prot_length = 333 #aa
Prot_length = 300 #aa
ke_prot = 18 #aa/s @ doubling rate of 40min. BIND: 100059
Ribosomes = 26.3e3 #ribosomes/cell @ doubling rate of 40 min. BIND: 100059
Ribosome_concentration = Ribosomes*(1/6.02e23)*(1/((1-Water)*Cell_weight))*1e9 #nmol/gDW
avg_elongation_constant = ke_prot/Char_prot_length*3600 #1/h
Prot_elongation_constant = avg_elongation_constant*(Prot_length/Char_prot_length) #1/h

#assume that k_a is very small
tau_tl = Prot_elongation_constant/k_I_tl
V_max_prot = Prot_elongation_constant*Ribosome_concentration #nmol/(gDW*hr)

K_X =  0.575 #nmol/gDW according to prelim 1
r_L = V_max_prot*(K_X/(tau_tl*Tl_sat_coefficient)) #1/hr
K_L = (r_L/K_X)/(dilution+Prot_deg) #unitless

K_p = 1.5 #arbitrary value greater than 1

u_array = zeros(99)
p_star = zeros(99) #nmol/gDW
p_star_Kp = zeros(99)
p = zeros(99)
p_Kp = zeros(99)

j = 1

for u in 0.01:0.01:0.99
    p_star[j] = K_L*K_X*u #steady-state concentration, unbounded
    max_flux = V_max_prot*(p_star[j]/char_enzyme_con) #maximum protein flux nmol/(gDW*hr)

    p[j] = (max_flux/(dilution+Prot_deg))*u #bounded concentration nmol/gDW

#For polysome amplification being greater than 1
    p_star_Kp[j] = K_L*K_p*K_X*u
    max_flux_Kp = V_max_prot*(p_star_Kp[j]/char_enzyme_con)

    p_Kp[j] = (max_flux_Kp/(dilution+Prot_deg))*u
    u_array[j] = u
    global j = j+1
end


plot(u_array,p, label = "Kp = 1")
plot!(u_array, p_Kp, xlabel = "u", ylabel = "P concetration [nmol/gDW]", label = "Kp = 1.5")
