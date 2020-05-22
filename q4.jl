import Plots
import Roots
using Roots
using Plots
F6P = 0.1*1000 #uM
ATP = 2.3*1000 #uM
PFK = 0.12 #uM
K_F6P = 0.11*1000 #uM
K_ATP = 0.42*1000 #uM
k_cat = 0.4 * 3600 #1/s

r1 = k_cat*PFK*(F6P/(K_F6P+F6P))*(ATP/(K_ATP+ATP))

AMP = [0.000, 0.055, 0.093, 0.181, 0.405, 0.990] #mM
Rate = [3.003, 6.30, 29.761, 52.002, 60.306, 68.653] #uM/h
Confidence = [0.59,1.20,5.7,10.2,11.8,13.3] #95% interval +/-

#The activity function for this model can be assumed to be u(x) = (W1+W2*f)/(1+W1+W2*f)
#Where f = x^n/(K^n + x^n), the binding function
#When AMP = 0 the binding function is 0 so the only value for u(I) should be leaky expression
#W1/(1+W1)
W1 = (r1/Rate[1] - 1)^-1

plot(AMP,Rate, ribbon=(Confidence,Confidence), label = "actual data", lw = 3)

#When AMP = 0.990 mM the system is near saturation therefore nearly all of the inducer can be considered bound
#Then u(x) = (W1+W2)/(1+W1+W2) = 0.99
f(W2) = ((W1+W2)/(1+W1+W2))-0.99
W2 = fzero(f, 1)

print("W1 is ")
println(W1)
print("W2 is ")
println(W2)

calc_rate = zeros(length(Rate))


#Values of K and n were guessed to find the best seeming fit for the data
K = 0.5 #mM
n = 3.25

for i in 1:1:length(Rate)
f = AMP[i]^n/(K^n + AMP[i]^n)
calc_rate[i] = r1*((W1+f*W2)/(1+W1+f*W2))


end
print("K is ")
println(K)
print("n is ")
println(n)
print(calc_rate)


plot!(AMP,calc_rate, xlabel = "AMP concentration [mM]", ylabel = "overall rate [uM/hr]", label = "Calc rate", lw = 3)
