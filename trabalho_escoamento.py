#------------------------Fase Gás------------------------#
import numpy as np
import math

Tpc = 343.33
Ppc = 666.4
T = 194 #°F
P = 3810 #Psi
Tpr = T / Tpc
Ppr = P / Ppc
Mg = 30
R = 8.314

# Fator de Compressibilidade Isotérmico
a_z = 1.39 * ( (Tpr - 0.92) ** 1/2) - 0.36 * Tpr - 0.101
b_z = (0.62 - 0.23 * Tpr) * Ppr + ( (0.066 / (Tpr - 0.86) ) - 0.037) * Ppr ** 2 + (0.32 * Ppr ** 6 / (10** (9 * (Tpr - 1) ) ) )
c_z = 0.132 - 0.32 * np.log10(Tpr)
d_z = 10 ** (0.3106 - 0.49 * Tpr + 0.1824 * Tpr ** 2)

z = a_z + ( (1 - a_z) / (2.71 ** b_z)) + c_z * Ppr ** d_z   # euler ????????????????????????????????????


# Massa Específica do Gás
rho_g = P * Mg / (z * R * T)    # T em K Mg =????????????????????????????????????

# Viscosidade do Gás
x_ug = 3.47 + (1588 / T) + 0.0009 * Mg   # T em Rankine
y_ug = 1.66378 - 0.04679 * x_ug
a_ug = x_ug * rho_g ** y_ug
k_ug = ( 0.807 * ( Tpr ** 0.618) - 0.357 * np.exp( - 0.449 * Tpr ) + 0.34 * np.exp( - 4.058 * Tpr ) + 0.018) / (0.9490 * ( Tpc / ( (Mg ** 3) * (Ppc ** 4))) ** 1/6) # T em Rankine

ug = 10 ** - 4 * k_ug * np.exp (a_ug)

# Fator Volume Formação
bg = 14.7 / 60 * ( (z * T) / P)   # Psc = 14.7 psia e Tsc = 60 F

# Compressibilidade Isotérmica
dz_dp = 1 /(( ( 0.62 - 0.23 * Tpr ) * Ppr + ( 0.066 / (Tpr  - 0.86) - 0.37 ) * Ppr ** 2 + ( 0.32 / (10 ** (9 * Tpr - 9) )) * Ppr ** 6 ) * 2.71 ** (( ( 0.62 - 0.23 * Tpr ) + (( 0.132 / (Tpr - 0.86) ) - 0.074) * Ppr ) + ( 1.92 / (10 ** (9 * Tpr - 9) ) ) * Ppr ** 5 )) + ( 10 ** ( 0.3106 - 0.49 * Tpr + 0.1824 * Tpr ** 2) * ( Ppr ** (((10 ** (0.3106 - 0.49 * Tpr + 0.1824 * Tpr ** 2) ) - 1 ) * ( 0.132 - 0.32 * np.log10 (Tpr) ) ) ) )
cpr = 1 / Ppr - 1 / z * dz_dp

cg = cpr / Ppc

print(z)