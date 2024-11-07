'''
This python code is written by Kaur Shruti in an attempt to code the volume of change
which is depicted in Peter Mazur's Equation 19. 
'''

import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp
from . import _imaging as core
from . import _api, _version, cbook, docstring, rcsetup
from . import _ufuncs

'''
This part of the code focuses on deriving the permeability that is strongly 
depedent based off of the temperature. These are our parameters: 
k_0 = permeability constant that we know of from a given temperature constant 
T = Temperature 
T_0 = Known Temperature
b = temperature coefficient

The equation that I make use here is Equation (15):
k = k_{g} * e^(b * (T - T_{g}))
'''
def permeability(T, k_initial, b, T0):
    return k_initial * np.exp(b * (T - T0))

'''
This part is mainly coming from the equation (19)
'''
def dVdT(T, V, k_initial, A, v_initial, n2, Lf, R, b, T_initial):
    V = V[0] # First volume

    # Derivative of the initial value 
    dV_dT = np.gradient(V) 

    # Second derivative 
    d2V_dT2 = np.gradient(dV_dT) 

    # Left-hand side
    lhs = T * np.exp(b * (T_initial - T)) * d2V_dT2
    lhs -= ((b * T + 1) * np.exp(b * (T_initial - T)) - (A * R * k_initial * n2) / (B * (V + n2 * v_initial)) * (T**2 / V)) * dV_dT
    
    # Right-hand side
    rhs = (Lf * A * k_initial) / (B * v_initial)
    
    # Differential equation
    dV_dt = rhs - lhs
    
    return dV_dt

'''
These are the values that we have set based off of Peter Mazur's paper
called Kinetics of Water Loss and looking into the Table 1 named "Variables,
Parameters, and Constants in 
'''
L_f = 333.5 # Molar heat of fusion 
R = 8.314 # Gas content 
k_0 = 0.3 # Permeability constant 
b = 0.0325  # Temperature coefficient 
n_2 = 1.04e-10 # Solute molarity 
V_initial= 18-6 # Molar volume of the water 
A = 1e-5 # Surface area of our cell
B = 10 # Cooling rate (K/min)

# Hard coding the values for the initial and final temperature 
T_initial = 273 # Initital temperature of K
T_final = 240 # Final temperature of K
T_list_values = np.linspace(T_initial, T_final, 1000)

V0 = [V_initial]  # Starting volume

# Solve the ODE
solution = solve_ivp(dVdT, [T_list_values[0], T_list_values[-1]], V0, args=(k_0, A, V_initial, n_2, L_f, R, b, T_initial), t_eval=T_list_values)

# Extract the results
T_vals = solution.t
V_vals = solution.y[0]

# Output the result
for T, V in zip(T_vals, V_vals):
    print(f"Temperature: {T:.2f} K, Volume: {V:.4e} m^3")
    


