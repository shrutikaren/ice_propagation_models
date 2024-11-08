
'''
This python code is written by Kaur Shruti in an attempt to code the volume of change
which is depicted in Peter Mazur's Equation 19. 
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def permeability(T, k_initial, b, T0):
    return k_initial * np.exp(b * (T - T0))

def dVdT(T, V, k_initial, A, v_initial, n2, Lf, R, b, T_initial):
    k_T = permeability(T, k_initial, b, T_initial)
    
    dV_dT = - (k_T * A * R * (T**2) * n2) / (B * V[0] * (V[0] + n2 * v_initial))
    
    lhs = T * np.exp(b * (T_initial - T)) * dV_dT
    rhs = (Lf * A * k_initial) / (B * v_initial)
    
    dV_dt = rhs - lhs
    return np.array([dV_dt])  # Ensure output is a 1D NumPy array

# Constants and parameters
L_f = 333.5       
R = 8.314         
k_0 = 0.3         
b = 0.0325        
n_2 = 1.04e-100   
V_initial = 88    
A = 1e-5          
B = 10            

T_initial = 272.15
T_final = 262.15
T_list_values = np.linspace(T_initial, T_final, 1000)

# Initial volume as a simple 1D array
V0 = [V_initial]

# Solve the ODE
solution = solve_ivp(
    dVdT,
    [T_list_values[0], T_list_values[-1]], V0,
    args=(k_0, A, V_initial, n_2, L_f, R, b, T_initial),
    t_eval=T_list_values,
    method='BDF',
    atol=1e-6, rtol=1e-6
)

# Extract and display results
T_vals = solution.t
V_vals = solution.y[0]

for T, V in zip(T_vals, V_vals):
    print(f"Temperature: {T:.2f} K, Volume: {V:.4e} m^3")

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(T_vals, V_vals, label="Volume vs. Temperature", color="blue")
plt.xlabel("Temperature (K)")
plt.ylabel("Volume of Water (m³)")
plt.title("Volume Change of Water in Cell with Temperature")
plt.legend()
plt.grid()
plt.gca().invert_xaxis()  # Optional: invert x-axis if you want decreasing T
plt.show()