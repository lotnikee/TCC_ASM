import numpy as np 

### eV to Joules conversion rate 
conversion = 1.602177e-19

### Vibrational frequencies constants (m^-1)
vib_freq_O2 = 158000
vib_freq_CO = 214300
vib_freq_CO2 = [234900, 138800, 66700, 66700]

### Moments of inertia constants (J s^2)
O2_inertia = 1.22e-27 * conversion
CO_inertia = 9.09e-28 * conversion
CO2_inertia = 4.46e-27 * conversion

### Symmetry numbers 
O2_symmetry = 2
CO_symmetry = 1
CO2_symmetry = 2

### Masses of species (kg)
m_O2 = 5.3137e-26
m_CO = 4.6495e-26
m_CO2 = 7.3068e-26

### Initial energies (J)
energy_CO2 = -2.77 * conversion

### Partial pressures (Pa)
p_CO = 6.0e-04
p_O2 = 1.2e-03

### Other constants 
boltzmann_joules = 1.380649e-23                     ### J/K
boltzmann_ev = 8.617333e-05                         ### eV/K
temperature = 900                                   ### K
c = 2.998e+08                                       ### m/s
h_ev = 4.135668e-15                                 ### eV*s
h_joules = 6.626070e-34                             ### J*s 
h_bar_joules = 6.582120e-16 * conversion            ### eV*s
standard_pressure = 1e+05                           ### Pa

### Energy exponent 
energy_exponent = np.exp(-energy_CO2 / (boltzmann_joules * temperature))

### Calculation of the mass contribution 
mass_contribution = (m_CO2 / ((m_O2 ** 0.5) * m_CO)) ** (3/2)

### Calculate the constant contributions 
constant_numerator = h_joules ** 3
constants_denominator = (2 * np.pi) ** (3/2) * (boltzmann_joules * temperature) ** (5/2)
constants = (constant_numerator / constants_denominator) ** (1/2)

### Calculate the vibrational partition function values for the different species 
def vibrational_partitions(temperature, h_joules, c, boltzmann_joules):
    vib_CO2_list = []
    for freq in vib_freq_CO2:
        exponent = (-c * h_joules * freq) / (boltzmann_joules * temperature)
        q_vib_CO2 = np.exp((-0.5 * c * h_joules * freq) / (boltzmann_joules * temperature)) / (1 - np.exp(exponent))
        vib_CO2_list.append(q_vib_CO2)
    return np.prod(vib_CO2_list)

q_vib_CO2 = vibrational_partitions(temperature, h_joules, c, boltzmann_joules)
q_vib_CO = np.exp((-0.5 * c * h_joules * vib_freq_CO) / (boltzmann_joules * temperature)) / (1 - np.exp((-c * h_joules * vib_freq_CO) / (boltzmann_joules * temperature)))
q_vib_O2 = np.exp((-0.5 * c * h_joules * vib_freq_O2) / (boltzmann_joules * temperature)) / (1 - np.exp((-c * h_joules * vib_freq_O2) / (boltzmann_joules * temperature)))

### Calculate the rotational partition function values for the different species 
rotational_constant = (2 * boltzmann_joules * temperature) / (h_bar_joules ** 2)
q_rot_CO = rotational_constant * (CO_inertia / CO_symmetry)
q_rot_O2 = rotational_constant * (O2_inertia / O2_symmetry)
q_rot_CO2 = rotational_constant * (CO2_inertia / CO2_symmetry)

### Calculate the combined vibrational and rotational partition function 
rot_vib_contribution_numerator = q_vib_CO2 * q_rot_CO2 * energy_exponent
rot_vib_contribution_denominator = (q_vib_CO * q_rot_CO) * ((q_rot_O2 * q_vib_O2) ** (1/2))
rot_vib_contribution = rot_vib_contribution_numerator / rot_vib_contribution_denominator

### Calculation of the CO2 partial pressure 
p_CO2 = constants * mass_contribution * rot_vib_contribution * p_CO * (p_O2 ** 0.5) / (standard_pressure ** 0.5)

### Print calculation values to check if they are correct
print(f'Energy exponent = {energy_exponent:.5g}')
print(f'Mass contribution = {mass_contribution:.5g}')
print(f'Constants contribution = {constants:.5g}'"\n")

print(f"Vibrational partition function of CO2: {q_vib_CO2:.5g}")
print(f'Vibrational partition function value CO = {q_vib_CO:.5g}')
print(f'Vibrational partition function value O2 = {q_vib_O2:.5g}'"\n")

print(f"Rotational partition function of CO2: {q_rot_CO2:.5g}")
print(f'Rotational partition function value CO = {q_rot_CO:.5g}')
print(f'Rotational partition function value O2 = {q_rot_O2:.5g}'"\n")

print(f'Rotational and vibrational contributions = {rot_vib_contribution:.5g}')
print(f'Partial pressure CO2 = {p_CO2:.5g}')

