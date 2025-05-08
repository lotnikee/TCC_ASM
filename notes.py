import numpy as np 

### Vibrational frequencies constants (cm^-1)
vib_freq_O2 = 1580
vib_freq_CO = 2143
vib_freq_CO2 = [2349, 1388, 667, 667]

### Moments of inertia constants (eV s^2)
O2_inertia = 1.22e-27
CO_inertia = 9.09e-28
CO2_inertia = 4.46e-27

### Symmetry numbers 
O2_symmetry = 2
CO_symmetry = 1
CO2_symmetry = 2

### Masses of species (kg)
m_O2 = 5.3137e-26
m_CO = 4.6495e-26
m_CO2 = 7.3068e-26

### Other constants 
boltzmann_ev = 8.617e-05    ### eV/K
temperature = 900           ### K
c = 2.998e+10               ### m/s
h_ev = 4.1357e-15           ### eV*s
pressure = 0.30             ### bar
partial_CO = 6.0e-09        ### bar 
partial_O2 = 1.2e-08        ### bar
h_bar_ev = 6.582e-16        ### eV*s

### Energy exponent 
energy_exponent = np.exp(-(-2.77) / (boltzmann_ev * temperature))
import numpy as np 

### Vibrational frequencies constants (cm^-1)
vib_freq_O2 = 1580
vib_freq_CO = 2143
vib_freq_CO2 = [2349, 1388, 667, 667]

### Moments of inertia constants (eV s^2)
O2_inertia = 1.22e-27
CO_inertia = 9.09e-28
CO2_inertia = 4.46e-27

### Symmetry numbers 
O2_symmetry = 2
CO_symmetry = 1
CO2_symmetry = 2

### Masses of species (kg)
m_O2 = 5.3137e-26
m_CO = 4.6495e-26
m_CO2 = 7.3068e-26

### Other constants 
boltzmann_ev = 8.617e-05    ### eV/K
temperature = 900           ### K
c = 2.998e+10               ### m/s
h_ev = 4.1357e-15           ### eV*s
pressure = 0.30             ### bar
partial_CO = 6.0e-09        ### bar 
partial_O2 = 1.2e-08        ### bar
h_bar_ev = 6.582e-16        ### eV*s

### Energy exponent 
energy_exponent = np.exp(-(-2.77) / (boltzmann_ev * temperature))

### Calculation of the mass contribution 
mass_contribution = (m_CO2 ** 1.5) / ((m_CO ** 1.5) * (m_O2 ** 0.75))

### Calculate the constant contributions 
constants_contribution_pressure = 1 / np.sqrt((boltzmann_ev * temperature) / pressure)
constants_contribution_h = 1 / ((2 * np.pi * boltzmann_ev * temperature) / (h_ev ** 2)) ** (3/4)
constants_contribution = constants_contribution_pressure * constants_contribution_h

### Calculate the vibrational partition function values for the different species 
def vibrational_partitions(temperature, h_J, c, boltzmann_J):
    vib_CO2_list = []
    for freq in vib_freq_CO2:
        exponent = (-c * h_J * freq) / (boltzmann_J * temperature)
        q_vib_CO2 = 1 / (1 - np.exp(exponent))
        vib_CO2_list.append(q_vib_CO2)
    return np.prod(vib_CO2_list)

q_vib_CO2 = vibrational_partitions(temperature, h_ev, c, boltzmann_ev)
q_vib_CO = 1 / (1 - np.exp((-c * h_ev * vib_freq_CO) / (boltzmann_ev * temperature)))
q_vib_O2 = 1 / (1 - np.exp((-c * h_ev * vib_freq_O2) / (boltzmann_ev * temperature)))

### Calculate the rotational partition function values for the different species 
rotational_constant = (2 * boltzmann_ev * temperature) / (h_bar_ev ** 2)
q_rot_CO = rotational_constant * (CO_inertia / CO_symmetry)
q_rot_O2 = rotational_constant * (O2_inertia / O2_symmetry)
q_rot_CO2 = rotational_constant * (CO2_inertia / CO2_symmetry)

### Calculate the combined vibrational and rotational partition function 
rot_vib_contribution_numerator = q_vib_CO2 * q_rot_CO2 * energy_exponent
rot_vib_contribution_denominator = (q_vib_CO * q_rot_CO) * ((q_rot_O2 * q_vib_O2) ** (1/2))
rot_vib_contribution = rot_vib_contribution_numerator / rot_vib_contribution_denominator

### Calculation of the CO2 partial pressure 
partial_CO2 = constants_contribution * mass_contribution * rot_vib_contribution * partial_CO * np.sqrt(partial_O2)

### Print calculation values to check if they are correct
print(f'Energy exponent = {energy_exponent:.5g}')
print(f'Mass contribution = {mass_contribution:.5g}')
print(f'Constants contribution = {constants_contribution:.5g}'"\n")

print(f"Vibrational partition function of CO2: {q_vib_CO2:.5g}")
print(f'Vibrational partition function value CO = {q_vib_CO:.5g}')
print(f'Vibrational partition function value O2 = {q_vib_O2:.5g}'"\n")

print(f"Rotational partition function of CO2: {q_rot_CO2:.5g}")
print(f'Rotational partition function value CO = {q_rot_CO:.5g}')
print(f'Rotational partition function value O2 = {q_rot_O2:.5g}'"\n")

print(f'Rotational and vibrational contributions = {rot_vib_contribution:.5g}')
print(f'Partial pressure CO2 = {partial_CO2:.5g}')


### Calculation of the mass contribution 
mass_contribution = (m_CO2 ** 1.5) / ((m_CO ** 1.5) * (m_O2 ** 0.75))

### Calculate the constant contributions 
constants_numerator = (h_ev) ** (3/2)
constants_denominator = (2 * np.pi) ** (3/4) * ((boltzmann_ev * 900) ** (5/4))
constants_contribution = constants_numerator / constants_denominator

### Calculate the vibrational partition function values for the different species 
def vibrational_partitions(temperature, h_J, c, boltzmann_J):
    vib_CO2_list = []
    for freq in vib_freq_CO2:
        exponent = (-c * h_J * freq) / (boltzmann_J * temperature)
        q_vib_CO2 = 1 / (1 - np.exp(exponent))
        vib_CO2_list.append(q_vib_CO2)
    return np.prod(vib_CO2_list)

q_vib_CO2 = vibrational_partitions(temperature, h_ev, c, boltzmann_ev)
q_vib_CO = 1 / (boltzmann_ev * temperature) / (1 - np.exp((-c * h_ev * vib_freq_CO) / (boltzmann_ev * temperature)))
q_vib_O2 = 1 / (boltzmann_ev * temperature) / (1 - np.exp((-c * h_ev * vib_freq_O2) / (boltzmann_ev * temperature)))

### Calculate the rotational partition function values for the different species 
rotational_constant = (2 * boltzmann_ev * temperature) / (h_bar_ev ** 2)
q_rot_CO = rotational_constant * (CO_inertia / CO_symmetry)
q_rot_O2 = rotational_constant * (O2_inertia / O2_symmetry)
q_rot_CO2 = rotational_constant * (CO2_inertia / CO2_symmetry)

### Calculate the combined vibrational and rotational partition function 
rot_vib_contribution_numerator = q_vib_CO2 * q_rot_CO2 * energy_exponent
rot_vib_contribution_denominator = (q_vib_CO * q_rot_CO) * ((q_rot_O2 * q_vib_O2) ** (1/2))
rot_vib_contribution = rot_vib_contribution_numerator / rot_vib_contribution_denominator

### Calculation of the CO2 partial pressure 
partial_CO2 = constants_contribution * mass_contribution * rot_vib_contribution * partial_CO * np.sqrt(partial_O2)

### Print calculation values to check if they are correct
print(f'Energy exponent = {energy_exponent:.5g}')
print(f'Mass contribution = {mass_contribution:.5g}')
print(f'Constants contribution = {constants_contribution:.5g}'"\n")

print(f"Vibrational partition function of CO2: {q_vib_CO2:.5g}")
print(f'Vibrational partition function value CO = {q_vib_CO:.5g}')
print(f'Vibrational partition function value O2 = {q_vib_O2:.5g}'"\n")

print(f"Rotational partition function of CO2: {q_rot_CO2:.5g}")
print(f'Rotational partition function value CO = {q_rot_CO:.5g}')
print(f'Rotational partition function value O2 = {q_rot_O2:.5g}'"\n")

print(f'Rotational and vibrational contributions = {rot_vib_contribution:.5g}')
print(f'Partial pressure CO2 = {partial_CO2:.5g}')

