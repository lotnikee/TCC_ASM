import numpy as np 

### conversion factor eV to Joules
conversion = 1.602177e-19

### Masses of species (kg)
m_O2 = 5.3137e-26
m_CO = 4.6495e-26
m_CO2 = 7.3068e-26

### Useful constants 
boltzmann_joules = 1.380649e-23     ### boltzmann constant in J/K
temperature = 900                   ### temperature in K
c = 2.998e+08                       ### speed of light in m/s
h_joules = 6.626070e-34             ### Plank constant in J*s 
h_bar_joules = 1.054572e-34         ### Plank constant in J*s
standard_pressure = 101325          ### standard pressure in Pa
pressure = 0.30e+05                 ### total pressure in Pa
A = 1.5e-20                         ### area of the catalytic site in m^2

### Partial pressure constants (Pa)
P_CO = 6.0e-4              
P_O2 = 1.2e-3              

### Vibrational partition function values taken from Q2
vib_CO2 = 0.04639
vib_CO = 0.18639
vib_O2 = 0.3074

### Rotational partition function values taken from Q2
rot_CO2 = 798.39
rot_CO = 325.45
rot_O2 = 218.39

### Vibrational partition function constants (m^-1)
O_adsorbed = [40500, 36400, 36400]
CO_adsorbed = [172500, 33400, 29700, 29700, 15400, 15400]
O2_transition_dissociation = [99400, 32500, 19400, 17300, 11400]
CO_transition_oxidation = [189600, 54400, 42200, 37000, 35000, 29000, 13400, 9500]

### Energy values
reaction_energy = -2.77 * conversion
CO_adsorption_energy = -2.38 * conversion
O2_transition_adsorption_energy = 0.12 * conversion
O2_adsorption_energy = -1.14 * conversion
CO_transition_oxidation_energy = 1.23 * conversion
CO_oxidation_energy = 0.75 * conversion

### Define a function to calculate the 2D translational partition function value 
def trans_2D(species_mass):
    return (A * 2 * np.pi * species_mass * boltzmann_joules * temperature) / (h_joules ** 2)

### Call the different 2D translational values 
trans_CO = trans_2D(m_CO)
trans_O2 = 2 * trans_2D(m_O2)
trans_CO2 = 2 * trans_2D(m_CO2)

### Define a function for the calculation of the vibrational partition function values for the different species 
def vibrational_partitions(vib_species):
    vib_list = []
    for freq in vib_species:
        exponent = (-c * h_joules * freq) / (boltzmann_joules * temperature)
        q_vib = np.exp(-0.5 * exponent / (1 - np.exp(exponent)))
        vib_list.append(q_vib)
    return np.prod(vib_list)

### Call the different vibrational partition function values 
vib_O_adsorbed = vibrational_partitions(O_adsorbed)
vib_CO_adsorbed = vibrational_partitions(CO_adsorbed)
vib_O2_transition = vibrational_partitions(O2_transition_dissociation)
vib_CO_transition = vibrational_partitions(CO_transition_oxidation)

### Define a function to determine the CO adsorption reaction rate
def CO_ads():
    CO_ads_numerator = A * P_CO
    CO_ads_denominator = np.sqrt(2 * np.pi * m_CO * boltzmann_joules * temperature)
    return CO_ads_numerator / CO_ads_denominator

### Define a function to determine the CO desorption reaction rate
def CO_des():
    CO_des_numerator = (boltzmann_joules * temperature * vib_CO * rot_CO * trans_CO) * np.exp((CO_adsorption_energy / (boltzmann_joules * temperature)))
    CO_des_denominator = h_joules * vib_CO_adsorbed
    return CO_des_numerator / CO_des_denominator

### Define a function to determine the O2 adsorption reaction rate 
def O2_ads():
    O2_ads_numerator = (vib_O2_transition * 2 * A * P_O2) * np.exp(- O2_transition_adsorption_energy / (boltzmann_joules * temperature))
    O2_ads_denominator = (vib_O2 * rot_O2 * trans_O2) * np.sqrt(2 * np.pi * m_O2 * boltzmann_joules * temperature)
    return O2_ads_numerator / O2_ads_denominator

### Define a function to determine the O2 desorption reaction rate 
def O2_des():
    O2_des_numerator = (boltzmann_joules * temperature * vib_O2_transition) * np.exp((-O2_transition_adsorption_energy - O2_adsorption_energy) / (boltzmann_joules * temperature))
    O2_des_denominator = h_joules * (vib_O_adsorbed ** 2)
    return O2_des_numerator / O2_des_denominator

### Define a function to determine the CO oxidation reaction rate
def CO_oxi():
    CO_oxi_numerator = (boltzmann_joules * temperature * vib_CO_transition) * np.exp((-CO_transition_oxidation_energy) / (boltzmann_joules * temperature))
    CO_oxi_denominator = h_joules * vib_CO_adsorbed * vib_O_adsorbed
    return CO_oxi_numerator / CO_oxi_denominator

### Define a function to determine the CO reversed oxidation reaction rate
def CO_oxi_rev():
    CO_oxi_rev_numerator = (vib_CO_transition * 2 * A * P_CO) * np.exp((-CO_transition_oxidation_energy - CO_oxidation_energy) / (boltzmann_joules * temperature))
    CO_oxi_rev_denominator = (vib_CO2 * rot_CO2 * trans_CO2) * np.sqrt(2 * np.pi * m_CO2 * boltzmann_joules * temperature)
    return CO_oxi_rev_numerator / CO_oxi_rev_denominator

### Call all the reaction rates 
k_CO_ads = CO_ads()
k_CO_des = CO_des()
k_O2_ads = O2_ads()
k_O2_des = O2_des()
k_CO_oxi = CO_oxi()
k_CO_oxi_rev = CO_oxi_rev()

### Print all the reaction rate constants 
print("CO adsorption rate constant      = ", k_CO_ads)
print("CO desorption rate constant      = ", k_CO_des, "\n")
print("O2 adsorption rate constant      = ", k_O2_ads)
print("O2 desorption rate constant      = ", k_O2_des, "\n")
print("CO oxidation rate constant      = ", k_CO_oxi)
print("CO reverse oxidation rate constant      = ", k_CO_oxi_rev, "\n")

### Determine final expression
rate_constant = ((k_O2_ads / k_O2_des) ** 0.5) * (k_CO_ads / k_CO_des) * (k_CO_oxi / k_CO_oxi_rev)

### Print final expression 
print("Final reaction rate expresion   =", rate_constant)

