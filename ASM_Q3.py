import numpy as np 

### Masses of species (kg)
m_O2 = 5.3137e-26
m_CO = 4.6495e-26
m_CO2 = 7.3068e-26

### Useful constants 
boltzmann_ev = 8.617e-05    ### eV/K
boltzmann_J = 1.381e-23     ### J/K
temperature = 900           ### K
c = 2.998e+10               ### m/s
h_ev = 4.1357e-15           ### eV*s
h_J = 6.626e-34             ### J*s
pressure = 0.30             ### bar
h_bar_ev = 6.582e-16        ### eV*s
A = 1.5e-20                 ### m^2

### Pressure constants 
P_CO = 6.0e-4               ### bar
P_O2 = 1.2e-3               ### bar

### Vibrational partition function values 
vib_CO2 = 2.6717
vib_CO = 1.0336
vib_O2 = 1.0869

### Rotational partition function values 
rot_CO2 = 798.39
rot_CO = 325.44
rot_O2 = 218.39

### Vibrational partition function constants
O_surface = [405, 364, 364]
CO_surface = [1725, 334, 297, 297, 154, 154]
O2_dissociation = [994, 325, 194, 173, 114]
CO_oxidation = [1896, 544, 422, 370, 350, 290, 134, 95]

### Energy values
reaction_energy = -2.77
CO_adsorption_energy = -2.38
O2_transition_energy = 0.12
O2_adsorption_energy = -2.28
CO_transition_energy = -2.29
CO_oxidation_energy = 0.75

### Define a function to calculate the 2D translational partition function value 
def trans_2D(species_mass, boltzmann_J, temperature, h_J, A):
    return ((A * np.pi * species_mass * boltzmann_J * temperature) / (h_J ** 2))

### Call the different 2D translational values 
trans_CO = trans_2D(m_CO, boltzmann_J, temperature, h_J, A)
trans_O2 = trans_2D(m_O2, boltzmann_J, temperature, h_J, A)
trans_CO2 = trans_2D(m_CO2, boltzmann_J, temperature, h_J, A)

### Define a function for the calculation of the vibrational partition function values for the different species 
def vibrational_partitions(vib_species, temperature, h_J, c, boltzmann_J):
    vib_list = []
    for freq in vib_species:
        exponent = (-c * h_J * freq) / (boltzmann_J * temperature)
        q_vib = 1 / (1 - np.exp(exponent))
        vib_list.append(q_vib)
    return np.prod(vib_list)

### Call the different vibrational partition function values 
vib_O_surface = vibrational_partitions(O_surface, temperature, h_J, c, boltzmann_J)
vib_CO_surface = vibrational_partitions(CO_surface, temperature, h_J, c, boltzmann_J)
vib_O2_dissociation = vibrational_partitions(O2_dissociation, temperature, h_J, c, boltzmann_J)
vib_CO_oxidation = vibrational_partitions(CO_oxidation, temperature, h_J, c, boltzmann_J)

### Define a function to determine the CO adsorption reaction rate
def CO_ads(boltzmann_J, temperature, P_CO, A, m_CO):
    CO_ads_numerator = A * P_CO
    CO_ads_denominator = np.sqrt(2 * np.pi * m_CO * boltzmann_J * temperature)
    return CO_ads_numerator / CO_ads_denominator

### Define a function to determine the CO desorption reaction rate
def CO_des(boltzmann_J, temperature, h_J, vib_CO, rot_CO, boltzmann_ev, CO_adsorption_energy, trans_CO):
    CO_des_numerator = (boltzmann_J * temperature) * ((vib_CO * rot_CO * trans_CO) / vib_CO_surface) * np.exp((CO_adsorption_energy / (boltzmann_ev * temperature)))
    CO_des_denominator = (h_J * vib_CO_surface)
    return CO_des_numerator / CO_des_denominator

### Define a function to determine the O2 adsorption reaction rate 
def O2_ads(vib_O2_dissociation, vib_O2, rot_O2, trans_O2, boltzmann_ev, A, P_O2, m_O2, temperature, boltzmann_J):
    O2_ads_numerator = (vib_O2_dissociation * 2 * A * P_O2 * np.exp(-O2_transition_energy / (boltzmann_ev * temperature)))
    O2_ads_denominator = vib_O2 * rot_O2 * trans_O2 * np.sqrt(2 * np.pi * m_O2 * boltzmann_J * temperature)
    return O2_ads_numerator / O2_ads_denominator

###










