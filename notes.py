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
P_CO = 6.0e-9              ### bar
P_O2 = 1.2e-8              ### bar

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
def trans_2D(species_mass):
    return ((A * np.pi * species_mass * boltzmann_ev * temperature) / (h_ev ** 2))

### Call the different 2D translational values 
trans_CO = trans_2D(m_CO)
trans_O2 = trans_2D(m_O2)
trans_CO2 = trans_2D(m_CO2)

### Define a function for the calculation of the vibrational partition function values for the different species 
def vibrational_partitions(vib_species, temperature, h_ev, c, boltzmann_ev):
    vib_list = []
    for freq in vib_species:
        exponent = (-c * h_ev * freq) / (boltzmann_ev * temperature)
        q_vib = 1 / (1 - np.exp(exponent))
        vib_list.append(q_vib)
    return np.prod(vib_list)

### Call the different vibrational partition function values 
vib_O_surface = vibrational_partitions(O_surface, temperature, h_ev, c, boltzmann_ev)
vib_CO_surface = vibrational_partitions(CO_surface, temperature, h_ev, c, boltzmann_ev)
vib_O2_dissociation = vibrational_partitions(O2_dissociation, temperature, h_ev, c, boltzmann_ev)
vib_CO_oxidation = vibrational_partitions(CO_oxidation, temperature, h_ev, c, boltzmann_ev)

### Define a function to determine the CO adsorption reaction rate
def CO_ads():
    CO_ads_numerator = A * P_CO
    CO_ads_denominator = np.sqrt(2 * np.pi * m_CO * boltzmann_ev * temperature)
    return CO_ads_numerator / CO_ads_denominator

### Define a function to determine the CO desorption reaction rate
def CO_des():
    CO_des_numerator = (boltzmann_ev * temperature) * (vib_CO * rot_CO * trans_CO) * np.exp((CO_adsorption_energy / (boltzmann_ev * temperature)))
    CO_des_denominator = (h_ev * (vib_CO_surface ** 2))
    return CO_des_numerator / CO_des_denominator

### Define a function to determine the O2 adsorption reaction rate 
def O2_ads():
    O2_ads_numerator = (vib_O2_dissociation * 2 * A * P_O2 * np.exp(-O2_transition_energy / (boltzmann_ev * temperature)))
    O2_ads_denominator = 2 * vib_O2 * rot_O2 * trans_O2 * np.sqrt(2 * np.pi * m_O2 * boltzmann_ev * temperature)
    return O2_ads_numerator / O2_ads_denominator

### Define a function to determine the O2 desorption reaction rate 
def O2_des():
    O2_des_numerator = boltzmann_J * temperature * vib_O2_dissociation * np.exp((-O2_transition_energy - O2_adsorption_energy) / (boltzmann_ev * temperature))
    O2_des_denominator = h_ev * (vib_O_surface ** 2)
    return O2_des_numerator / O2_des_denominator

### Define a function to determine the CO oxidation reaction rate
def CO_oxi():
    CO_oxi_numerator = boltzmann_ev * temperature * vib_CO_oxidation * np.exp((-CO_transition_energy) / (boltzmann_ev * temperature))
    CO_oxi_denominator = h_ev * vib_CO_surface * vib_O_surface
    return CO_oxi_numerator / CO_oxi_denominator

### Define a function to determine the CO reversed oxidation reaction rate
def CO_oxi_rev():
    CO_oxi_rev_numerator = vib_CO_oxidation * 2 * A * P_CO * np.exp((-CO_transition_energy - CO_oxidation_energy) / (boltzmann_ev * temperature))
    CO_oxi_rev_denominator = 2 * vib_CO2 * rot_CO2 * trans_CO2 * np.sqrt(2 * np.pi * m_CO2 * boltzmann_ev * temperature)
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

