import pandas as pd
import matplotlib.pyplot as plt

### Load KMC output
df = pd.read_csv('/Users/lotburgstra/Desktop/TCC_ASM/ASM_Q4/simulation_1/specnum_output.txt', delim_whitespace=True, comment='#')

### Compute coverages
total_sites = 40 * 25
df['Coverage_CO*'] = df['CO*'] / total_sites
df['Coverage_O*']  = df['O*']  / total_sites

### Plot CO* and O* coverage versus time
plt.figure(figsize=(8, 5))
plt.plot(df['Time'], df['Coverage_CO*'], label='CO* coverage')
plt.plot(df['Time'], df['Coverage_O*'],  label='O* coverage')
plt.xlabel('Time')
plt.ylabel('Coverage (fraction of sites)')
plt.title('Surface Coverages vs. Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

### Plot gas‐phase counts (O2, CO, CO2) versus time
plt.figure(figsize=(8, 5))
plt.plot(df['Time'], df['O2'],  label='O2')
plt.plot(df['Time'], df['CO'],  label='CO')
plt.plot(df['Time'], df['CO2'], label='CO2')
plt.xlabel('Time')
plt.ylabel('Number of Gas Molecules')
plt.title('Gas‐Phase Molecule Counts vs. Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
