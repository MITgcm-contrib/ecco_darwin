import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('FCO2.dat',
                     dtype=float,
                     delimiter='\t')
# save time (s)
time = data[:,0]
# remove time column from outputs
data = data[:,1:]
data[data == 0] = np.nan

plt.figure(figsize=(10,6))
plt.imshow(data, cmap='RdBu', aspect='auto')
# Set tick labels
#months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D','J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

plt.xticks(range(0,np.size(data,1),5), range(0,np.size(data,1),5), rotation=90)
plt.yticks(range(0,np.size(data,0),1488*2), months, rotation=90)
plt.xlim(0, 130)
plt.xlabel("Distance from the mouth (km)")
plt.ylabel("Month")
# Add colorbar
#plt.colorbar(label= '[NO3] (mmol/m^3)')
plt.colorbar(label= 'Salinity (PSU)')


plt.savefig('FCO2_2.png')
plt.show()
plt.close()