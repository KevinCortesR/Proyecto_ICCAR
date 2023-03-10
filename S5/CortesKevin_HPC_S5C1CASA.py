import numpy as np
import matplotlib.pylab as plt

# Adquisición de datos y definición de las variables

info = np.genfromtxt('CortesKevin_HPC_S5C1.dat', delimiter='\t')
Dist_Uni = info[1:,0]
Dist_Gauss = info[1:,1]

# Generación de los histogramas

fig, (Uni, Gauss) = plt.subplots(2,1)
Uni.hist(Dist_Uni,bins=31)
Uni.set_title('Distribución uniforme')
Uni.grid()
Gauss.hist(Dist_Gauss,bins=200)
Gauss.set_title('Distribución gaussiana')
Gauss.grid()
fig.tight_layout()
plt.savefig("histogramas.pdf")
