import numpy as np
import matplotlib.pyplot as plt

# Adquisición de datos y definición de las variables:

xcr = np.genfromtxt('Celdasx_malla.dat')
ycr = np.genfromtxt('Celdasy_malla.dat')
xfr = np.genfromtxt('Carasx_malla.dat')
yfr = np.genfromtxt('Carasy_malla.dat')
T = np.genfromtxt('Temperatura.dat')
info = np.genfromtxt('Residuo_norm.dat', delimiter = '\t')
it_V = info[:,0]
Rn_V = info[:,1]

# Posprocesamiento:

z = np.zeros((len(ycr),len(xcr)))
cont = 0

for j in range(len(ycr)):
  for i in range(len(xcr)):
    z[j,i] = T[cont]
    cont += 1

Dx = np.zeros((len(xcr)))
Dy = np.zeros((len(ycr)))

for i in range(len(Dx)):
  Dx[i] = xfr[i+1] - xfr[i]

for i in range(len(Dy)):
  Dy[i] = yfr[i+1] - yfr[i]

# Generación de las gráficas:

plt.figure(figsize = (10, 8))

for i in range(len(xfr)):
  plt.axvline(xfr[i], linewidth = '0.5')

for i in range(len(yfr)):
  plt.axhline(yfr[i], linewidth = '0.5')

plt.axvline(0,c = "k")
plt.axvline(xfr[len(xfr) - 1],c = "k")
plt.axhline(0,c = "k")
plt.axhline(yfr[len(yfr) - 1],c = "k")
xcr_aux = np.zeros((len(xcr) - 2))
ycr_aux = np.zeros((len(ycr) - 2))

for i in range(len(xcr_aux)):
  xcr_aux[i] = xcr[i+1]

for i in range(len(ycr_aux)):
  ycr_aux[i] = ycr[i+1]

plt.scatter(xcr_aux, ycr[len(ycr) - 1] * np.ones((1,len(xcr_aux))), s = 5, c = "r", label = "Celdas N (Tipo 1)")
plt.scatter(xcr_aux, ycr[0] * np.ones((1,len(xcr_aux))), s = 5, c = "b", label = "Celdas S (Tipo 2)")
plt.scatter(xcr[0] * np.ones((1,len(ycr_aux))), ycr_aux, s = 5, c = "g", label = "Celdas W (Tipo 3)")
plt.scatter(xcr[len(xcr) - 1] * np.ones((1,len(ycr_aux))), ycr_aux, s = 5, c = "m", label = "Celdas E (Tipo 4)")
plt.scatter(xcr[len(xcr) - 1], ycr[len(ycr) - 1], s = 5, c = "yellow", label = "Celda NE (Tipo 5)")
plt.scatter(xcr[0], ycr[len(ycr) - 1], s = 5, c = "violet", label = "Celda NW (Tipo 6)")
plt.scatter(xcr[0], ycr[0], s = 5, c = "orange", label = "Celda SW (Tipo 7)")
plt.scatter(xcr[len(xcr) - 1], ycr[0], s = 5, c = "c", label = "Celda SE (Tipo 8)")

for i in range(1, len(xcr) - 1):
  if (i == 1):
    plt.scatter(xcr[i] * np.ones(len(ycr_aux)), ycr_aux, s = 5, c = "k", label = "Celdas interiores (Tipo 9)")
  else:
    plt.scatter(xcr[i] * np.ones(len(ycr_aux)), ycr_aux, s = 5, c = "k")

plt.title(f"Malla con x = {len(xcr)} celdas, y = {len(ycr)} celdas")
plt.legend(bbox_to_anchor = (1.04, 1), borderaxespad = 0)
plt.xlabel("Longitud [m]")
plt.ylabel("Longitud [m]")
plt.savefig("Malla.pdf", format = 'pdf', bbox_inches='tight')

print(f"\nNúmero de caras en x con refinamiento = {len(xfr)} y tamaño de celda promedio = {np.mean(Dx):.4f}\n")
print(f"Vector de caras en x:\n{xfr}\n")
print(f"Número de caras en y con refinamiento = {len(yfr)} y tamaño de celda promedio = {np.mean(Dy):.4f}\n")
print(f"Vector de caras en y:\n{yfr}\n")

plt.figure(figsize = (12, 8))
plt.contourf(xcr, ycr, z, 30, cmap = 'inferno')
plt.colorbar(label = 'Temperatura [°C]', ticks = np.linspace(np.min(T), np.max(T), 15))
CS = plt.contour(xcr, ycr, z, 30, colors = 'k', linewidths = 0.5, linestyles = '--')
plt.clabel(CS, inline = 1, fontsize = 9)
plt.title(f"Distribución de temperatura en la placa en °C", y = 1.01)
plt.xlabel('Longitud [m]')
plt.ylabel('Longitud [m]')
plt.savefig("Dist_T.pdf", format = 'pdf')

plt.figure(figsize = (15,6))
plt.plot(it_V, Rn_V, "-r")
plt.title('Residual normalizado Vs. Iteraciones')
plt.xlabel('Iteraciones [-]')
plt.ylabel('Residual normalizado [-]')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Residuales.pdf", format = 'pdf')

print(f"Última iteración: {it_V[len(it_V) - 1]:n} con residual normalizado: {Rn_V[len(Rn_V) - 1]:.10f}\n")
