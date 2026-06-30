import matplotlib.pyplot as plt
import csv

X = []
Uinit = []
Ufinal = []

# Leer condición inicial
with open('initialCondition.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        if lines: # Evita errores si hay líneas en blanco al final
            X.append(float(lines[0]))
            Uinit.append(float(lines[1]))

# Leer estado final
with open('final.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        if lines:
            Ufinal.append(float(lines[1]))
        
# Generar la gráfica
plt.plot(X, Uinit, '.-', color='k', label="Init")
plt.plot(X, Ufinal, '.-', color='r', label="Final")
plt.legend()
plt.grid()
plt.title("Simulación de Ondas 1D - HPC MPI")
plt.xlabel("Posición (X)")
plt.ylabel("Amplitud (U)")
plt.show()