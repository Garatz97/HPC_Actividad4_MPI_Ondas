import matplotlib.pyplot as plt

# Datos de Strong Scaling (200.000 puntos fijos)
procesos = [1, 2, 4, 8]
tiempos = [477.9, 576.4, 793.1, 1480.0]
# Speedup = Tiempo_1_proceso / Tiempo_n_procesos
speedup = [tiempos[0]/t for t in tiempos]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Gráfica 1: Tiempo vs Procesos
ax1.plot(procesos, tiempos, marker='o', linestyle='-', color='b', linewidth=2)
ax1.set_title('Evolución del Tiempo (Strong Scaling)')
ax1.set_xlabel('Número de Procesos')
ax1.set_ylabel('Tiempo de Cómputo (s)')
ax1.grid(True)

# Gráfica 2: Speedup vs Procesos
ax2.plot(procesos, speedup, marker='s', linestyle='--', color='r', linewidth=2)
ax2.axhline(1, color='gray', linestyle=':', label='Ideal (Speedup=1)')
ax2.set_title('Speedup (200k puntos fijos)')
ax2.set_xlabel('Número de Procesos')
ax2.set_ylabel('Speedup (T1 / Tn)')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig('analisis_strong_scaling.png', dpi=300)
print("Gráfica guardada como 'analisis_strong_scaling.png'")