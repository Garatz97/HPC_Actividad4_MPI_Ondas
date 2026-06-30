import matplotlib.pyplot as plt

# Datos extraídos de tus ejecuciones reales
procesos = [1, 2, 4]
tiempos = [28.08, 116.4, 793.1] # Medias calculadas de los logs
errores = [0.01987, 0.01405, 0.009934]

fig, ax1 = plt.subplots(figsize=(8, 5))

# Eje Y principal: Tiempo de ejecución
color_tiempo = 'tab:blue'
ax1.set_xlabel('Número de Procesos (Carga por proceso constante: 50k pts)')
ax1.set_ylabel('Tiempo de Cómputo (s)', color=color_tiempo)
ax1.plot(procesos, tiempos, marker='o', color=color_tiempo, linewidth=2, label='Tiempo (s)')
ax1.tick_params(axis='y', labelcolor=color_tiempo)
ax1.set_xticks(procesos)
ax1.grid(True)

# Eje Y secundario: Error L2
ax2 = ax1.twinx()  
color_error = 'tab:red'
ax2.set_ylabel('Error Numérico Global L2', color=color_error)
ax2.plot(procesos, errores, marker='s', linestyle='--', color=color_error, linewidth=2, label='Error L2')
ax2.tick_params(axis='y', labelcolor=color_error)

fig.suptitle('Análisis de Escalabilidad Débil (Weak Scaling)')
fig.tight_layout()
plt.show()