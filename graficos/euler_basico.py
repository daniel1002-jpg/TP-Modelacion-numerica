import numpy as np
import matplotlib.pyplot as plt

# Parámetros
m = 200  # masa [kg]
k = 25000  # constante elástica [N/m]
lambda_ = 0  # amortiguación [N·s/m]
c = 0.1  # excitación del terreno [m]
dt = 0.005  # paso de tiempo [s]
t_max = 5  # tiempo máximo de simulación [s]

# Número de pasos de tiempo
num_pasos = int(t_max / dt)

# Inicialización de vectores de posición (u) y velocidad (v)
u = np.zeros(num_pasos)  # posición de la carrocería
v = np.zeros(num_pasos)  # velocidad de la carrocería

# Condiciones iniciales
u[0] = 0  # posición inicial
v[0] = 0  # velocidad inicial

# Resolver usando un método básico de Euler
for n in range(num_pasos - 1):
    # Aceleración calculada a partir de las fuerzas elástica y de amortiguación
    a = (k / m) * (c - u[n]) + (lambda_ / m) * (0 - v[n])  # (0 - v[n]) es c'
    
    # Actualización de la velocidad y la posición
    v[n + 1] = v[n] + a * dt  # velocidad
    u[n + 1] = u[n] + v[n + 1] * dt  # posición

# Graficar el resultado
t = np.arange(0, t_max, dt)
plt.plot(t, u)
plt.title("Posición de la carrocería sin amortiguación")
plt.xlabel("Tiempo (s)")
plt.ylabel("Posición (m)")
plt.grid(True)
plt.savefig("Posición de la carrocería sin amortiguación")
