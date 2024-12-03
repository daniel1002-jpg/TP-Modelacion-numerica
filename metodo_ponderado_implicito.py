import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
# from random import seed

# Siempre seteamos la seed de aleatoridad para que los # resultados sean reproducibles
# seed(12345)
# np.random.seed(12345)

# sns.set_theme()

# Función de exitación externa c(t)
def c(t):
	if t >= 1 and t < 1.2:
		return float(t - 1.0)
	elif (t >= 1.2 and t < 1.4):
		return float(0.2)
	elif (t >= 1.4 and t < 1.6):
		return float(1.6 - t)
	return float(0)

# Derivada función de variable de exitación externa c'(t)
def c_prima(t):
	if (t < 1):
		return 0
	if (t < 1.2 and t > 1):
		return 1
	if (t < 1.4 and t > 1.2):
		return 0
	if (t < 1.6 and t > 1.4):
		return -1
	return 0

# Solución analítica de la ecuación diferencial
def sol_analitica(t):
	return 0.1 - 0.1 * np.cos(((k/m)**0.5)*t)

# Función del sistema sin amortiguación
def f_sin_amortiguar(y, v):
	return ((k/m) * (c_inicial - y)) + ((l_inicial/m) * (0 - v))

# Función del sistema
def f(y, c, v, c_prima, t, l):
	return ((k/m) * (c(t) - y)) + ((l/m) * (c_prima(t) - v))

# Parámetros del sistema
padron = 109932
m = padron / 200
k = 25000
l_inicial = 0
c_inicial = 0.1
t_final = 5
dt = 0.005
# dt = 0.05

def metodo_ponderado_implicito(beta, f, A_inversa, termino_indep, t):
	cantidad = t_final/dt + 1
	# Inicialización de variables
	u2 = np.zeros(int(cantidad*2))
	u2.shape=(2, int(cantidad))
	# u2 = np.zeros_like(t)
	u = np.zeros_like(t)
	v = np.zeros_like(t)

	# Condiciones iniciales (ajustar según el problema)
	u[0] = 0
	v[0] = 0

	# Bucle principal para la integración numérica
	for n in range(1, len(t)-1):
		if beta == 1:
		#      Método ponderado implícito
			aux = np.array(A_inversa @ u2[:,n] + termino_indep)
			u2[:,(n+1)] = aux
		else:
			u[n+1] = u[n] + dt * (beta * v[n+1] + (1-beta) * v[n])
			v[n+1] = v[n] + dt * ((beta * f(u[n+1], v[n+1])) + ((1 - beta)*(f(u[n], v[n]))))
		# if beta == 0:
		# 	u2[n+1] = u2[n] + dt *v[n]
		# 	v[n+1] = v[n] + dt * f(u2[n], v[n])
		# else:
		# #      Método ponderado implícito
		# 	aux = np.array(A_inversa @ u[:,n] + termino_indep)
		# 	u[:,(n+1)] = aux

	if beta == 1:
		return u2[1], t

	return u, t

def solucion_sistema_amortiguado(beta, f, extras, l):
	# Vector de tiempo
	t = np.arange(0, t_final+dt, dt)

	# Inicialización de variables
	u = np.zeros_like(t)
	v = np.zeros_like(t)

	# Condiciones iniciales (ajustar según el problema)
	u[0] = 0
	v[0] = 0

	# Bucle principal para la integración numérica
	for n in range(1, len(t)-1):
		# Método ponderado implícito
		u[n+1] = u[n] + dt * (beta * v[n+1] + (1-beta) * v[n])
		v[n+1] = v[n] + dt * ((beta * f(u[n+1], extras[0], v[n+1], extras[1], n, l)) + ((1 - beta)*(f(u[n], extras[0], v[n], extras[1], n, l))))
	return u, t

def grafico(t, u, beta):
	name = 'solución obtenida con beta = ' + str(beta) + '.png'
	# Graficar la solución
	# ax: tuple[plt.Axes, plt.Axes]
	ax: plt.Axes
	# fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
	fig, ax = plt.subplots()
	ax.plot(t, u, label='aproximación')
	ax.plot(t, sol_analitica(t), 'r--', label='solucion analítica')
	ax.set_title(f'y(t) con beta = {beta} y paso = {dt}')
	ax.set_xlabel('Tiempo (s)')
	ax.set_ylabel('y')
	ax.grid(True)
	plt.savefig(name)

	error_name = f'Error con beta = {str(beta)}.png'

	paso = 0
	aproximacion = np.copy(u)
	for i in range(len(t)):
		aproximacion[i] = (sol_analitica(paso) - aproximacion[i])
		paso += dt
	ax: plt.Axes
	# fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
	fig, ax = plt.subplots()
	ax.plot(t, aproximacion)
	# ax.plot(t, sol_analitica(t), 'r--', label='solucion analítica')
	ax.set_title(f'e(t) con beta = {beta} y paso = {dt}')
	ax.set_xlabel('Tiempo (s)')
	ax.set_ylabel('Error')
	ax.grid(True)
	plt.savefig(error_name)

def grafico_sistema_amortiguado(t, u, h, beta):
	name = 'solución del sistema amortiguado.png'
	# Graficar la solución
	ax: plt.Axes
	#        fig, ax = plt.subplot()
	fig, ax = plt.subplots()
	plt.plot(t, u, label='aproximación')
	plt.title(f'y(t) con beta = {beta} y paso = {h}')
	plt.xlabel('Tiempo (s)')
	plt.ylabel('y')
	plt.grid(True)
	# plt.savefig(name)
	plt.savefig(name)

if __name__ == "__main__":
	# Elección del sistema sin amortiguación
	betas = [0, 0.25, 0.5, 0.75, 1]
	h = dt
	for beta in betas:
		# Matriz inversa para la solución sin amortiguación
		A_inversa = []
		termino_indep = []
		if beta != 0:
			divisor = ((m + h**2)* k * beta**2) + h*beta*l_inicial
			A_inversa = np.array([[(m+(h*beta*l_inicial))/divisor , -1*(h*m*beta)/divisor], 
						[(h*k*beta)/divisor , m/divisor]])
			termino_indep = np.array([[h*(1-beta) , h*(1-beta)*(k/m*c_inicial + l_inicial/m*c_prima(0))]])

		# Vector de tiempo
		t = np.arange(0, t_final+dt, dt)

		u, t = metodo_ponderado_implicito(beta, f_sin_amortiguar, A_inversa, termino_indep, t)
		print("soulciones:", u)
		grafico(t, u, beta)

	# Prueba del sistema amortiguado con beta = 0.5 y dt = 0.005
	l = 750
	extras = [c, c_prima]
	u2, t2 = solucion_sistema_amortiguado(beta, f, extras, l)
	print("soulciones:", u2)
	grafico_sistema_amortiguado(t2, u2, 0.005, 0.5)

	#        u1, t = metodo_ponderado_implicito(0)
	#        print("soulciones:", u1)
	#        grafico(t, u1, 0)
	#        u2, t = metodo_ponderado_implicito(0.25)
	#        print("soulciones:", u2)
	#        grafico(t, u2, 0.25)
	#        u3, t = metodo_ponderado_implicito(0.5)
	#        print("soulciones:", u3)
	#        grafico(t, u3, 0.5)
	#        u4, t = metodo_ponderado_implicito(0.75)
	#        print("soulciones:", u4)
	#        grafico(t, u4, 0.75)
	#        u5, t = metodo_ponderado_implicito(1)
	#        print("soulciones:", u5)
	#        grafico(t, u5, 1)