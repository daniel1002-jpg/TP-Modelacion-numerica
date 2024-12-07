import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# Parámetros del sistema
padron = 109932
m = padron / 200
k = 25000
l_inicial = 0
c_inicial = 0.1
t_final = 5
dt = 0.005

# Analizar compresión del sistema
def evaluar_compresion(beta, h, f, l, k, t):
	u, v, t = solucion_sistema_amortiguado(beta, h, f, l, k)
	# print("soulciones:", u)
	grafico_sistema_amortiguado(t, u, v, h, beta, k, l)

	return u, v, t


# Calcular el orden de convergencia
# def calcular_orden(beta, f, hs):
# 	errores = []
# 	for h in hs:
# 		t = np.arange(0, t_final + h, h)

# 		u, _ = metodo_ponderado_implicito(beta, h, f, t)
# 		error = max(np.abs(sol_analitica(t) - u))
# 		errores.append(error)

# 	# Calcular orden de convergencia
# 	# print("Errores:", errores)
# 	ordenes = []
# 	for i in range(len(hs) - 1):
# 		orden = np.log(errores[i+1] / errores[i]) / np.log(hs[i] / hs[i+1])
# 		ordenes.append(orden)

# 	# print("Ordenes:", ordenes)
# 	return errores, ordenes

# Obtener período promedio de oscilación
def obtener_periodo_promedio(aproximacion, t):
	# Encontrar picos (máximos)
	picos, _ = find_peaks(aproximacion)

	# Calcular los tiempos de los picos
	tiempos = t[picos]

	# Calcular período entre picos consecutivos
	T = np.diff(tiempos)

	# Calcualr período promedio
	T_promedio = np.mean(T)

	return T_promedio

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
	return ((k/m) * (c_inicial - y)) + ((l_inicial/m) * (c_prima(0) - v))

# Función del sistema
def f(y,v, t, l):
	return ((k/m) * (c(t) - y)) + ((l/m) * (c_prima(t) - v))

# Obtener matriz inversa para resolver sistema: Ax = b
def calcular_A_inversa(beta, h, l, k):
	divisor = (m + ((h**2) * k * (beta**2)) + (h * beta * l))

	return np.array([[((m + (h * beta * l)) / divisor) , ((h * m * beta) / divisor)], 
			      [(((-h * k * beta)) / divisor) , (m / divisor)]])
	
# Obtnerer término independiente para resolver el sistema
def calcular_termino_indep(beta, f, u, v):
	return np.array([h * ((1 - beta) * v), 
			(h * beta) * (k * (c_inicial/m) +  0 * (0/m)) + h * (1 - beta) * f(u, v)])

# Obtnerer término independiente para resolver el sistema amortiguado
def calcular_termino_indep(beta, f, u, v, t, l, k):
	return np.array([h * ((1 - beta) * v), 
			(h * beta) * (k * (c(t)/m) +  l * (c_prima(t)/m)) + h * (1 - beta) * f(u, v, t, l)])

def euler_explicito(t, h, f):
	# Inicialización de variables
	u = np.zeros_like(t)
	v = np.zeros_like(t)

	# Condiciones iniciales
	u[0] = 0
	v[0] = 0

	# Bucle principal del método
	for n in range(len(t)-1):
		u[n+1] = u[n] + h * v[n]
		v[n+1] = v[n] + h * f(u[n], v[n])

	return u, t


def euler_implicito(t, h, f):
	# Numéro de pasos de tiempo
	cantidad = len(t)
	
	# Inicialización de variables
	u = np.zeros((2, cantidad))
	u[:, 0] = [0, 0]

	# Obtención de A inversa
	A_inversa = calcular_A_inversa(1, h, 0, k)

	# Obtención del término independiente
	termino_indep = calcular_termino_indep(1, f, 0, 0)
	
	# Bucle principal del método
	for n in range(cantidad - 1):
		aux = A_inversa @ u[:,n] + termino_indep
		u[:,(n+1)] = aux
	
	# print(f"resulatdos: u: {u[0]}, v: {u[1]}")
	return u[0], t

def metodo_ponderado_implicito(beta, h, f, t):
	# Inicialización de variables
	u = np.zeros_like(t)
	v = np.zeros_like(t)

	# Condiciones iniciales
	u[0] = 0
	v[0] = 0

	uv_n = np.array([u[0], v[0]])

	if beta == 1:
		return euler_implicito(t, h, f)
	if beta == 0:			
		return euler_explicito(t, h, f)
	else:
		# Calculo A inviversa para resolver el sistema
		A_inversa = calcular_A_inversa(beta, h, 0, k)

		# Bucle principal para la integración numérica
		for n in range(len(t) - 1):
			# Cálculo el término independiente para cada iteración
			termino_indep = calcular_termino_indep(beta, f, u[n], v[n])
			
			b = uv_n + termino_indep

			uv_np1 = A_inversa @ b

			u[n+1], v[n+1] = uv_np1
			uv_n = uv_np1

	return u, t

def solucion_sistema_amortiguado(beta, h, f, l, k):
	# Vector de tiempo
	t = np.arange(0, t_final+dt, dt)

	# Inicialización de variables
	u = np.zeros_like(t)
	v = np.zeros_like(t)

	# Condiciones iniciales
	u[0] = 0
	v[0] = 0

	uv_n = np.array([u[0], v[0]])

	# Calculo A inviversa para resolver el sistema
	A_inversa = calcular_A_inversa(beta, h, l, k)

	# Bucle principal para la integración numérica
	i = 0
	for n in range(len(t) - 1):
		# Cálculo el término independiente para cada iteración
		termino_indep = calcular_termino_indep(beta, f, u[n], v[n], t[n], l, k)
		
		b = uv_n + termino_indep

		uv_np1 = A_inversa @ b

		u[n+1], v[n+1] = uv_np1
		uv_n = uv_np1
		i += dt
	
	return u, v, t

def grafico(t, u, beta):
	name = 'solución obtenida con beta = ' + str(beta) + '.png'
	# Graficar la solución
	# ax: plt.Axes
	ax: tuple[plt.Axes]
	fig, ax = plt.subplots(1, 2, figsize=(10, 5))
	ax[0].plot(t, u, label='aproximación')
	ax[0].plot(t, sol_analitica(t), 'r--', label='solucion analítica')
	ax[0].set_title(f'y(t) con beta = {beta} y paso = {dt}')
	ax[0].set_xlabel('Tiempo (s)')
	ax[0].set_ylabel('y')
	ax[0].grid(True, linestyle='-', linewidth=0.5, alpha=0.5)
	ax[0].autoscale(enable=True, axis="y")
	ax[0].margins(y=0.1)
	# plt.savefig(name)

	error_name = f'Error con beta = {str(beta)}.png'

	paso = 0
	errores = (sol_analitica(t) - u)
	errores = np.copy(u)
	aproximacion = np.copy(u)
	for i in range(int((t_final/dt)+1)):
		errores[i] = (sol_analitica(paso)- aproximacion[i])
		paso += dt

	# Calculo Error máximo
	error_maximo = max(np.abs(sol_analitica(t) - u))

	print(f"Error máximo: {error_maximo}\n")
	name = f"Resultados obtenidos con beta = {beta}.png"
	# Graficar error de la aproximación
	# ax: plt.Axes
	# fig, ax = plt.subplots()
	ax[1].plot(t, errores, label="error")
	ax[1].set_title(f'e(t) con beta = {beta} y paso = {dt}')
	ax[1].set_xlabel('Tiempo (s)')
	ax[1].set_ylabel('Error')
	ax[1].grid(True, linestyle='-', linewidth=0.5, alpha=0.5)
	# ax[1].autoscale(enable=True, axis="both")
	ax[1].margins(y=0.1)

	plt.tight_layout()
	plt.savefig(name)

def grafico_sistema_amortiguado(t, u, v, h, beta, k, l):
	# name = 'solución del sistema amortiguado optimizado.png'
	name = 'acelaración del sistema amortiguado optimizado.png'
	
	# Graficar la solución
	ax: plt.Axes
	fig, ax = plt.subplots()
	# plt.plot(t, u, label='aproximación')

	exitaciones = [c(paso) for paso in t]

	# Obtener máxima compresión
	aproximacion = np.copy(u)
	minimo = 0
	for i in range(len(t)):
		aproximacion[i] = (u[i]- float(c(t[i])))
		if aproximacion[i] < minimo:
			minimo = aproximacion[i]

	print(f"Máxima compresión obtenida con k = {k} y lambda = {l}: {minimo}")

	# plt.plot(t, exitaciones, "r-", label="amortiguación")
	# print(f"acelarión: {v}")
	# print(f"exitación externa: {exitaciones}")
	plt.plot(t, v, "r-", label='aceleración')
	plt.title(f'y\'\' (t) Crank-Nicolson')
	# plt.title("y(t) del amortiguador")
	plt.xlabel('Tiempo (s)')
	plt.ylabel('y\'\'')
	# plt.ylabel('compresión')
	plt.grid(True)
	plt.savefig(name)

def maximaCompresion(beta, h,k,l, t):
	# def f2(y,c,t,y_prim,c_prim):
	# 	m = 104351/200
	# 	res =  (k/m)*(c(t)-y) + (lam/m) * (c_prim(t) - y_prim)
	# 	return res

	paso = 0
	aproximada, _, _ = evaluar_compresion(beta, h, f, l, k, t)
	minimo = float("inf")
	for j in range (len(t)-1):
		aproximada[j] = (aproximada[j] - c(paso))
		if aproximada[j] < minimo :
			minimo = aproximada[j]
		paso += h
	return minimo

if __name__ == "__main__":
	# Elección del sistema sin amortiguación
	# betas = [0, 0.25, 0.5, 0.75, 1]
	# h = dt
	# beta = 1
	# l = 0
	# cp = 0
	
	# Vector de tiempos
	# t = np.arange(0, t_final+dt, dt)
	# for beta in betas:

	# 	print(f"paso usado: {dt}")
	# 	u, t = metodo_ponderado_implicito(beta, h, f_sin_amortiguar, t)
	# 	print(f"soulciones para beta = {beta}: {u}")
	# 	grafico(t, u, beta)

	# 	# Periodo de oscialción
	# 	T = obtener_periodo_promedio(u, t)
	# 	print("Período de oscilación (promedio):", T)

		# Orden de cada caso
		# hs = [0.02, 0.01, 0.005, 0.0025]  # Diferentes pasos
		# errores, ordenes = calcular_orden(beta, f_sin_amortiguar, hs)

		# for i, h in enumerate(hs[:-1]):
		# 	print(f"h = {h}, Error = {errores[i]:.5e}, Orden de Convergencia = {ordenes[i]:.2f}\n\n")
	

	# Prueba del sistema amortiguado con beta = 0.5 y dt = 0.005
	l = 750
	k = 25000
	beta = 0.5
	h = 0.005
	t = np.arange(0, t_final+h, h)
	# # extras = [c, c_prima]
	# u2, v2, t2 = solucion_sistema_amortiguado(beta, h, f, l)
	# print("soulciones:", u2)
	# grafico_sistema_amortiguado(t2, u2, h, beta, k, l)

	# Búsqueda de lambda y k óptimos
	# Primer método
	# for k_test in [2500, 5000, 10000, 15000, 25000, 30000, 45000, 50000]:
	# 	evaluar_compresion(beta, h, f, l, k_test, t)

	# for l_test in [500, 750, 1000, 1250, 1500, 1750, 2000, 2500]:
	# 	evaluar_compresion(beta, h, f, l_test, 10000, t)

	# for l_test in [500, 750, 1000, 1250, 1500, 1750, 2000, 2500]:
	# 	evaluar_compresion(beta, h, f, l_test, 25000, t)	

	# evaluar_compresion(beta, h, f, l, k, t)
	# evaluar_compresion(beta, h, f, l, 2500, t)
	# evaluar_compresion(beta, h, f, l, 5000, t)
	# evaluar_compresion(beta, h, f, l, 10000, t)
	# evaluar_compresion(beta, h, f, l, 15000, t)
	# evaluar_compresion(beta, h, f, l, 40000, t)


	# Segundo método
	# maximaComp = -0.05
	# compFinal = float("inf")
	# minimaPonderacion = 10000000000
	# lamElecto = 150
	# kElecto = 25000
	# lamV = np.arange(150,12000,50)
	# kV = np.arange(2500,100000,1000)
	# for i in range(int(lamV.shape[0])-1):
	# 	for j in range(int(kV.shape[0])-1):
	# 		act = 0
	# 		act = maximaCompresion(beta, h, int(kV[j]), int(lamV[i]), t)
	# 		if act >= maximaComp:
	# 			ponderacionActual = lamV[i]/750 + kV[j]/25000
	# 			if ponderacionActual < minimaPonderacion :
	# 				minimaPonderacion = ponderacionActual
	# 				lamElecto= lamV[i]
					
	# 				kElecto = kV[j]
	# 				compFinal = act
	# 		print(f"Máxima compresión k = {kV[j]} lamda = {lamV[i]}: {act}")
	# print('K elegido = '+ str(kElecto) +' y lambda electo = ' + str(lamElecto) + ' con compresion = ' + str(compFinal))
	kElecto = 62500
	lamElecto = 8400
	t=np.arange(0,t_final+h,h)

	evaluar_compresion(0.5, h, f, lamElecto, kElecto, t)