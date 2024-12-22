# TP: Análisis de la suspensión de un vehículo: Sistema oscilatorio amortiguado

---

## Objetivo del trabajo
- El objetivo de este trabajo es analizar un modelo matemático que simula un sistema de amortiguación de un vehículo con un oscilador amortiguado. Para eso se utilizará el método ponderado implícito variando los datos del problema para hallar un resultado que ajuste mejor a las condiciones impuestas en el enunciado.
---
## Gráficos
La carpeta gráficos contiene los gráficos de las soluciones y errores de truncamiento obtenidos con los métodos numéricos pedidos en el enunciado
 
### Generación de gráficos
   En caso de que se quisiera probar el sistema con diferentes datos a los usados en el código principal se necesitarán tener las siguientes dependencias:
  - **Python 3.12+**
  - `matplotlib`
  - `numpy`
---
## Ejecución 
Para probar el sistema simulado basta con correr la siguiente línea de código `python metodo_ponderado_implicito.py` en la terminal, o también puede ejecutarse desde el IDE que se utilizando. 
Para probar la simulación se puede descomentar la sección que se desee ejecutar en el main.
- Secciones:
  - **Sistema sin amortiguación:** este sección ejecutara el método numérico usando los distintos valores para `beta` y los de pasos de tiempo (h) cargados. una vez ejecutados todos los casos se generarán automáticamente los gráficos correspondientes con los valores de `beta` y `h` usados.
  - **Sistema amortiguado:** está sección, a diferencia de la descrita arriba, ejecutará el método numérico utilizando los valores de `beta` y `h` definidos cómo óptimos para resolución del problema, y después generará automáticamente el grafico respectivo a la simulación probada.
  - **Obtención de `lambda` y `k` óptimos:** Esta sección evaluará la compresión del sistema con distintos valores de `k` y `lambda` ,dentro de un rango definido para cada variable, hasta encontrar la combinación que se ajuste mejor a las condiciones indicadas en el enunciado. Luego se generará un gráfico con el resultado final usando los datos determinados a lo largo del análisis.
