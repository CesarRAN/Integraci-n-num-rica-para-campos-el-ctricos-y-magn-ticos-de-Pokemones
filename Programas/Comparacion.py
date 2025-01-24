# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:44:12 2024

@author: serpi
"""

import numpy as np
import matplotlib.pyplot as plt

# Constantes
epsilon_0 = 8.854187817e-12  # Permisividad eléctrica en el vacío (F/m)
R = 1.0                      # Radio de la esfera (m)
rho = 1e-6                   # Densidad volumétrica de carga (C/m^3)

# Punto de evaluación del campo eléctrico
r_eval = np.linspace(0.1, 2.0, 50)  # Distancias desde el centro (m)

# Función analítica
def campo_analitico(r):
    E = np.zeros_like(r)
    dentro = r <= R
    fuera = r > R
    E[dentro] = (rho * r[dentro]) / (3 * epsilon_0)
    E[fuera] = (rho * R**3) / (3 * epsilon_0 * r[fuera]**2)
    return E

# Integración numérica
def campo_numerico(r_eval, N):
    # Crear una malla dentro de la esfera
    x = np.linspace(-R, R, N)
    y = np.linspace(-R, R, N)
    z = np.linspace(-R, R, N)
    X, Y, Z = np.meshgrid(x, y, z)
    r_i = np.sqrt(X**2 + Y**2 + Z**2)
    
    # Filtrar puntos dentro de la esfera
    inside = r_i <= R
    X, Y, Z = X[inside], Y[inside], Z[inside]
    r_i = r_i[inside]
    
    # Volumen diferencial
    dx = x[1] - x[0]
    dV = dx**3
    
    # Campo eléctrico numérico en cada punto de evaluación
    E_numerico = []
    for r in r_eval:
        # Posición del punto de evaluación
        rx, ry, rz = r, 0, 0  # Considerar puntos en el eje x
        r_vec = np.array([rx - X, ry - Y, rz - Z])
        r_mag = np.sqrt(r_vec[0]**2 + r_vec[1]**2 + r_vec[2]**2)
        
        # Evitar divisiones por cero
        r_mag[r_mag == 0] = np.inf
        
        # Suma de contribuciones de todos los elementos de volumen
        E = (1 / (4 * np.pi * epsilon_0)) * np.sum(
            (rho * r_vec / r_mag**3) * dV,
            axis=1
        )
        E_numerico.append(np.linalg.norm(E))
    
    return np.array(E_numerico)

Error = []
finura = [10,20,50,100]



for i in finura:
    # Graficar los resultados
    plt.figure(figsize=(10, 10), dpi=800)
    # Valores de r para comparar los campos
    r_values = np.linspace(0, 2 * R, 100)
    #campo_analitico_vals = np.array([campo_analitico(r) for r in r_values])
    #campo_numerico_vals = np.array([campo_numerico(r,i) for r in r_values])
    campo_analitico_vals = campo_analitico(r_values)
    campo_numerico_vals = campo_numerico(r_values,i)
    # Evitamos divisiones por cero: ignoramos valores analíticos iguales a cero
    diferencias_relativas = np.where(
        campo_analitico_vals != 0,
        np.abs(campo_analitico_vals - campo_numerico_vals) / campo_analitico_vals * 100,
        0  # Si campo_analitico_vals es 0, la diferencia relativa se considera 0
    )
    
    # Calculamos el promedio de las diferencias relativas (en porcentaje)
    promedio_diferencia_porcentual = np.mean(diferencias_relativas)
    Error.append(promedio_diferencia_porcentual)
    print(f"El promedio de la diferencia relativa entre los campos es: {promedio_diferencia_porcentual:.2f}%")
    
    plt.plot(r_values, campo_analitico_vals, label="Analítico", linewidth=2)
    plt.plot(r_values, campo_numerico_vals, 'o', label="Numérico", markersize=5, color='g')
    plt.axvline(R, color='k', linestyle='--', label="Radio de la esfera")
    plt.xlabel("Distancia desde el centro (m)")
    plt.ylabel("Campo eléctrico (N/C)")
    if i !=20:
        plt.title(f"Campo eléctrico: solución analítica vs numérica d: {i}" )
    else:
        a=15
        plt.title(f"Campo eléctrico: solución analítica vs numérica d: {a}" )
    plt.legend()
    plt.grid()
    plt.show()
