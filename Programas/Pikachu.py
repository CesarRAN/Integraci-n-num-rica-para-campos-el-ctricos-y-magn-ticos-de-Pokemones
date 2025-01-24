# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:53:45 2024

@author: serpi
"""

import numpy as np
import matplotlib.pyplot as plt

# Constantes
k = 8.99e9  # Constante de Coulomb en N·m²/C²
q_total = 1e-9  # Carga total de la media esfera en Coulombs (1 nC por ejemplo)
R = 1.0  # Radio de la esfera (en metros)
# Función para definir el campo vectorial

def campo_vectorial(x, y, z):
    """
    Devuelve las componentes del campo vectorial en función de las coordenadas (x, y, z).
    Modifica las expresiones según el campo que desees representar.
    """
    condicion_1 = np.sqrt((x-3)**2 + y**2 + (2*z**2)) < 1
    condicion_2 = np.sqrt((x+3)**2 + y**2 + (2*z**2)) < 1
    
    
    f_x = np.where(condicion_1, 1, np.where(condicion_2, -1, 0))  # (0 si condición, -y si no)
    f_y = np.where(condicion_1, 1, np.where(condicion_2, -1, 0))    # (0 si condición, x si no)
    f_z = np.where(condicion_1, 1, np.where(condicion_2, -1, 0))    # (0 si condición, z si no)
    return f_x, f_y, f_z

def campo_vectorial_1(x, y, z):
    """
    Devuelve las componentes del campo vectorial en función de las coordenadas (x, y, z).
    Modifica las expresiones según el campo que desees representar.
    """
    if np.sqrt((x-3)**2 + y**2 + (2*z**2)) < 1:
        return 1
    if np.sqrt((x+3)**2 + y**2 + (2*z**2)) < 1:
        return -1
    else:
        return 0

def campo_vectorial_2(x, y, z, v):
    """
    Devuelve las componentes del campo vectorial en función de las coordenadas (x, y, z).
    Modifica las expresiones según el campo que desees representar.
    """
    condicion_1 = np.sqrt((x-3)**2 + y**2 + (2*z**2)) < 1
    condicion_2 = np.sqrt((x+3)**2 + y**2 + (2*z**2)) < 1
    
    f_x = np.where(condicion_1, v, np.where(condicion_2, -v, 0))  # (0 si condición, -y si no)
    f_y = np.where(condicion_1, v, np.where(condicion_2, -v, 0))    # (0 si condición, x si no)
    f_z = np.where(condicion_1, v, np.where(condicion_2, -v, 0))    # (0 si condición, z si no)
    return f_x, f_y, f_z

def calcular_campo_electrico(x,y,z, densidad, rho_func):
    """
    Calcula el campo eléctrico a partir de una densidad de carga.

    Parameters:
    - rango: Tupla de tres valores (min, max) para el rango de los ejes x, y, z.
    - densidad: El número de puntos en cada dimensión (tamaño de la grilla).
    - rho_func: Función que devuelve la densidad de carga en un punto dado (x, y, z).

    Returns:
    - Ex, Ey, Ez: Componentes del campo eléctrico en la grilla de puntos.
    """
    # Constantes
    epsilon_0 = 8.854187817e-12  # Permitivdad del vacío en F/m (faradios por metro)
    k_e = 1 / (4 * np.pi * epsilon_0)  # Constante de Coulomb
    
    # Crear una grilla de puntos en el rango especificado
    X, Y, Z = np.meshgrid(x, y, z)
    
    # Inicializar las componentes del campo eléctrico
    Ex = np.zeros_like(X)
    Ey = np.zeros_like(Y)
    Ez = np.zeros_like(Z)
    
    # Integración para calcular el campo eléctrico en cada punto
    for i in range(densidad):
        for j in range(densidad):
            for k in range(densidad):
                # Coordenadas del punto de observación
                r = np.array([X[i,j,k], Y[i,j,k], Z[i,j,k]])

                # Iteramos sobre todos los puntos en la grilla de carga
                for xi in range(densidad):
                    for yi in range(densidad):
                        for zi in range(densidad):
                            # Coordenadas del punto de carga
                            r_prime = np.array([x[xi], y[yi], z[zi]])
                            # Vector de desplazamiento
                            r_vec = r - r_prime
                            r_magn = np.linalg.norm(r_vec)
                            
                            if r_magn != 0:
                                # Densidad de carga en el punto de carga
                                rho = rho_func(x[xi], y[yi], z[zi])
                                
                                # Contribución del campo eléctrico en este punto de carga
                                Ex[i,j,k] += k_e * rho * r_vec[0] / r_magn**3
                                Ey[i,j,k] += k_e * rho * r_vec[1] / r_magn**3
                                Ez[i,j,k] += k_e * rho * r_vec[2] / r_magn**3
    
    return Ex, Ey, Ez


def calcular_campo_magnetico(rho_func, v, x, y, z,densidad):
    """
    Calcula el campo magnético generado por una densidad de carga en movimiento.
    
    Parámetros:
    - rho_func: función que recibe (x, y, z) y devuelve la densidad de carga ρ en ese punto.
    - v: velocidad de la densidad de carga (array de 3 elementos: vx, vy, vz).
    - X, Y, Z: grilla de puntos 3D donde se evalúa el campo magnético.
    - rango: rango [min, max] de la región de integración.
    - densidad: densidad de puntos en la grilla de integración.
    
    Retorna:
    - Bx, By, Bz: componentes del campo magnético en cada punto de la grilla (X, Y, Z).
    """

    X, Y, Z = np.meshgrid(x, y, z)
    
    # Coordenadas de integración aplanadas
    x_prime = X.ravel()
    y_prime = Y.ravel()
    z_prime = Z.ravel()
    
    # Vector de corriente J = ρ * v
    Jx = rho_func(x_prime, y_prime, z_prime,v[0])
    Jy = rho_func(x_prime, y_prime, z_prime,v[1])
    Jz = rho_func(x_prime, y_prime, z_prime,v[2])
    
    # Constante μ0 / (4π)
    mu0 = 4 * np.pi * 1e-7
    prefactor = mu0 / (4 * np.pi)
    
    # Inicializar componentes del campo magnético
    Bx = np.zeros_like(X)
    By = np.zeros_like(Y)
    Bz = np.zeros_like(Z)
    
    # Evaluar el campo magnético en cada punto de la grilla (X, Y, Z)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            for k in range(X.shape[2]):
                # Punto de observación
                rx, ry, rz = X[i, j, k], Y[i, j, k], Z[i, j, k]
                
                # Diferencias r - r'
                dx = rx - x_prime
                dy = ry - y_prime
                dz = rz - z_prime
                r_squared = dx**2 + dy**2 + dz**2
                r_cubed = np.sqrt(r_squared)**3
                
                # Evitar división por cero
                r_cubed[r_cubed == 0] = np.inf
                
                # Producto cruz J × (r - r')
                cross_x = Jy * dz - Jz * dy
                cross_y = Jz * dx - Jx * dz
                cross_z = Jx * dy - Jy * dx
                
                # Contribución al campo magnético
                Bx[i, j, k] = np.sum(prefactor * cross_x / r_cubed)
                By[i, j, k] = np.sum(prefactor * cross_y / r_cubed)
                Bz[i, j, k] = np.sum(prefactor * cross_z / r_cubed)
    
    return Bx, By, Bz




# Función para graficar el campo vectorial
def graficar_campo_vectorial(campo_func, rango=(-2, 2), densidad=10, R=1):
    """
    Grafica el campo vectorial en 3D.
    
    :param campo_func: Función que define las componentes del campo vectorial.
    :param rango: Tupla con los límites del espacio (xmin, xmax).
    :param densidad: Número de puntos por eje.
    """
    # Crear una grilla de puntos en el rango especificado
    x = np.linspace(rango[0], rango[1], densidad)
    y = np.linspace(rango[0], rango[1], densidad)
    z = np.linspace(rango[0], rango[1], densidad)
    X, Y, Z = np.meshgrid(x, y, z)
    V = [0,1,0]
    # Calcular las componentes del campo en cada punto de la grilla
    # F_x, F_y, F_z = campo_func(campo_vectorial_2,V,x,y,z,densidad)
    F_x, F_y, F_z = campo_func(x,y,z,densidad,campo_vectorial_1)
    
    # Graficar el campo vectorial
    fig = plt.figure(figsize=(15, 15), dpi=150)
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(X, Y, Z, F_x, F_y, F_z, length=1, normalize=True, color='green', linewidth=0.6)

    # Configurar los ejes
    ax.set_xlim(rango)
    ax.set_ylim(rango)
    ax.set_zlim(rango)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Campo Electrica, Malla: {densidad}',fontsize=30)
    
    # Crear la malla de puntos en coordenadas esféricas
    theta = np.linspace(0, np.pi/2, 300)  # Ángulo polar
    phi = np.linspace(0, 2 * np.pi, 300)  # Ángulo azimutal
    
    # Crear las matrices de coordenadas
    theta, phi = np.meshgrid(theta, phi)

    # Parametrización esférica para obtener las coordenadas cartesianas
    x = R * np.sin(theta) * np.cos(phi)+3
    y = R * np.sin(theta) * np.sin(phi)
    z = R/5 * np.cos(theta)

    # Graficar la esfera
    ax.plot_surface(x, y, z, color='b', alpha=0.6, rstride=5, cstride=5, linewidth=0)
    
    # Crear la malla de puntos en coordenadas esféricas
    theta = np.linspace(np.pi/2,np.pi, 300)  # Ángulo polar
    phi = np.linspace(0, 2 * np.pi, 300)  # Ángulo azimutal
    
    # Crear las matrices de coordenadas
    theta, phi = np.meshgrid(theta, phi)

    # Parametrización esférica para obtener las coordenadas cartesianas
    x = R * np.sin(theta) * np.cos(phi) +3
    y = R * np.sin(theta) * np.sin(phi)
    z = R/5 * np.cos(theta)

    # Graficar la esfera
    ax.plot_surface(x, y, z, color='b', alpha=0.6, rstride=5, cstride=5, linewidth=0)

    # Crear la malla de puntos en coordenadas esféricas
    theta = np.linspace(0, np.pi/2, 300)  # Ángulo polar
    phi = np.linspace(0, 2 * np.pi, 300)  # Ángulo azimutal
    
    # Crear las matrices de coordenadas
    theta, phi = np.meshgrid(theta, phi)

    # Parametrización esférica para obtener las coordenadas cartesianas
    x = R * np.sin(theta) * np.cos(phi)-3
    y = R * np.sin(theta) * np.sin(phi)
    z = R/3 * np.cos(theta)

    # Graficar la esfera
    ax.plot_surface(x, y, z, color='r', alpha=0.6, rstride=5, cstride=5, linewidth=0)
    
    # Crear la malla de puntos en coordenadas esféricas
    theta = np.linspace(np.pi/2,np.pi, 300)  # Ángulo polar
    phi = np.linspace(0, 2 * np.pi, 300)  # Ángulo azimutal
    
    # Crear las matrices de coordenadas
    theta, phi = np.meshgrid(theta, phi)

    # Parametrización esférica para obtener las coordenadas cartesianas
    x = R * np.sin(theta) * np.cos(phi) -3
    y = R * np.sin(theta) * np.sin(phi)
    z = R/3 * np.cos(theta)

    # Graficar la esfera
    ax.plot_surface(x, y, z, color='r', alpha=0.6, rstride=5, cstride=5, linewidth=0)


    # Mostrar la gráfica
    plt.show()

# Usar las funciones para graficar el campo
graficar_campo_vectorial(calcular_campo_electrico, rango=(-5, 5), densidad=10)
