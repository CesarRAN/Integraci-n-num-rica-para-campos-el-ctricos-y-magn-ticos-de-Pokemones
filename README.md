# Integracion-numerica-para-campos-electricos-y-magneticos-de-Pokemones
Este proyecto estudia los campos eléctricos y magnéticos generados por distribuciones de carga en movimiento
Los pokemones considerados son:

- **Electrode**: representado como una esfera definida por la ecuación  
  
  $$x^2 + y^2 + z^2 = 1$$
  - Para modelar el campo eléctrico y magnético, se considerará una distribución de carga uniforme positiva en una mitad de la esfera y una distribución de carga negativa en la otra mitad.
  - En este caso, se obtendrá analíticamente el campo eléctrico \(\mathbf{E}\) y el campo magnético \(\mathbf{B}\) y se comparará la solución analítica con la obtenida por simulación numérica, suponiendo que el pokemon se mueve a una velocidad constante \(\mathbf{v}\).

  
  **Figura 1. Electrode**  
![electrode](https://github.com/user-attachments/assets/ed4239af-b81b-45a3-bd4c-6f82cf0f10a3)

  **Figura 2. Densidad de carga de Electrode**  
![GEOGEBRA_1](https://github.com/user-attachments/assets/7b595226-ead1-430d-b151-1f39f44a86ff)

  
  $$\rho(\vec{r}) = 
  \begin{cases}
  Q & \text{si } r \in [0,R) \text{ y } \theta \in \left[0,\frac{\pi}{2}\right] \\
  -Q & \text{si } r \in [0,R) \text{ y } \theta \in \left[\frac{\pi}{2},\pi\right] \\
  0 & \text{si otro caso}
  \end{cases}$$
  

- **Pikachu**: su carga estará localizada en sus mejillas rojas, modeladas por dos elipsoides definidos por las ecuaciones  
  
  $$(x - 5)^2 + y^2 + 5z^2 = 1$$  
  
  $$(x + 5)^2 + y^2 + 5z^2 = 1$$  
  - En este caso, el campo eléctrico \(\mathbf{E}\) y el campo magnético \(\mathbf{B}\) serán obtenidos exclusivamente mediante simulación numérica, suponiendo nuevamente un movimiento a velocidad constante \(\mathbf{v}\).

  **Figura 3. Pikachu**  
![pika](https://github.com/user-attachments/assets/1fa839db-6b22-46ea-ba32-8f180b7c5a34)

  **Figura 4. Densidad de carga de Pikachu**  
![Pikachu_electrico](https://github.com/user-attachments/assets/4884d960-3736-4a51-9687-c63d833301ae)

  
  $$\rho(\vec{r}) = 
  \begin{cases}
  Q & \text{si } (x + 5)^2 + y^2 + 5z^2 < 1 \\
  -Q & \text{si } (x - 5)^2 + y^2 + 5z^2 < 1 \\
  0 & \text{si otro caso}
  \end{cases}$$

---

# Metodología

1. **Electroid**  
   - Se obtendrá una solución analítica para los campos \(\mathbf{E}\) y \(\mathbf{B}\) de acuerdo con la configuración de carga descrita.  
   - Se realizará una simulación numérica del campo eléctrico y magnético de Electroid, considerando su movimiento a una velocidad \(\mathbf{v}\).  
   - Finalmente, se compararán los resultados de la simulación numérica con la solución analítica para verificar la precisión de la simulación.

2. **Pikachu**  
   - Se implementará una simulación numérica para calcular el campo eléctrico y magnético generados por la configuración de carga en las mejillas, modeladas por los dos elipsoides mencionados.  
   - Se asumirá un movimiento constante del pokemon a velocidad \(\mathbf{v}\).
