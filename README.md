# Integracion-numerica-para-campos-electricos-y-magneticos-de-Pokemones
Este proyecto estudia los campos eléctricos y magnéticos generados por distribuciones de carga en movimiento
Los pokemones considerados son:

- **Electrode**: representado como una esfera definida por la ecuación  
  
  $$x^2 + y^2 + z^2 = 1$$
  - Para modelar el campo eléctrico y magnético, se considerará una distribución de carga uniforme positiva en una mitad de la esfera y una distribución de carga negativa en la otra mitad.
  - En este caso, se obtendrá analíticamente el campo eléctrico $\mathbf{E}$ y el campo magnético $\mathbf{B}$ y se comparará la solución analítica con la obtenida por simulación numérica, suponiendo que el pokemon se mueve a una velocidad constante $\mathbf{v}$.

  
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
  - En este caso, el campo eléctrico $\mathbf{E}$ y el campo magnético $\mathbf{B}$ serán obtenidos exclusivamente mediante simulación numérica, suponiendo nuevamente un movimiento a velocidad constante $\mathbf{v}$.

  **Figura 3. Pikachu**  
![pika](https://github.com/user-attachments/assets/1fa839db-6b22-46ea-ba32-8f180b7c5a34)

  **Figura 4. Densidad de carga de Pikachu**  
![GEOGEBRA_2](https://github.com/user-attachments/assets/db76d504-27de-4bfd-891f-2ba3170986a0)


  
$$\rho(\vec{r}) = 
  \begin{cases}
  Q & \text{si } (x + 5)^2 + y^2 + 5z^2 < 1 \\
  -Q & \text{si } (x - 5)^2 + y^2 + 5z^2 < 1 \\
  0 & \text{si otro caso}
  \end{cases}$$

---

# Metodología

1. **Electroid**  
   - Se obtendrá una solución analítica para los campos $\mathbf{E}$ y $\mathbf{B}$ de acuerdo con la configuración de carga descrita.  
   - Se realizará una simulación numérica del campo eléctrico y magnético de Electroid, considerando su movimiento a una velocidad $\mathbf{v}$.  
   - Finalmente, se compararán los resultados de la simulación numérica con la solución analítica para verificar la precisión de la simulación.
  

2. **Pikachu**  
   - Se implementará una simulación numérica para calcular el campo eléctrico y magnético generados por la configuración de carga en las mejillas, modeladas por los dos elipsoides mencionados.  
   - Se asumirá un movimiento constante del pokemon a velocidad $\mathbf{v}$.


---

# Algoritmo para la integración numérica del campo eléctrico y magnético

## Algoritmo 1: Cálculo del campo eléctrico mediante discretización

**Datos de entrada:**  
- Volumen de integración $V$, dividido en $N$ celdas pequeñas con volúmenes $\Delta V_i$ y centros $\vec{r}_i$.

**Resultado:**  
- Campo eléctrico $\vec{E}(\vec{r})$ en un punto $\vec{r}$.

1. Para cada celda $i = 1 \dots N$:
   - Aproximar $\rho(\vec{r}')$ como constante dentro de la celda:  
$$\rho(\vec{r}') \approx \rho(\vec{r}_i).$$

2. Inicializar $\vec{E} \gets \vec{0}$.

3. Para cada celda $i = 1 \dots N$:
   - Calcular el vector de diferencia:  
$$\vec{d}_i \gets \vec{r} - \vec{r}_i. $$
   - Calcular la norma:  
$$d_i \gets |\vec{d}_i|.$$
   - Actualizar $\vec{E}$:  
$$\vec{E} \gets \vec{E} + \frac{\rho(\vec{r}_i)}{4 \pi \epsilon_0} \frac{\vec{d}_i}{d_i^3} \Delta V_i.$$




## Algoritmo 2: Cálculo del campo magnético mediante discretización (Ley de Biot-Savart)

**Datos de entrada:**  
- Volumen $V$ dividido en $N$ celdas con volúmenes $\Delta V_i$ y centros $\vec{r}_i$.

**Resultado:**  
- Campo magnético $\vec{B}(\vec{r})$ en un punto $\vec{r}$.


1. Para cada celda $i = 1 \dots N$:
   - Aproximar $\rho(\vec{r}') \approx \rho(\vec{r}_i)$ (constante dentro de la celda).
   - Aproximar $\vec{v}(\vec{r}') \approx \vec{v}(\vec{r}_i)$ (constante dentro de la celda).

2. Inicializar $\vec{B} \gets \vec{0}$.

3. Para cada celda $i = 1 \dots N$:
   - Calcular el vector de diferencia:  
$$\vec{d}_i \gets \vec{r} - \vec{r}_i$$
   - Calcular la norma:  
$$d_i \gets |\vec{d}_i|.$$
   - Calcular el producto cruzado:  
$$\vec{C}_i \gets \vec{v}(\vec{r}_i) \times \frac{\vec{d}_i}{d_i^3}.$$
   - Actualizar $\vec{B}$:  
$$\vec{B} \gets \vec{B} + \frac{\mu_0}{4\pi} \, \rho(\vec{r}_i) \, \vec{C}_i \, \Delta V_i.$$
