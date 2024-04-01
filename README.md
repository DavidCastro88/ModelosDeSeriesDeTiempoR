# **ModelosDeSeriesDeTiempoR**

![image](https://github.com/DavidCastro88/ModelosDeSeriesDeTiempoR/assets/91480088/d30e6fe1-6646-4077-a632-41154951725f)

**Componentes de una serie:**
- Tendencia
- Estacionalidad
- Ciclos (No se pueden modelar)
- Error

**Combinación de las componentes:**
- Aditiva =  Tendencia + Estacionalidad + Error
- Multiplicativo
    - Completamente multiplicativo: Tendencia x Estacionalidad x exp(Error)  =>  log(Tendencia) + log(Estacionalidad) + (Error)
    - Parcialmente multiplicativo: Tendencia x Estacionalidad + exp(Error)

## Pasos para el análisis y modelación de una serie temporal
![image](https://github.com/DavidCastro88/ModelosDeSeriesDeTiempoR/assets/91480088/cd0cc62c-7283-4421-a1d2-ef74ffcc5f2c)


##  **Modelos**
### **1. Modelos globales:**

**NOTA**: Si no tienen componente estacional solo se trabaja con la componente de tendencia.

**Para series aditivas:**
- Polinomial de grado p, estacional con indicadoras nivel de referencia último periodo del año (Modelo Aditivo).
- Polinomial de grado p, estacional con trigonométricas en k frecuencias Fj, (Modelo aditivo).
  
**Para series multiplicativas:**

- Log polinomial de grado p, estacional con indicadoras nivel de referencia último periodo del año (Modelo completamente multiplicativo)
- Log polinomial de grado p, estacional con trigonométricas en k frecuencias Fj, (Modelo completamente multiplicativo).
- Exponencial polinomial grado p, estacional con indicadoras nivel de referencia último periodo del año (Modelo parcialmente multiplicativo).
- Exponencial polinomial de grado p, estacional con trigonométricas en k frecuencias Fj, (Modelo parcialmente multiplicativo).
     
### **2. Modelos locales:**
**NOTA**: Si no tienen componente estacional solo se trabaja con la componente de tendencia.

**Para series aditivas:**
- Modelo de tendencia local de cambio de nivel y de pendientes que evolucionan lentamente con t, aditiva a factor estacional con efectos   estacionales que evolucionan lentamente en el tiempo. **(Suavizamiento Exponencial Holt Winters aditivo)**
- Modelo de tendencia local lineal aditiva a factor estacional con efectos estacionales constantes en el tiempo. **(Descomposición aditiva &
loess lineal)**
- Modelo de tendencia local cuadrática aditiva a factor estacional con efectos estacionales constantes en el tiempo. **(Descomposición aditiva & loess cuadrático)**

**Para series multiplicativas:**
- Modelo de tendencia local de cambio de nivel y de pendientes que evolucionan lentamente con t, multiplicativa a factor estacional con efectos estacionales que evolucionan lentamente en el tiempo: **(Suavizamiento Exponencial Holt Winters multiplicativo)**
- Modelo de tendencia local lineal multiplicativa a factor estacional con efectos estacionales constantes en el tiempo: **(Descomposición mutiplicativa & loess lineal)**
- Modelo de tendencia local cuadrática multiplicativa a factor estacional con efectos estacionales constantes en el tiempo: **(Descomposición multiplicativa & loess cuadrático)**

#### **NOTA**: Orden de prioridad para elegir el modelo
1. Validación de suspuestos de los residuos de ajuste:
    - Media cero. (Grafica Residuales vs Tiempo bien centrada)
    - Varianza constante. (Que no hayan patrones tipo cono o rombo en los Residuos vs ajustados)
    - Independencia de los errores (Que no haya ciclos ni rachas en los Residuos vs Tiempo)
    - Carencia de ajuste (Que no hayan patrones tipo U, W, V en los Residuos vs Ajustados)
2. Calidad de pronóstico, MSE, MAE, MAPE, RMSE
3. Calidad de ajuste, AIC y BIC.
   
### **3. Modelos para el error**
**NOTA**: Utilizar función armasubsets(), para hallar los valores de p,q.

Tanto para series aditivas o multiplicativas usar los mejores modelos globales, y modelar el error con modelos estructurales AR (Modelo Autorregresivo) y MA (Modelo de Media Móvil).

Validar que los residuos de ajuste (at) sean un Ruido Blanco.
- Independencia (Durbin Watson -  Ljung Box)
- Normalidad (Shapiro Wilk)
- Autocorrelación (ACF - PACF)
- Media cero

![image](https://github.com/DavidCastro88/ModelosDeSeriesDeTiempoR/assets/91480088/fdea56e6-e99e-4aeb-a04f-96e6b13ae40a)

![image](https://github.com/DavidCastro88/ModelosDeSeriesDeTiempoR/assets/91480088/9716ea25-ee71-4682-85c5-14f8bd15d254)

**Funciones para identificar ARMA´s en R**

- auto.arima() de la librería forecast.
- autoarmafit() de la librearía timesac
- SelectModel() de la librearía FitAR
- Recomendada armasubsets() de la librearia TSA

### **4. Modelos ARIMA(p,d,q)(P,D,Q)[s]**:

    Pasos:
        1. Estabilizar varianza, si fuese necesario (en series multiplicativas se modela en la escala logaritmo natural).
        2. Si existe tendencia y/o raız unitaria regular, aplicar  ́ ∇d, variando orden d progresivamente desde 1, hasta que el nivel de la serie filtrada sea estable y la ACF muestral indique proceso ergodico en la parte regular (en los primeros k la ACF muestral decae rapido).
        3. Si hay componente St periodica exacta (o casi exacta) y/o con raıces unitarias estacionales, de periodo s, aplicar ∇Ds , variando el orden D progresivamente desde 1, hasta que no se observe el patron periódico exacto y  ́la ACF muestral indique proceso ergodico en la parte estacional (en los k múltiplos de  ́ s la ACF muestral decae rapido).  ́
        4. Una vez hallados d y D verifique que el proceso filtrado es estacionario, es decir: tiene media constante, varianza constante y la ACF parte regular y estacional, indica que el proceso es ergodico.
        5. Sobre la serie debidamente diferenciada por tendencia y/o estacionalidad, aplique metodos para identificar el modelo ARMA estacional estacionario: ACF-PACF, funcion R armasubsets, para luego definir el ARIMA sobre la serie sin diferenciar.
        6. Si los ordenes d y D son distintos de cero, los modelos ARIMA seran ́sin deriva, pero si uno de estos ordenes es cero, tal vez sea necesario incluir una deriva en los modelos identificados y para resolver esto, analice la media muestral de la serie debidamente diferenciada: si esta es muy pequena no se incluye la deriva pero si pudiera ser significativa, se incluye la deriva.
        7. NOTA: La funcion R auto.arima podemos aplicarla directamente sobre la serie sin haberla diferenciado; esta funcion tiene la capacidad de realizar pruebas de raıces unitarias regular y estacionales y por tanto puede identificar modelos ARIMA(p,d,q)(P,D,Q)[s] sin deriva y con deriva.
        
![image](https://github.com/DavidCastro88/ModelosDeSeriesDeTiempoR/assets/91480088/0281d5d4-d911-4b9f-b241-0b79eb474ad4)


