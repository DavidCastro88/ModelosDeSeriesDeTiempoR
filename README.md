# **ModelosDeSeriesDeTiempoR**

![image](https://github.com/DavidCastro88/ModelosDeSeriesDeTiempoR/assets/91480088/d30e6fe1-6646-4077-a632-41154951725f)

Componentes de una serie:
- Tendencia
- Estacionalidad
- Ciclos (No se pueden modelar)
- Error

Combinación de las componentes:
- Aditiva =  Tendencia + Estacionalidad + Error
- Multiplicativo
    - Completamente multiplicativo: Tendencia x Estacionalidad x exp(Error)  =>  log(Tendencia) + log(Estacionalidad) + (Error)
    - Completamente multiplicativo: Tendencia x Estacionalidad + exp(Error)

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
