#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(forecast)
#library(FitAR)
library(stats)
library(car)
library(carData)
library(TSA)
library(fANCOVA)
library(TSA)
library(pdR)
library(timsac);library(lmtest) 
library(uroot)
library(prophet)
library(MLmetrics)
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")

#Data Embalses
Data=read.table(file.choose(),header=T,sep=";")

obtener_serie <- function(nombre_serie,anio_inicio,mes_inicio){
  serie=Data[Data$serie == nombre_serie, ]
  rownames(serie) <- NULL
  yt=ts(serie$y,freq=365,start=c(anio_inicio,mes_inicio))
  return(yt)
}

aplicar_filtros_suavizador <- function(serie,nombre,anio_inicio,mes_inicio,anio_pron,mes_pron) {
  m=5
  k=24
  n=length(serie)-m
  t=1:n
  yt_train=ts(serie[t],frequency=365,start=c(anio_inicio,mes_inicio))
  t_test=(n+1):length(serie)
  yt_test=ts(serie[t_test],freq=1,start=c(anio_pron,mes_pron,27))
  w5=c(-0.073,0.294,0.558,0.294,-0.073) #Vector de pesos para Filtros Henderson
  
  #Ajuste de los modelos
  mod1 = HoltWinters(yt_train,beta=FALSE,gamma=FALSE)
  mod2 <- filter(yt_train, rep(1/(k+1),k+1), sides = 1,circular = T)
  mod3 <- filter(yt_train, rep(1/(2*k+1),2*k+1), sides = 2,circular = T)
  mod4 <-  filter(yt_train,w5,"convolution",2,circular=T)
  ythatmod=ts(fitted(mod1)[1:length(yt_train)],freq=365,start=c(anio_inicio,mes_inicio))
  
  #Predicciones
  predSES=predict(mod1,n.ahead=m,prediction.interval=T,level=0.95)
  predMMSU = rep(tail(mod2, 1), 5)
  predMMSB = rep(tail(mod3, 1), 5)
  predFH = rep(tail(mod4,1), 5)
  ytPron1 = ts(predSES[,1],freq=1,start=c(anio_pron,mes_pron,27))
  ytPron2 = ts(predMMSU,freq=1,start=c(anio_pron,mes_pron,27))
  ytPron3 = ts(predMMSB,freq=1,start=c(anio_pron,mes_pron,27))
  ytPron4 = ts(predFH,freq=1,start=c(anio_pron,mes_pron,27))
  
  #Graficas de ajuste 
  win.graph()
  par(mfrow = c(2, 2))
  plot(yt_train, ylab = "Yt", main = "SES")
  lines(ythatmod,col=2)
  legend("topleft",legend=c("Real","Ajustado"),col=c(1,2,4),lty=1,lwd=2)
  plot(yt_train, ylab = "Yt", main = "Media Movil Unilateral")
  lines(mod2,col=2)
  legend("topleft",legend=c("Real","Ajustado"),col=c(1,2,4),lty=1,lwd=2)
  plot(yt_train, ylab = "Yt", main = "Media Movil Bilateral")
  lines(mod3,col=2)
  legend("topleft",legend=c("Real","Ajustado"),col=c(1,2,4),lty=1,lwd=2)
  plot(yt_train, ylab = "Yt", main = "Filtro Henderson")
  lines(mod4,col=2)
  legend("topleft",legend=c("Real","Ajustado"),col=c(1,2,4),lty=1,lwd=2)
  mtext(paste("Ajuste de modelos Locales: Serie",nombre), side = 3, line = -1, outer = TRUE,font=2)
  
  #Grafica de Pronosticos
  win.graph()
  plot(family="A",yt_test,ylim=c(min(yt_test,ytPron1,ytPron2,ytPron3,ytPron4),
                                 max(yt_test,ytPron1,ytPron2,ytPron3,ytPron4)+200),type="b",pch=19,lwd=2,xaxt="n",ylab = "Yt", main = "Valores Predichos vs Valores Reales")
  axis(1,at=time(yt_test),labels=c("2023-12-27","2023-12-28","2023-12-29","2023-12-30","2023-12-31"))
  lines(ytPron1,lty=2,col="blue",type="b",pch=2,lwd=2)
  lines(ytPron2,lty=3,col="orange",type="b",pch=3,lwd=2)
  lines(ytPron3,lty=4,col="darkgreen",type="b",pch=4,lwd=2)
  lines(ytPron4,lty=5,col="darkred",type="b",pch=5,lwd=2)
  legend("topleft",legend=c("Real","SES", "Media Movil Unilateral", "Media Movil Bilateral","Filtro de Henderson"),
         lty=c(1:5),pch=c(19,2:5),col=c("black","blue","orange","darkgreen","darkred"),lwd=2,cex=0.8)
  print('Metricas en Entrenamiento')
  #Tabulando medidas de pronosticos, RMSE, MAE, MAPE, MSE PRONOSTICO
  precision.puntualesT=rbind(accuracy(ythatmod,yt_train), accuracy(mod2,yt_train), accuracy(mod3,yt_train), accuracy(mod4,yt_train))[,c(2,3,5)]
  tabla.precisionT=cbind(precision.puntualesT)
  rownames(tabla.precisionT)=c("SES","MMSU","MMSB","FH")
  mseT <- c(MSE(ythatmod,yt_train),MSE(mod2,yt_train),MSE(mod3,yt_train),MSE(mod4,yt_train))
  tabla.precisionT <- cbind(tabla.precisionT, MSE = mseT)
  print(tabla.precisionT)
  
  print('Metricas en Validación')
  #Tabulando medidas de pronosticos, RMSE, MAE, MAPE, MSE VALIDACIÓN
  precision.puntuales=rbind(accuracy(ytPron1,yt_test), accuracy(ytPron2,yt_test), accuracy(ytPron3,yt_test), accuracy(ytPron4,yt_test))[,c(2,3,5)]
  tabla2.precision=cbind(precision.puntuales)
  rownames(tabla2.precision)=c("SES","MMSU","MMSB","FH")
  mse <- c(MSE(ytPron1,yt_test),MSE(ytPron2,yt_test),MSE(ytPron3,yt_test),MSE(ytPron4,yt_test))
  tabla2.precision <- cbind(tabla2.precision, MSE = mse)
  print(tabla2.precision)
  

}

generar_grafica_componentes<-function(serie){
  win.graph()
  descom= decompose(serie)
  plot(descom)
}

#Obtenemos las series
serie1 = obtener_serie("ALBAN",2022,1)
serie2= obtener_serie("BETANIA",2022,1)
serie3 = obtener_serie("CALIMA",2022,1)
serie4= obtener_serie("CHIVOR",2022,1)
serie5= obtener_serie("EL QUIMBO",2022,1)
serie6= obtener_serie("GUATRON",2022,1)
serie7= obtener_serie("GUAVIO",2022,1)
serie8= obtener_serie("JAGUAS",2022,1)
serie9= obtener_serie("LA TASAJERA",2022,1)
serie10= obtener_serie("MIEL I",2022,1)
serie11= obtener_serie("PAGUA",2022,1)
serie12= obtener_serie("PLAYAS",2022,1)
serie13= obtener_serie("PORCE II",2022,1)
serie14= obtener_serie("PORCE III",2022,1)
serie15= obtener_serie("PRADO",2022,1)
serie16= obtener_serie("SALVAJINA",2022,1)
serie17= obtener_serie("SAN CARLOS",2022,1)
serie18= obtener_serie("SOGAMOSO",2022,1)


#Generamos:
 #1. Graficas de ajuste de los modelos de suavizamiento y medias móviles
 #2. Graficas de pronóstico con los últimos 5 días
 #3. Medidas de pronostico en entrenamiento y validación
aplicar_filtros_suavizador(serie1,"Alban",2022,1,2023,12)
aplicar_filtros_suavizador(serie2,"BETANIA",2022,1,2023,12)
aplicar_filtros_suavizador(serie3,"CALIMA",2022,1,2023,12)
aplicar_filtros_suavizador(serie4,"CHIVOR",2022,1,2023,12)
aplicar_filtros_suavizador(serie5,"EL QUIMBO",2022,1,2023,12)
aplicar_filtros_suavizador(serie6,"GUATRON",2022,1,2023,12)
aplicar_filtros_suavizador(serie7,"GUAVIO",2022,1,2023,12)
aplicar_filtros_suavizador(serie8,"JAGUAS",2022,1,2023,12)
aplicar_filtros_suavizador(serie9,"LA TASAJERA",2022,1,2023,12)
aplicar_filtros_suavizador(serie10,"MIEL I",2022,1,2023,12)
aplicar_filtros_suavizador(serie11,"PAGUA",2022,1,2023,12)
aplicar_filtros_suavizador(serie12,"PLAYAS",2022,1,2023,12)
aplicar_filtros_suavizador(serie13,"PORCE II",2022,1,2023,12)
aplicar_filtros_suavizador(serie14,"PORCE III",2022,1,2023,12)
aplicar_filtros_suavizador(serie15,"PRADO",2022,1,2023,12)
aplicar_filtros_suavizador(serie16,"SALVAJINA",2022,1,2023,12)
aplicar_filtros_suavizador(serie17,"SAN CARLOS",2022,1,2023,12)
aplicar_filtros_suavizador(serie18,"SOGAMOSO",2022,1,2023,12)



generar_grafica_componentes(serie1)
generar_grafica_componentes(serie2)
generar_grafica_componentes(serie3)
generar_grafica_componentes(serie4)
generar_grafica_componentes(serie5)
generar_grafica_componentes(serie6)
generar_grafica_componentes(serie7)
generar_grafica_componentes(serie8)
generar_grafica_componentes(serie9)
generar_grafica_componentes(serie10)
generar_grafica_componentes(serie11)
generar_grafica_componentes(serie12)
generar_grafica_componentes(serie13)
generar_grafica_componentes(serie14)
generar_grafica_componentes(serie15)
generar_grafica_componentes(serie16)




