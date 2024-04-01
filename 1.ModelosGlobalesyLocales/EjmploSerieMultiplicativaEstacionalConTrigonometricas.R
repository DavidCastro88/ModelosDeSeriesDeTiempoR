#LIBRERIAS Y PAQUETES ----------------------------------------------------------

rm(list=ls(all=TRUE))
library(TSA)
library(forecast)
library(fANCOVA)
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")

#Leer anex-EMMET-TotalNacional-oct2023-Elaboracion de Bebidas.csv, columna 5: Produccion nominal
Datos7=read.table(file.choose(),header=T,sep=";",skip=14,dec=",",colClasses=c(rep("NULL",4),"numeric",rep("NULL",6)))
Datos7=ts(Datos7,freq=12,start=c(2001,1))
plot(Datos7)

#DEFINIENDO VARIABLES Y CREACION DE LA MATRIZ NO MOVER -------------------------
m=12 
n=length(Datos7)-m

#PARA EL AJUSTE

t=1:n
t2=t^2
t3=t^3
t4=t^4
yt=ts(Datos7[t],freq=m,start=c(2001,1))

#Funciones trigonometricas en las frecuencias j/12, j=1,2,3,4,5
sen1=sin(pi*t/6)
cos1=cos(pi*t/6)
sen2=sin(pi*t/3)
cos2=cos(pi*t/3)
sen3=sin(pi*t/2)
cos3=cos(pi*t/2)
sen4=sin(2*pi*t/3)
cos4=cos(2*pi*t/3)
cos5=sin(5*pi*t/6)

#Matriz de diseño en el ajuste
X1=data.frame(t,t2,t3,t4,sen1,cos1,sen2,cos2,sen3,cos3,sen4,cos4,cos5)
head(X1) 

#PARA LOS PRONOSTICOS

tnuevo=(n+1):length(Datos7)
t2nuevo=tnuevo^2
t3nuevo=tnuevo^3
t4nuevo=tnuevo^4
ytnuevo=ts(Datos7[tnuevo],freq=m,start=c(2022,11))

#Funciones trigonometricas en los periodos de pronostico
sen1n=sin(pi*tnuevo/6)
cos1n=cos(pi*tnuevo/6)
sen2n=sin(pi*tnuevo/3)
cos2n=cos(pi*tnuevo/3)
sen3n=sin(pi*tnuevo/2)
cos3n=cos(pi*tnuevo/2)
sen4n=sin(2*pi*tnuevo/3)
cos4n=cos(2*pi*tnuevo/3)
cos5n=cos(5*pi*tnuevo/6) 

#Matriz de diseño en los pronosticos
X1nuevo=data.frame(t=tnuevo,t2=t2nuevo,t3=t3nuevo,t4=t4nuevo,sen1=sen1n,cos1=cos1n,sen2=sen2n,cos2=cos2n,sen3=sen3n,cos3=cos3n,sen4=sen4n,cos4=cos4n,cos5=cos5n)
head(X1nuevo) 

#--ANALISIS DESCRIPTIVO --------------------------------------------------------

#Graficando la serie 
win.graph()
plot(Datos7, ylab="Datos7",main="Gráfica de la serie")

#Realizando descompose
win.graph()
descom.multiplicativo= decompose(Datos7,type = "multiplicative")
plot(descom.multiplicativo)

#Graficando la serie en escala logaritmica 
win.graph()
plot(log(Datos7), ylab="Datos7",main="Gráfica de la serie en escala logarítmica")

#Graficando la tendencia 
win.graph()
Tt.log=decompose(log(Datos7))$trend
plot(Tt.log,ylim=c(min(log(Datos7)),max(log(Datos7))),main="Gráfica de la tendencia")

#Graficando boxplot 
win.graph()
boxplot(log(Datos7)~cycle(Datos7),main="Gráfico de boxplots")

#periodograma 
win.graph()
x=diff(log(Datos7))  
plot(x,ylab=expression(log(Y[t])-log(Y[t-1])));abline(h=mean(x))
periodogram(x,main="Periodograma");abline(v=c(1:6)/12,col=2,lty=2) 

#-----MODELO 1:COMPLETAMENTE MULTIPLICATIVO------------------------------------------------------

mod1=lm(log(yt)~.,data=X1)
summary(mod1)

#Calculo valores ajustados del modelo 1
ythatmod1=ts(exp(fitted(mod1))*exp(summary(mod1)$sigma^2/2),freq=m,start=start(yt))

#Calculo de los criterios AIC y BIC en modelo 1
nparmod1=length(coef(mod1)[coef(mod1)!=0]);nparmod1 
res.orig.mod1=yt-ythatmod1
Criterios1= exp.crit.inf.resid(residuales=res.orig.mod1,n.par=nparmod1);Criterios1

#Pronosticos del modelo 1 en la escala original
pronos1=exp(predict(mod1,newdata=X1nuevo,interval="prediction",level=0.95))*exp(summary(mod1)$sigma^2/2)
pronos1=ts(pronos1,freq=m,start=start(ytnuevo))
pronos1
ytpron1=pronos1[,1]

#precision pronosticos puntuales modelo 1
accuracy(ytpron1,ytnuevo)

#precision pronosticos por I.P modelo 1
amplcobmod1=amplitud.cobertura(real=ytnuevo,LIP=pronos1[,2],LSP=pronos1[,3]);amplcobmod1

#-----MODELO 2:PARCIALMENTE MULTIPLICATIVO------------------------------------------------------

#beta: parametros para la tendencia 
#delta: parametros para la estacionalidad (INDICADORAS)
#"alfa1","gamma1","alfa2","gamma2","alfa3","gamma3","alfa4","gamma4","alfa5","gamma5","gamma6"(TRIGONOMETRICAS)

parametros.mod2=c(paste0("beta",0:4),"alfa1","gamma1","alfa2","gamma2","alfa3","gamma3","alfa4","gamma4","gamma5")
parametros.mod2
mod2=regexponencial(respuesta=yt,data=X1,names.param=parametros.mod2)
summary(mod2)

#Calculo valores ajustados del modelo 2
ythatmod2=ts(fitted(mod2),freq=m,start=start(yt))

#Calculo de los criterios AIC y BIC en modelo 2
nparmod2=length(coef(mod2)[coef(mod2)!=0]);nparmod2          
Criterios2=exp.crit.inf.resid(residuales=residuals(mod2),n.par=nparmod2); Criterios2

#Pronosticos del modelo 2 en la escala original, solo son de tipo puntual por ser modelo no lineal
pronos2=predict(mod2,newdata=X1nuevo,interval="prediction",level=0.95)
ytpron2=ts(pronos2,freq=m,start=start(ytnuevo))
ytpron2 

#precision pronosticos puntuales modelo 2
accuracy(ytpron2,ytnuevo) 

#-----MODELO 3:SEHW ------------------------------------------------------------

mod3=SuavizamientoEstacional(yt,seasonal="multiplicative",h=m,optim.start = c(alpha = 0.3, beta = 0.01, gamma = 0.1))
str(mod3) 

#Calculo de AIC y BIC 
s=m
p3=(s-1)+2
Criterios3=exp.crit.inf.resid(residuales=residuals(mod3),n.par=p3);Criterios3
MSE3=mod3$MSE 

#Pronosticos del modelo 3 
pronos3=mod3$forecast
pronos3
ytpron3=pronos3[,1] #solo los pronosticos puntuales del suavizamiento

#precision pronosticos puntuales modelo 3
accuracy(ytpron3,ytnuevo) 

#Precision pronosticos por I.P Modelo 3
amplcobmod3=amplitud.cobertura(real=ytnuevo,LIP=pronos3[,2],LSP=pronos3[,3]);amplcobmod3

#-----MODELO 4:DLC(aicc) -------------------------------------------------------

mod4=Descomp.Loess(serie.ajuste=yt,h=m,tipo.descomp="multiplicative",grado=2,criterio="aicc")
str(mod4) 

#Calculo AIC y BIC 
Criterios4=exp.crit.inf.resid(residuales=residuals(mod4),n.par=mod4$p);Criterios4
MSE4=mod4$MSE 

#Pronosticos del modelo 4 
pronos4=mod4$tablapron

#Precision pronosticos puntuales
ytpron4=mod4$ytpron
accuracy(ytpron4,ytnuevo)

#Graficando St estimada por el filtro de descomposicion
win.graph()
plot(mod4$St,ylab=expression(hat(S)[t]), main="Gráfica de la estimación de la estacionalidad") 

#--TABLAS DE LAS MEDIDAS DE AJUSTE Y PRONOSTICOS DE LOS 4 MODELOS --------------

Modelo1=summary(mod1)$coefficients
Modelo2=summary(mod2)$coefficients

tabla.parametros.globales=cbind(Modelo1,Modelo2)
colnames(tabla.parametros.globales)=c("Modelo 1","Modelo 2")
tabla.parametros.globales

#Tabulando medidas de ajuste
tabla1.criterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tabla1.criterios)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla1.criterios

#tabla resumen de pronosticos 
tabla.pronosticos=cbind(pronos1,ytpron2,pronos3,ytpron4)
colnames(tabla.pronosticos)=c("Modelo 1","Modelo 2","Modelo 3","Modelo4")
tabla.pronosticos

#Tabulando medidas de pronosticos
precision.puntuales=rbind(accuracy(ytpron1,ytnuevo), accuracy(ytpron2,ytnuevo), accuracy(ytpron3,ytnuevo), accuracy(ytpron4,ytnuevo))[,c(2,3,5)]
precision.intervalos=rbind(amplcobmod1,c(NA,NA),amplcobmod3,c(NA,NA))
tabla2.precision=cbind(precision.puntuales,precision.intervalos)
rownames(tabla2.precision)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla2.precision

#--GRAFICAS VALORES AJUSTADOS PARA LOS 4 MODELOS -------------------------------

#MODELO 1
win.graph()
plot(Datos7,ylab="Datos7",main="Modelo 1")
lines(ythatmod1,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 1"),lty=1,col=c(1,2))

#MODELO 2
win.graph()
plot(Datos7,ylab="Datos7",main="Modelo 2")
lines(ythatmod2,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 2"),lty=1,col=c(1,2)) 

#MODELO 3
win.graph()
plot(Datos7,ylab="Datos7",main="Modelo 3")
lines(fitted(mod3),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 3"),lty=1,col=c(1,2))

#MODELO 4**
win.graph()
plot(Datos7,ylab="Datos7",main="Modelo 4")
lines(fitted(mod4),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 4"),lty=1,col=c(1,2))

#**Graficando la serie desestacionalizada y su ajuste loess para la estacionalidad
win.graph()
plot(mod4$ytd,ylab=NA,main="Gráfica del ajuste loess sobre la serie desestacionalizada")
lines(mod4$Tt,col=2,lwd=2)
legend("topleft",legend=c("Serie ajustada estacionalmente","Tendencia LOESS"),col=c(1,2),lty=1)

#--GRAFICAS RESIDUOS DE AJUSTE PARA LOS 4 MODELOS ------------------------------

#MODELO 1
win.graph()
plot.ts(residuals(mod1),ylab="Residuales",main="Residuales Modelo 1")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 1")
#Residuos vs. ajustados
win.graph()
plot(fitted(mod1),residuals(mod1),ylab="Residuos vs. ajustados",main="Residuos vs. ajustados Modelo 1")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 1")

#MODELO 2
win.graph()
plot.ts(residuals(mod2),ylab="Residuales",main="Residuales Modelo 2")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 2")
#Residuos vs. ajustados
win.graph()
plot(fitted(mod2),residuals(mod2),ylab="Residuos vs. ajustados",main="Residuos vs. ajustados Modelo 2")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 2")

#MODELO 3
win.graph()
plot(residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuales",main="Residuales Modelo 3")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
legend("topleft",legend="Modelo 3")
#Residuos vs. ajustados
win.graph()
plot(as.numeric(fitted(mod3)),residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuos vs. ajustados",main="Residuos vs. ajustados Modelo 3")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
legend("topleft",legend="Modelo 3")

#MODELO 4 
win.graph()
plot(residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuales",main="Residuales Modelo 4")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)
legend("topleft",legend="Modelo 4")
#Residuos vs. ajustados
win.graph()
plot(as.numeric(fitted(mod4)),residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuos vs. ajustados",main="Residuos vs. ajustados Modelo 4")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)
legend("topleft",legend="Modelo 4")

#--GRAFICA COMPARATIVA DE LOS PRONOSTICOS DE LOS 4 MODELOS ---------------------

win.graph()
plot(ytnuevo,type="b",ylab="Datos7",col=1,pch=19,ylim=c(min(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4),max(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)+25),lwd=2,xaxt="n", main="Gráfica comparativa de los pronósticos de los 4 modelos")
axis(1,at=time(ytnuevo),labels=c("22-XI","22-XII","23-I","23-II","23-III","23-IV","23-V","23-VI","23-VII","23-VIII","23-IX","23-X"),cex.axis=0.7)
lines(ytpron1,col=2,pch=1,type="b",lwd=2)
lines(ytpron2,col=3,pch=2,type="b",lwd=2)
lines(ytpron3,col=4,pch=3,type="b",lwd=2)
lines(ytpron4,col=5,pch=4,type="b",lwd=2)
legend("topleft",legend=c("Real","Modelo 1","Modelo 2","Modelo 3","Modelo 4"),pch=c(19,1:4),col=c(1:5),lwd=2)


#RESUMEN PROGRAMACION-----------------------------------------------------------

#Grafica calidad de ajuste
win.graph()
layout(cbind(c(1,3),c(2,4)))
plot(Datos7, ylab="Datos7",main="Modelo 1");lines(ythatmod1,col=2,lwd=1);legend("topleft",legend=c("Modelo original","Ajuste Modelo 1"),lty=1,col=1:2)
plot(Datos7, ylab="Datos7",main="Modelo 2");lines(ythatmod2,col=2,lwd=1);legend("topleft",legend=c("Modelo original","Ajuste Modelo 2"),lty=1,col=1:2)
plot(Datos7, ylab="Datos7",main="Modelo 3");lines(fitted(mod3),col=2,lwd=1);legend("topleft",legend=c("Modelo original","Ajuste Modelo 3"),lty=1,col=1:2)
plot(Datos7, ylab="Datos7",main="Modelo 4");lines(fitted(mod4),col=2,lwd=1);legend("topleft",legend=c("Modelo original","Ajuste Modelo 4"),lty=1,col=1:2)

#Analisis de residuales
win.graph()
layout(cbind(c(1,3),c(2,4)))
plot.ts(residuals(mod1),main="Modelo 1");abline(h=0)
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=2, col=2)
plot.ts(residuals(mod2),main="Modelo 2");abline(h=0)
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=2, col=2)
plot(residuals(mod3),main="Modelo 3");abline(h=0)
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=2,col=2)
plot(residuals(mod4),main="Modelo 4");abline(h=0)
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=2,col=2)

win.graph()
layout(cbind(c(1,3),c(2,4)))
plot(fitted(mod1),residuals(mod1),main="Modelo 1");abline(h=0)
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=2, col=2)
plot(fitted(mod2),residuals(mod2),main="Modelo 2");abline(h=0)
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=2, col=2)
plot(as.numeric(mod3$fitted),residuals(mod3),main="Modelo 3");abline(h=0)
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=2, col=2)
plot(as.numeric(mod4$fitted),residuals(mod4),main="Modelo 4");abline(h=0)
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=2, col=2)

#Exportacion tablas al directorio de trabajo------------------------------------

write.csv2(tabla.parametros.globales,file="tablamod1ymod2trabajo1.csv",row.names = TRUE)
write.csv2(mod3$coefficients,file="tablamod3trabajo1.csv",row.names = TRUE)
write.csv2(mod4$deltasi,file="tablamod4trabajo1.csv",row.names = TRUE)
write.csv2(tabla1.criterios,file="tablacriteriostrabajo1.csv",row.names = TRUE)
write.csv2(tabla.pronosticos,file="tablapronosticostrabajo1.csv",row.names = TRUE)
write.csv2(pronos4,file="tablapronosticosmodelo4trabajo1.csv",row.names = TRUE)
write.csv2(tabla2.precision,file="tablamedidasprecisiondepronosticostrabajo1.csv",row.names = TRUE)


