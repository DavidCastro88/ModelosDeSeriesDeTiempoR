#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(TSA)
library(forecast)
library(fANCOVA)
source("https://raw.githubusdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==ercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loess.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-d-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")

#Leer anex-EMMET-TotalNacional-oct2023-Fabricacion de Calzado.csv, columna 7: Ventas nominal
Datos21=read.table(file.choose(),header=T,sep=";",skip=14,dec=",",colClasses=c(rep("NULL",6),"numeric",rep("NULL",4)))
Datos21=ts(Datos21,freq=12,start=c(2001,1))
Datos21
plot(Datos21)

#DEFINIENDO VARIABLES Y CREACION DE LA MATRIZ NO MOVER -------------------------
m=12
n=length(Datos21)-m

#PARA EL AJUSTE

t=1:n
t2=t^2
t3=t^3
t4=t^4
yt=ts(Datos21[t],freq=m,start=c(2001,1))
mes=seasonaldummy(yt)

#Matriz de diseño en el ajuste
X1=data.frame(t,mes)
head(X1) 

X2=data.frame(t,t2,t3,t4,mes)
head(X2) 

#PARA LOS PRONOSTICOS

tnuevo=(n+1):length(Datos21) 
t2nuevo=tnuevo^2
t3nuevo=tnuevo^3
t4nuevo=tnuevo^4
ytnuevo=ts(Datos21[tnuevo],freq=m,start=c(2022,11)) 
mesnuevo=seasonaldummy(yt,h=m)

#Matriz de diseño en los pronosticos
X1nuevo=data.frame(t=tnuevo,mesnuevo)
head(X1nuevo) 

X2nuevo=data.frame(t=tnuevo,t2=t2nuevo,t3=t3nuevo,t4=t4nuevo,mesnuevo)
head(X2nuevo) 

#--ANALISIS DESCRIPTIVO --------------------------------------------------------

#Graficando la serie 
win.graph()
plot(Datos21, ylab="Datos21",xlab="Time")

#Realizando descompose
win.graph()
descom.aditivo= decompose(Datos21,type = "additive")
plot(descom.aditivo)

#Graficando la tendencia 
win.graph()
Tt=decompose(Datos21)$trend
plot(Tt,ylim=c(min(Datos21),max(Datos21)))

#Graficando boxplot 
win.graph()
boxplot(Datos21~cycle(Datos21),names=month.abb)

#periodograma 
win.graph()
periodogram(diff(Datos21),lwd=3)
abline(h=0)
abline(v=c(2:6)/12,col=2,lty=2)

#-----MODELO 1:POLINOMIAL-GRADO1 LINEAL------------------------------------------------------

mod1=lm((yt)~.,data=X1)
summary(mod1)

#Calculo valores ajustados del modelo 1
ythatmod1=ts(fitted(mod1),freq=m,start=start(yt))

#Calculo de los criterios AIC y BIC en modelo 1
nparmod1=length(coef(mod1)[coef(mod1)!=0]);nparmod1 
Criterios1= exp.crit.inf.resid(residuales=residuals(mod1),n.par=nparmod1);Criterios1

#Pronosticos del modelo 1 
pronos1=predict(mod1,newdata=X1nuevo,interval="prediction",level=0.95)
pronos1=ts(pronos1,freq=m,start=start(ytnuevo))
pronos1
ytpron1=pronos1[,1]

#precision pronosticos puntuales modelo 1
accuracy(ytpron1,ytnuevo)

#precision pronosticos por I.P modelo 1
amplcobmod1=amplitud.cobertura(real=ytnuevo,LIP=pronos1[,2],LSP=pronos1[,3]);amplcobmod1

#-----MODELO 2:POLINOMIAL-GRADO 4------------------------------------------------------

mod2=lm((yt)~.,data=X2)
summary(mod2)

#Calculo valores ajustados del modelo 2
ythatmod2=ts(fitted(mod2),freq=m,start=start(yt))

#Calculo de los criterios AIC y BIC en modelo 2
nparmod2=length(coef(mod2)[coef(mod2)!=0]);nparmod2 
Criterios2= exp.crit.inf.resid(residuales=residuals(mod2),n.par=nparmod2);Criterios2

#Pronosticos del modelo 2 
pronos2=predict(mod2,newdata=X2nuevo,interval="prediction",level=0.95)
pronos2=ts(pronos2,freq=m,start=start(ytnuevo))
pronos2
ytpron2=pronos2[,1]

#precision pronosticos puntuales modelo 1
accuracy(ytpron2,ytnuevo)

#precision pronosticos por I.P modelo 1
amplcobmod2=amplitud.cobertura(real=ytnuevo,LIP=pronos2[,2],LSP=pronos2[,3]);amplcobmod2


#-----MODELO 3:SEHW ------------------------------------------------------------

mod3=SuavizamientoEstacional(yt,seasonal="additive",h=m,beta=1e-5)
str(mod3) 

#Calculo de AIC y BIC 

p3=(m-1)+2
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


#-----MODELO 4:DLL(GCV) -------------------------------------------------------

mod4=Descomp.Loess(serie.ajuste=yt,h=m,tipo.descomp="additive",grado=1,criterio="gcv")
str(mod4) 

#Calculo AIC y BIC 
Criterios4=exp.crit.inf.resid(residuales=residuals(mod4),n.par=mod4$p);Criterios4
MSE4=mod4$MSE 

#Pronosticos del modelo 4 
pronos4=mod4$tablapron
pronos4

#Precision pronosticos puntuales
ytpron4=mod4$ytpron
accuracy(ytpron4,ytnuevo)

#--TABLAS DE LAS MEDIDAS DE AJUSTE Y PRONOSTICOS DE LOS 4 MODELOS --------------

Modelo1=summary(mod1)$coefficients
Modelo2=summary(mod2)$coefficients

tabla.parametros.globales=rbind(Modelo1,Modelo2)
rownames(tabla.parametros.globales)=c("Modelo 1","Modelo 2")
tabla.parametros.globales

#Tabulando medidas de ajuste
tabla1.criterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tabla1.criterios)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla1.criterios

#tabla resumen de pronosticos 
tabla.pronosticos=cbind(pronos1,pronos2,pronos3,ytpron4)
colnames(tabla.pronosticos)=c("Modelo 1","Modelo 2","Modelo 3","Modelo4")
tabla.pronosticos

#Tabulando medidas de pronosticos
precision.puntuales=rbind(accuracy(ytpron1,ytnuevo), accuracy(ytpron2,ytnuevo), accuracy(ytpron3,ytnuevo), accuracy(ytpron4,ytnuevo))[,c(2,3,5)]
precision.intervalos=rbind(amplcobmod1,amplcobmod2,amplcobmod3,c(NA,NA))
tabla2.precision=cbind(precision.puntuales,precision.intervalos)
rownames(tabla2.precision)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla2.precision


#--GRAFICAS VALORES AJUSTADOS PARA LOS 4 MODELOS -------------------------------

#MODELO 1
win.graph()
plot(Datos21, ylab="Datos21")
lines(ythatmod1,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 1"),lty=1,col=c(1,2))

#MODELO 2
win.graph()
plot(Datos21, ylab="Datos21")
lines(ythatmod2,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 2"),lty=1,col=c(1,2)) 

#MODELO 3
win.graph()
plot(Datos21, ylab="Datos21")
lines(fitted(mod3),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 3"),lty=1,col=c(1,2))

#MODELO 4
win.graph()
plot(Datos21, ylab="Datos21")
lines(fitted(mod4),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 4"),lty=1,col=c(1,2))

#**Graficando la serie desestacionalizada y su ajuste loess para la estacionalidad
win.graph()
plot(mod4$ytd,ylab=NA)
lines(mod4$Tt,col=2,lwd=2)
legend("bottomleft",legend=c("Serie ajustada estacionalmente","Tendencia LOESS"),col=c(1,2),lty=1)

#Graficando St estimada por el filtro de descomposicion
win.graph()
plot(mod4$St,ylab=expression(hat(S)[t])) 
#--GRAFICAS RESIDUOS DE AJUSTE PARA LOS 4 MODELOS ------------------------------

#MODELO 1
win.graph()
plot.ts(residuals(mod1),ylab="Residuales")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 1")
#Residuos vs. ajustados
win.graph()
plot(fitted(mod1),residuals(mod1),ylab="Residuos vs. ajustados")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 1")

#MODELO 2
win.graph()
plot.ts(residuals(mod2),ylab="Residuales")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 2")
#Residuos vs. ajustados
win.graph()
plot(fitted(mod2),residuals(mod2),ylab="Residuos vs. ajustados")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 2")

#MODELO 3
win.graph()
plot(residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuales")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
legend("topleft",legend="Modelo 3")
#Residuos vs. ajustados
win.graph()
plot(as.numeric(fitted(mod3)),residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuos vs. ajustados")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
legend("topleft",legend="Modelo 3")

#MODELO 4 
win.graph()
plot(residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuales")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)
legend("topleft",legend="Modelo 4")
#Residuos vs. ajustados
win.graph()
plot(as.numeric(fitted(mod4)),residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuos vs. ajustados")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)
legend("topleft",legend="Modelo 4")

#--GRAFICA COMPARATIVA DE LOS PRONOSTICOS DE LOS 4 MODELOS ---------------------

win.graph()
plot(ytnuevo,type="b",ylab="Datos21",col=1,pch=19,ylim=c(min(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4),max(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)+25),lwd=2,xaxt="n")
axis(1,at=time(ytnuevo),labels=c("22-XI","22-XII","23-I","23-II","23-III","23-IV","23-V","23-VI","23-VII","23-VIII","23-IX","23-X"),cex.axis=0.7)
lines(ytpron1,col="blue",pch=1,type="b",lwd=2)
lines(ytpron2,col="orange",pch=2,type="b",lwd=2)
lines(ytpron3,col="darkgreen",pch=3,type="b",lwd=2)
lines(ytpron4,col="darkred",pch=4,type="b",lwd=2)
legend("topleft",legend=c("Real","Modelo 1","Modelo 2","Modelo 3","Modelo 4"),pch=c(19,1:4),col=c("black","blue","orange","darkgreen","darkred"),lwd=2)


#Exportacion tablas al directorio de trabajo------------------------------------
write.csv2(tabla.parametros.globales,file="tablamod1ymod2trabajo1.csv",row.names = TRUE)
write.csv2(mod3$coefficients,file="tablamod3trabajo1.csv",row.names = TRUE)
write.csv2(mod4$deltasi,file="tablamod4trabajo1.csv",row.names = TRUE)
write.csv2(tabla1.criterios,file="tablacriteriostrabajo1.csv",row.names = TRUE)
write.csv2(tabla.pronosticos,file="tablapronosticostrabajo1.csv",row.names = TRUE)
write.csv2(pronos4,file="tablapronosticosmodelo4trabajo1.csv",row.names = TRUE)
write.csv2(tabla2.precision,file="tablamedidasprecisiondepronosticostrabajo1.csv",row.names = TRUE)
