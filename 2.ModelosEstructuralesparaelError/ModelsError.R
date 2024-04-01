#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(forecast)
library(lmtest);
library(TSA)
library("showtext")
windowsFonts(A = windowsFont("Times New Roman"))
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")

#Leer los datos
Datos16=read.table(file.choose(),header=T,sep=";",skip=9,dec=",",colClasses=c(rep("NULL",8),"numeric",rep("NULL",10)))
serie=ts(Datos16,freq=12,start=c(2003,1))

#DEFINIENDO VARIABLES Y CREACION DE LA MATRIZ
n=length(serie)-12
t=1:n
t2=t^2
t3=t^3
t4=t^4
t5=t^5
t6=t^6
yt=ts(serie[t],freq=12,start=c(2003,1))
sen1=sin(pi*t/3)
cos1=cos(pi*t/3)
sen2=sin(pi*t/2)
cos2=cos(pi*t/2)
sen3=sin(2*pi*t/3)
cos3=cos(2*pi*t/3)
sen4=sin(5*pi*t/6)
cos4=cos(5*pi*t/6)
X=data.frame(t,t2,t3,t4,t5,t6,sen1,cos1,sen2,cos2,sen3,cos3,sen4,cos4)

#Definiendo variables para pron?sticos
tnuevo=(n+1):length(serie) #?ndice de tiempo en los pron?sticos
t2nuevo=tnuevo^2
t3nuevo=tnuevo^3
t4nuevo=tnuevo^4
t5nuevo=tnuevo^5
t6nuevo=tnuevo^6
ytnuevo=ts(serie[tnuevo],freq=12,start=c(2018,7))
sen1n=sin(pi*tnuevo/3)
cos1n=cos(pi*tnuevo/3)
sen2n=sin(pi*tnuevo/2)
cos2n=cos(pi*tnuevo/2)
sen3n=sin(2*pi*tnuevo/3)
cos3n=cos(2*pi*tnuevo/3)
sen4n=sin(5*pi*tnuevo/6)
cos4n=cos(5*pi*tnuevo/6)
#matriz de predictores en el pron?stico
Xnuevo=data.frame(t=tnuevo,t2=t2nuevo,t3=t3nuevo,t4=t4nuevo,t5=t5nuevo,t6=t6nuevo,sen1=sen1n,cos1=cos1n,sen2=sen2n,cos2=cos2n,sen3=sen3n,cos3=cos3n,
                  sen4=sen4n,cos4=cos4n)

modeloglobal=lm(yt~.,data=X)
summary(modeloglobal)
modeloglobal$SSE

#EACF
eacf(serie,ar.max=36,ma.max=36)


#-----MODELO 1:AR(13)------------------------------------------------------
modelo1=Arima(yt,order=c(13,0,0),xreg=as.matrix(X),method="ML") 
k1=length(coef(modelo1)[coef(modelo1)!=0]);k1 
dfmodelo1=n-k1
coeftest(modelo1,df=dfmodelo1) 
summary(modelo1)

#-----MODELO 2:ARMA(2,10)------------------------------------------------------
#MODELO 2 ARMA(2,10)
modelo2=Arima(yt,order=c(2,0,10),fixed=c(rep(NA,12),coef(modeloglobal)),xreg=as.matrix(X),method="ML")
k2=length(coef(modelo2)[coef(modelo2)!=0]);k2 
dfmodelo2=n-k2
coeftest(modelo2,df=dfmodelo2) 
summary(modelo2)

#-----MODELO 2:ARMA (2,0) (1,0) [12]------------------------------------------------------
modelo3=Arima(yt,order=c(2,0,0),seasonal=list(order=c(1,0,0)),xreg=as.matrix(X),fixed=c(rep(NA,3),coef(modeloglobal)),method="ML") 
k3=length(coef(modelo3)[coef(modelo3)!=0]);k3 
dfmodelo3=n-k3
coeftest(modelo3,df=dfmodelo3) 
summary(modelo3)

#-----MODELO 2:AR(1,2,4,13)------------------------------------------------------
#ARMABUSSETS 18X18 m?todo ml renglon 2 
modelo4=Arima(yt,order=c(13,0,0),xreg=as.matrix(X),fixed= c(NA,NA,0,NA,rep(0,8),NA,rep(NA,15)),method="ML")
k4=length(coef(modelo4)[coef(modelo4)!=0]);k4 
dfmodelo4=n-k4
coeftest(modelo4,df=dfmodelo4) 
summary(modelo4)
#------------------------------------------------------------------
#------------------------------------------------------------------
#INDICADORES-PRUEBAS
#AIC-BIC
Criteriosmodelo1=exp.crit.inf.resid(residuales=residuals(modelo1),n.par=k1); Criteriosmodelo1
Criteriosmodelo2=exp.crit.inf.resid(residuales=residuals(modelo2),n.par=k2); Criteriosmodelo2 
Criteriosmodelo3=exp.crit.inf.resid(residuales=residuals(modelo3),n.par=k3); Criteriosmodelo3
Criteriosmodelo4=exp.crit.inf.resid(residuales=residuals(modelo4),n.par=k4); Criteriosmodelo4
#Ljun-Box
BP.LB.test(residuals(modelo1),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo2),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo3),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo4),maxlag=36,type="Ljung")
#Normalidad
shapiro.test(residuals(modelo1))
shapiro.test(residuals(modelo2))
shapiro.test(residuals(modelo3))
shapiro.test(residuals(modelo4))
#Pronosticos
predmodelo1=ts(as.data.frame(forecast(modelo1,xreg=as.matrix(Xnuevo),level=95)),freq=12, start=c(2018,7));predmodelo1
ytpronmodelo1=predmodelo1[,1]
predmodelo2=ts(as.data.frame(forecast(modelo2,xreg=as.matrix(Xnuevo),level=95)),freq=12, start=c(2018,7));predmodelo2
ytpronmodelo2=predmodelo2[,1]
predmodelo3=ts(as.data.frame(forecast(modelo3,xreg=as.matrix(Xnuevo),level=95)),freq=12, start=c(2018,7));predmodelo3
ytpronmodelo3=predmodelo3[,1]
predmodelo4=ts(as.data.frame(forecast(modelo4,xreg=as.matrix(Xnuevo),level=95)),freq=12, start=c(2018,7));predmodelo4
ytpronmodelo4=predmodelo4[,1]
#Medidas de pronostico
accuracy(ytpronmodelo1,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmodelo1[,2],LSP=predmodelo1[,3]) 
accuracy(ytpronmodelo2,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmodelo2[,2],LSP=predmodelo2[,3]) 
accuracy(ytpronmodelo3,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmodelo3[,2],LSP=predmodelo3[,3])
accuracy(ytpronmodelo4,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmodelo4[,2],LSP=predmodelo4[,3]) 


#------------------------------------------------------------------
#------------------------------------------------------------------
#GRAFICAS
#Graficas de ajuste:
#modelo 1
win.graph(width=5,height=5,pointsize=15)
ythat1=modelo1$fitted #Este objeto ya tiene fechas
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat1,col=2)
legend("topleft",legend=c("Datos","Ajuste Modelo 1"),lty=1,col=1:2,lwd=2, cex=0.8)
#modelo 2
win.graph(width=5,height=5,pointsize=8)
ythat2=modelo2$fitted #Este objeto ya tiene fechas
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat2,col=2)
legend("topleft",legend=c("Datos","Ajuste Modelo 2"),lty=1,col=1:2,lwd=2, cex=0.8)
#modelo 3
win.graph(width=5,height=5,pointsize=8)
ythat3=modelo3$fitted #este objeto ya queda con las fechas de la serie
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat3,col=2)
legend("topleft",legend=c("Datos","Ajuste Modelo 3"),lty=1,col=1:2,lwd=2, cex=0.8)
#modelo 4
win.graph(width=5,height=5,pointsize=8)
ythat4=modelo4$fitted #este objeto ya queda con las fechas de la serie
plot(serie,ylab='serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat4,col=2)
legend("topleft",legend=c("datos","Ajuste modelo 2"),lty=1,col=1:2,lwd=2, cex=0.8)
#Graficas de ACF
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo1)),ci.type="ma",lag.max=36,main="ACF Modelo 1",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo2)),ci.type="ma",lag.max=36,main="ACF Modelo 2",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo3)),ci.type="ma",lag.max=36,main="ACF Modelo 3",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo4)),ci.type="ma",lag.max=36,main="ACF Modelo 4",ci.col=2,family = "A")
#Graficas de PACF 
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo1)),lag.max=36,main="PACF Modelo 1",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo2)),lag.max=36,main="PACF Modelo 2",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo3)),lag.max=36,main="PACF Modelo 3",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo4)),lag.max=36,main="PACF Modelo 4",ci.col=2,family = "A")
#Graficas de Normalidad
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo1),main="Gr?fico de normalidad residuos de ajuste Modelo 1",family="A")
qqline(residuals(modelo1),col=2) 
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo2),main="Gr?fico de normalidad residuos de ajuste Modelo 2",family="A")
qqline(residuals(modelo2),col=2) 
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo3),main="Gr?fico de normalidad residuos de ajuste Modelo 3",family = "A")
qqline(residuals(modelo3),col=2)
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo4),main="Gr?fico de normalidad residuos de ajuste Modelo 4",family = "A")
qqline(residuals(modelo4),col=2) 
#Graficas de Residuos vs Tiempo
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo1),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo2),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo2$sigma2),2*sqrt(modelo2$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo3),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo4),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2)
#Graficas de Residuos vs valores ajustados 
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo1$fitted),residuals(modelo1),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2) 
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo2$fitted),residuals(modelo2),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo2$sigma2),2*sqrt(modelo2$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo3$fitted),residuals(modelo3),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2) 
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo4$fitted),residuals(modelo4),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2)
#Comparaci?on gr?afica de los pron?osticos y la serie real en periodos ex-post
win.graph(width=5,height=5,pointsize=10)
p=cbind(ytnuevo,ytpronmodelo1,ytpronmodelo2,ytpronmodelo3,ytpronmodelo4)
p
plot(main="Valores reales y pron?sticos ?ltimos 12 meses",family="A",p,plot.type="single",type="b",pch=c(19,1:4),col=c("black","blue","orange","darkgreen","darkred"),lty=1,lwd=2,ylab="Ytf",xaxt="n")
axis(1,at=time(p),labels=c("Jul-18","Ago-18","Sep-18","Oct-18","Nov-18","Dic-18","Ene-19","Feb-19","Mar-19","Abr-19","May-19","Jun-19"),family="A")
legend("topleft",legend=c("Serie Original","Pron?stico Modelo 1","Pron?stico Modelo 2","Pron?stico Modelo 3","Pron?stico Modelo 4"),
       pch=c(19,1:4),lty=c(1:4),col=c("black","blue","orange","darkgreen","darkred"),lwd=2,cex=0.8)

#Modelo con ajuste SEHW
ModeloSWHE=SuavizamientoEstacional(yt,seasonal="additive",h=12,beta=1e-5)
#Serie y su ajuste
win.graph(width=5,height=5,pointsize=8)
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(fitted(ModeloSWHE),col=2)
legend("topleft",legend=c("Datos","Ajuste Modelo 4"),lty=1,col=1:2,lwd=2, cex=0.8)
p3=(12-1)+2 #Aprox. del n?mero de par?metros del suavizamiento
Criterios3=exp.crit.inf.resid(residuales=residuals(ModeloSWHE),n.par=p3);Criterios3 
pronos3= ModeloSWHE$forecast;pronos3 
ytpron3=pronos3[,1] #s?lo los pron?sticos puntuales del suavizamiento 
accuracy(ytpron3,ytnuevo) #Precisi?n pron?sticos puntuales
#Precisi?n pron?sticos por I.P.
amplitud.cobertura(real=ytnuevo,LIP=pronos3[,2],LSP=pronos3[,3])
MSE3=ModeloSWHE$MSE;MSE3 
#residuales
win.graph(width=5,height=5,pointsize=10)
plot(ylab="Modelolocal",family = "A",residuals(ModeloSWHE),ylim=c(min(-2*sqrt(MSE3),residuals(ModeloSWHE)),max(2*sqrt(MSE3),residuals(ModeloSWHE))))
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),col=2)
legend("topleft",legend="Modelo 3",family = "A")
win.graph(width=5,height=5,pointsize=10)
plot(ylab="Modelolocal",family = "A",as.numeric(fitted(ModeloSWHE)),residuals(ModeloSWHE), ylim=c(min(-2*sqrt(MSE3),residuals(ModeloSWHE)),max(2*sqrt(MSE3),residuals(ModeloSWHE))))
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),col=2)
legend("topleft",legend="Modelo 3",family = "A") 
#ACF
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(ModeloSWHE)),ci.type="ma",lag.max=36,main="ACF Modelo local",ci.col=2,family="A")
#PACF
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(ModeloSWHE)),lag.max=36,main="PACF Modelo local",ci.col=2,family="A")
#Normalidad
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(ModeloSWHE),main="Gr?fico de normalidad residuos de ajuste Modelo local",family="A")
qqline(residuals(ModeloSWHE),col=2) 
shapiro.test(residuals(ModeloSWHE))
#MODELO 4
ajusteDLL1b=DBP.LB.test(residuals(modelo2),maxlag=36,type="Ljung")
BP.LB.test(residuals(ModeloSWHE),maxlag=36,type="Ljung")