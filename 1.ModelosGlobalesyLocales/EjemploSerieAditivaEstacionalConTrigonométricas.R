library(forecast)
install.packages("TSA")
library(TSA)
library(fANCOVA)
library("ggplot2")
library("showtext")

windowsFonts(A = windowsFont("Times New Roman"))
#LECTURA DE DATOS

rm(list=ls(all=TRUE))
#"EMCM Serie indices ventas nominales total nacional segun grupo mercancias-jun2019.csv" columna 9: Productos farmacéuticos y medicinales
Datos16=read.table(file.choose(),header=T,sep=";",skip=9,dec=",",colClasses=c(rep("NULL",8),"numeric",rep("NULL",10)))
serie=ts(Datos16,freq=12,start=c(2003,1))
win.graph(width=20,height=10) 
par(mfrow = c(1, 2))
plot(serie,xlab = "Time",ylab='Ventas',family = "A",cex.lab=20,cex.axis = 0.5)
plot(log(serie),xlab = "Time",ylab='log(Ventas)',family = "A")

# Comparación del patrón de la varianza del error resultante en la descomposición aditiva de la serie vs. el patrón de la varianza del error
#resultante en la descomposición multiplicativa de la serie
win.graph(width=20,height=10) 
par(mfrow = c(1, 2))
plot(decompose(serie,type='multiplicative')$random,ylab='Error',family = "A")
plot(decompose(serie,type='additive')$random,ylab='Error',family = "A")

#Componente De la tendencia
win.graph(width=20,height=10)
descom=decompose(serie,type="additive")$trend
plot(family = "A",descom,ylim=c(min(serie),max(serie)),ylab = "Tendencia")




#Determinación de la componente estacional
win.graph(width=20,height=10) 
par(mfrow = c(1, 2))
plot(serie,xlab = "Time",ylab='Yt',family = "A",cex.lab=20,cex.axis = 0.5)
boxplot(serie~cycle(serie),names=month.abb,family = "A",ylab='',xlab='')


win.graph(width=20,height=10)
estacionalidad=decompose(serie,type="additive")$seasonal
plot(estacionalidad)

#PERIODOGRAMA
x=diff(serie)
win.graph(width=10,height=5)
layout(cbind(c(1,1),c(2,2)))
ggplot(serie)
plot(x)
win.graph()
periodogram(x,lwd=3,family='A')
abline(h=0)
abline(v=c((2:6)/12),col=2,lty=2)

#AJUSTE DE LOS MODELOS PROPUESTOS CON VALIDACION CRUZADA modelos 1 y 2
rm(list=ls(all=TRUE))
source(file.choose()) 
source(file.choose()) 

#Definiendo periodos para ajuste y pronóstico
m=12 #longitud de la muestra en la validación cruzada 
n=length(serie)-m
t=1:n
t2=t^2
t3=t^3
t4=t^4
t5=t^5
t6=t^6

#Serie para ajustes
yt=ts(serie[t],freq=12,start=c(2003,1))

#Trigonométricas para ondas en las frecuencias j/12, para j=2, 3, 4, 5

sen1=sin(pi*t/3)
cos1=cos(pi*t/3)
sen2=sin(pi*t/2)
cos2=cos(pi*t/2)
sen3=sin(2*pi*t/3)
cos3=cos(2*pi*t/3)
sen4=sin(5*pi*t/6)
cos4=cos(5*pi*t/6)
#Matriz de diseño con las potencias de t en polinomio de grado 4 y trigonométricas para la estacionalidad, en el ajuste
X1=data.frame(t,t2,t3,t4,sen1,cos1,sen2,cos2,sen3,cos3,sen4,cos4)
#Matriz de diseño con las potencias de t en polinomio de grado 6 y trigonométricas para la estacionalidad, en el ajuste
X2=data.frame(t,t2,t3,t4,t5,t6,sen1,cos1,sen2,cos2,sen3,cos3,sen4,cos4)

#Valores de las variables en los pronósticos ex-post
tnuevo=(n+1):length(serie) #índice de tiempo en los pronósticos
t2nuevo=tnuevo^2
t3nuevo=tnuevo^3
t4nuevo=tnuevo^4
t5nuevo=tnuevo^5
t6nuevo=tnuevo^6

#Valores de la serie en períodos de pronósticos en este ejemplo inician en julio de 2011
ytnuevo=ts(serie[tnuevo],freq=12,start=c(2018,7))

#Funciones trigonométricas en los períodos de pronóstico
sen1n=sin(pi*tnuevo/3)
cos1n=cos(pi*tnuevo/3)
sen2n=sin(pi*tnuevo/2)
cos2n=cos(pi*tnuevo/2)
sen3n=sin(2*pi*tnuevo/3)
cos3n=cos(2*pi*tnuevo/3)
sen4n=sin(5*pi*tnuevo/6)
cos4n=cos(5*pi*tnuevo/6)

#Matriz de diseño con las potencias de t en polinomio de grado 4 y trigonométricas para la estacionalidad, en los pronósticos 
Xnuevo=data.frame(t=tnuevo,t2=t2nuevo,t3=t3nuevo,t4=t4nuevo,sen1=sen1n,cos1=cos1n,sen2=sen2n,cos2=cos2n,sen3=sen3n,cos3=cos3n,sen4=sen4n,cos4=cos4n)
modelo1=lm(yt~.,data=X1)
summary(modelo1) 

#Matriz de diseño con las potencias de t en polinomio de grado 6 y trigonométricas para la estacionalidad, en los pronósticos 
Xnuevo2=data.frame(t=tnuevo,t2=t2nuevo,t3=t3nuevo,t4=t4nuevo,t5=t5nuevo,t6=t6nuevo,sen1=sen1n,cos1=cos1n,sen2=sen2n,cos2=cos2n,sen3=sen3n,cos3=cos3n,sen4=sen4n,cos4=cos4n)
modelo2=lm(yt~.,data=X2)
summary(modelo2)

#Cálculo valores ajustados del modelo 1
ythatmodelo1=ts(fitted(modelo1),freq=12,start=start(yt))

#Calculo valores ajustados del modelo 2
ythatmodelo2=ts(fitted(modelo2),freq=12,start=start(yt))

#Cálculo de los criterios de información usando exp(C^*_n(p))
nparmodelo1=length(coef(modelo1)[coef(modelo1)!=0]); nparmodelo1 #número parámetros modelo1
Criterios1=exp.crit.inf.resid(residuales=residuals(modelo1),n.par=nparmodelo1);Criterios1

#Cálculo de los criterios de información usando exp(C^*_n(p))
nparmodelo2=length(coef(modelo2)[coef(modelo2)!=0]); nparmodelo2 #número parámetros modelo1
Criterios2=exp.crit.inf.resid(residuales=residuals(modelo2),n.par=nparmodelo2);Criterios2

#Gráfico de la serie y su ajuste con el modelo 1
plot(serie)
lines(ythatmodelo1,col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste modelo 1"),lty=1,col=c(1,2)) 

#Gráfico de la serie y su ajuste con el modelo 2
win.graph(width=5,height=5,pointsize=8)
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythatmodelo2,col=2)
legend("topleft",legend=c("Datos","Ajuste modelo 1"),lty=1,col=1:2,lwd=2, cex=0.8)
#LAS DOS JUNTAS
win.graph(width=20,height=10) 
par(mfrow = c(2, 2))
plot(serie,xlab = "Time",ylab='Yt',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythatmodelo1,col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste modelo 1            "),lty=1,col=c(1,2),cex=0.5)
plot(serie,xlab = "Time",ylab='Yt',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythatmodelo2,col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste modelo 2            "),lty=1,col=c(1,2),cex=0.5)

#pronósticos puntuales y por Intervalos del 95% Modelo 1
pronmodelo1=predict(modelo1,newdata=Xnuevo,interval="prediction",level=0.95)
pronmodelo1=ts(pronmodelo1,freq=12,start=start(ytnuevo)); pronmodelo1
ytpronmodelo1=pronmodelo1[,1]; ytpronmodelo1 #serie de los pronósticos puntuales
accuracy(ytpronmodelo1,ytnuevo) #Calculando exactitud de los pronósticos puntuales
win.graph()
plot(ytpronmodelo1)
#pronósticos puntuales y por Intervalos del 95% Modelo 1
pronmodelo2=predict(modelo2,newdata=Xnuevo2,interval="prediction",level=0.95)
pronmodelo2=ts(pronmodelo2,freq=12,start=start(ytnuevo)); pronmodelo2
ytpronmodelo2=pronmodelo2[,1]; ytpronmodelo2 #serie de los pronósticos puntuales
accuracy(ytpronmodelo2,ytnuevo) #Calculando exactitud de los pronósticos puntuales

#precisión pronósticos por I.P modelo 1
amplcobmodelo1=amplitud.cobertura(real=ytnuevo,LIP=pronmodelo1[,2],LSP=pronmodelo1[,3]);amplcobmodelo1

#precisión pronósticos por I.P modelo 2
amplcobmodelo2=amplitud.cobertura(real=ytnuevo,LIP=pronmodelo2[,2],LSP=pronmodelo2[,3]);amplcobmodelo2

#Gráficos de residuales
plot.ts(residuals(modelo1),ylim=c(min(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma),max(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma)))
abline(h=c(-2*summary(modelo1)$sigma,0,2*summary(modelo1)$sigma),col=2)

plot(fitted(modelo1),residuals(modelo1),ylim=c(min(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma),max(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma)))
abline(h=c(-2*summary(modelo1)$sigma,0,2*summary(modelo1)$sigma),col=2) 

#Gráficos de residuales
plot.ts(residuals(modelo2),ylim=c(min(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma),max(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma)))
abline(h=c(-2*summary(modelo2)$sigma,0,2*summary(modelo2)$sigma),col=2)

plot(fitted(modelo2),residuals(modelo2),ylim=c(min(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma),max(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma)))
abline(h=c(-2*summary(modelo2)$sigma,0,2*summary(modelo2)$sigma),col=2) 


#AJUSTE DE LOS MODELOS PROPUESTOS CON VALIDACION CRUZADA modelos 3 SHEW y 4 DLL(GCV)
#VARIABLES
source(file.choose())
m=12
n=length(serie)-m
t=1:n
s=12 #Longitud del periodo estacional
tnuevo=(n+1):length(serie)
#Serie para ajustes
yt=ts(serie[t],freq=12,start=c(2003,1))
#Valores de la serie en períodos de prognósticos ex-post
ytf=ts(serie[tnuevo],freq=12,start=c(2018,7)) 

#MODELO 3, SHEW 
suav=SuavizamientoEstacional(yt,seasonal="additive",h=m,beta=1e-5)
str(suav)
#Serie y su ajuste
win.graph()
plot(serie)
lines(fitted(suav),col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste H-W"),col=c(1,2),lty=1)

#Cálculo de AIC y BIC aproximados con exp(Cn*(p))
p3=(s-1)+2 #Aprox. del número de parámetros del suavizamiento
Criterios3=exp.crit.inf.resid(residuales=residuals(suav),n.par=p3);Criterios3 

MSE3=suav$MSE #MSE aproximado del ajuste total del Suavizamiento
MSE3 

win.graph()
par(mfrow = c(1, 2))
plot(residuals(suav),ylim=c(min(-2*sqrt(MSE3),residuals(suav)),max(2*sqrt(MSE3),residuals(suav))))
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),col=2)
legend("topleft",legend="Modelo 3")
plot(as.numeric(fitted(suav)),residuals(suav), ylim=c(min(-2*sqrt(MSE3),residuals(suav)),max(2*sqrt(MSE3),residuals(suav))))
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),col=2)
legend("topleft",legend="Modelo 3") 

#Pronósticos e I.P del 95% del suavizamiento
pronos3= suav$forecast
pronos3 
ytpron3=pronos3[,1] #sólo los pronósticos puntuales del suavizamiento 
accuracy(ytpron3,ytf) #Precisión pronósticos puntuales

#Precisión pronósticos por I.P.
amplitud.cobertura(real=ytf,LIP=pronos3[,2],LSP=pronos3[,3])


#MODELO 4, DLL (GCV)
ajusteDLL1b=Descomp.Loess(serie.ajuste=yt,h=m,tipo.descomp="additive",grado=1,criterio="gcv")
win.graph()
plot(ajusteDLL1b$St,ylab=expression(hat(S)[t])) #Gráfico estimación de St por el filtro de descomposición.
#Es el mismo obtenido en la primera ejecución de la función Descomp.Loess 

#Gráfico serie desestacionalizada y su ajuste loess
win.graph()
plot(ajusteDLL1b$ytd,family="A")
lines(ajusteDLL1b$Tt,col=2,lwd=2)
legend("topleft",legend=c("Serie ajustada estacionalmente","Tendencia LOESS lineal, criterio GCV         "),col=c(1,2),lty=1,cex=0.8)

#Gráfico de la serie y su ajuste final
win.graph()
plot(serie)
lines(fitted(ajusteDLL1b),col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste DLL(GCV)"),col=c(1,2),lty=1)

#Cálculo AIC y BIC aproximados versión exp(Cn*(p))
Criterios1b=exp.crit.inf.resid(residuales=residuals(ajusteDLL1b),n.par=ajusteDLL1b$p);Criterios1b


#Tabla con pronósticos de la tendencia, la estacionalidad y de la serie
ajusteDLL1b$tablapron 
ytpron4=ajusteDLL1b$ytpron #serie pron´osticos puntuales
#Precisión pronósticos puntuales
accuracy(ajusteDLL1b$ytpron,ytf) 

#Gráficos de residuos de ajuste
win.graph()
par(mfrow = c(1, 2))
plot(residuals(ajusteDLL1b),
     ylim=c(min(-2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b)),max(2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b))))
abline(h=c(-2*sqrt(ajusteDLL1b$MSE),0,2*sqrt(ajusteDLL1b$MSE)),col=2)
legend("topleft",legend="Modelo 4: DLL(GCV)")

plot(as.numeric(fitted(ajusteDLL1b)),residuals(ajusteDLL1b),
     ylim=c(min(-2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b)),max(2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b))))
abline(h=c(-2*sqrt(ajusteDLL1b$MSE),0,2*sqrt(ajusteDLL1b$MSE)),col=2)
legend("topleft",legend="Modelo 4: DLL(GCV)") 

#Comparación gráfica de efectos estacionales en los dos modelos locales

#extracción de las estimaciones en t=151 de los efectos estacionales seg´un Holt-Winters
deltasiHW=ts(coef(suav)[-c(1,2),],freq=1,start=1)
#Estimaciones de los efectos estacionales seg´un filtro de descomposici´on
deltasDescomp=ts(ajusteDLL1b$deltasi,freq=1,start=1)
#Gr´afico de los efectos estacionales estimados
win.graph()
plot(deltasiHW,lwd=3,ylim=c(min(deltasiHW,deltasDescomp),max(deltasiHW,deltasDescomp)+0.03),
     ylab="",xlab="Mes",family = "A")
lines(deltasDescomp,lty=2,lwd=3,col=2)
legend("topleft",legend=c("Efectos estacionales H-W en t=186","Efectos estacionales Filtro de descomposición"),
       col=1:2,lty=1:2,lwd=3,cex=0.8)

#Comparación del Ajuste de los 4 modelos
win.graph(width=10,height=10) 
par(mfrow = c(2, 2))
plot(serie,xlab = "Time",ylab='Yt',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythatmodelo1,col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste modelo 1            "),lty=1,col=c(1,2),cex=0.5)
plot(serie,xlab = "Time",ylab='Yt',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythatmodelo2,col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste modelo 2            "),lty=1,col=c(1,2),cex=0.5)
plot(serie,xlab = "Time",ylab='Yt',family = "A",cex.lab=20,cex.axis = 0.5)
lines(fitted(suav),col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste SEHW                "),col=c(1,2),lty=1,cex=0.5)
plot(serie,xlab = "Time",ylab='Yt',family = "A",cex.lab=20,cex.axis = 0.5)
lines(fitted(ajusteDLL1b),col=2,lwd=2)
legend("topleft",legend=c("Original","Ajuste DLL(GCV)            "),col=c(1,2),lty=1,cex=0.5)

#Comparación de los supuestos del error residuales vs tiempo
win.graph(width=10,height=10) 
par(mfrow = c(2, 2))
plot.ts(family = "A",residuals(modelo1),ylim=c(min(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma),max(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma)))
abline(h=c(-2*summary(modelo1)$sigma,0,2*summary(modelo1)$sigma),col=2)
legend("topleft",legend="Modelo 1",cex=0.5)
plot.ts(family = "A",residuals(modelo2),ylim=c(min(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma),max(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma)))
abline(h=c(-2*summary(modelo2)$sigma,0,2*summary(modelo2)$sigma),col=2)
legend("topleft",legend="Modelo 2",cex=0.5)
plot(family = "A",ylab='residuals(modelo 3)',residuals(suav),ylim=c(min(-2*sqrt(MSE3),residuals(suav)),max(2*sqrt(MSE3),residuals(suav))))
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),col=2)
legend("topleft",legend="Modelo 3: SEHW",cex=0.5)
plot(family = "A",ylab='residuals(modelo 4)',residuals(ajusteDLL1b),
     ylim=c(min(-2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b)),max(2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b))))
abline(h=c(-2*sqrt(ajusteDLL1b$MSE),0,2*sqrt(ajusteDLL1b$MSE)),col=2)
legend("topleft",legend="Modelo 4: DLL(GCV)",cex=0.5)

#Comparación de los supuestos del error residuales vs ajustados
win.graph(width=10,height=10) 
par(mfrow = c(2, 2))
plot(xlab='Valores ajustados',family = "A",fitted(modelo1),residuals(modelo1),ylim=c(min(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma),max(residuals(modelo1),-2*summary(modelo1)$sigma,2*summary(modelo1)$sigma)))
abline(h=c(-2*summary(modelo1)$sigma,0,2*summary(modelo1)$sigma),col=2) 
legend("topleft",legend="Modelo 1",cex=0.5)
plot(xlab='Valores ajustados',family = "A",fitted(modelo2),residuals(modelo2),ylim=c(min(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma),max(residuals(modelo2),-2*summary(modelo2)$sigma,2*summary(modelo2)$sigma)))
abline(h=c(-2*summary(modelo2)$sigma,0,2*summary(modelo2)$sigma),col=2) 
legend("topleft",legend="Modelo 2",cex=0.5)
plot(xlab='Valores ajustados',family = "A",as.numeric(fitted(suav)),residuals(suav), ylim=c(min(-2*sqrt(MSE3),residuals(suav)),max(2*sqrt(MSE3),residuals(suav))))
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),col=2)
legend("topleft",legend="Modelo 3",cex=0.5)
plot(xlab='Valores ajustados',family = "A",as.numeric(fitted(ajusteDLL1b)),residuals(ajusteDLL1b),
     ylim=c(min(-2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b)),max(2*sqrt(ajusteDLL1b$MSE),residuals(ajusteDLL1b))))
abline(h=c(-2*sqrt(ajusteDLL1b$MSE),0,2*sqrt(ajusteDLL1b$MSE)),col=2)
legend("topleft",legend="Modelo 4: DLL(GCV)",cex=0.5)

#Comparaci´on gr´afica de los pron´osticos y la serie real en periodos ex-post
win.graph()
plot(family="A",ytf,ylim=c(min(ytf,ytpronmodelo1,ytpronmodelo2,ytpron3,ytpron4),
                max(ytf,ytpronmodelo1,ytpronmodelo2,ytpron3,ytpron4)),type="b",pch=19,lwd=2,xaxt="n")
axis(1,at=time(ytf),labels=c("Jul-18","Ago-18","Sep-18","Oct-18","Nov-18","Dic-18","Ene-19","Feb-19","Mar-19","Abr-19","May-19","Jun-19"))
lines(ytpronmodelo1,lty=2,col="blue",type="b",pch=2,lwd=2)
lines(ytpronmodelo2,lty=3,col="orange",type="b",pch=3,lwd=2)
lines(ytpron3,lty=4,col="darkgreen",type="b",pch=4,lwd=2)
lines(ytpron4,lty=5,col="darkred",type="b",pch=5,lwd=2)
legend("topleft",legend=c("Real","Polinomial estacional grado 4","Polinomial estacional grado 6","Holt-Winters","Descomposición & LOESS"),
                          lty=c(1:5),pch=c(19,2:5),col=c("black","blue","orange","darkgreen","darkred"),lwd=2,cex=0.8)



