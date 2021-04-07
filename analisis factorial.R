library(psych)###paquete de análisis psicométrico
library(GPArotation)###paquete para obtener resultados de métodos de rotación
library(ggplot)##paquete para gráficos
library(ggplot2)## aqute para gráfico

##usaremos a la matriz de correlaciones "Thurstone"

det(Thurstone)###el valor debe reflejar una matriz con alta asociación lineal,
##por tanto se una matriz no singular, el det debe aproximarse a 0

cortest.bartlett(Thurstone)##test de esfericidad de Bartlett 
##Ho nula de que la matriz de correlaciones es igual a una
## matriz de identidad

KMO(Thurstone)##Kaiser-Meyer-Olkin Test de adecuación factorial su valor debe
##debe aproximarse a 1, es aceptable a partir de 0.7

##Porcentajes de variación explicada por cada componente

eigenfact=eigen(Thurstone)

Prop.Var=eigenfact$values/sum(eigenfact$values)*100
cumProp.var=cumsum(eigenfact$values/sum(eigenfact$values)*100)
porc=data.frame(Comp=1:9,Autovalor=round(eigenfact$values,3),Porc.Var=round(Prop.Var,3),Acum.Porc.Var=round(cumProp.var,3))

##variación total explicada

porc

##Cálculos de Comunalidades Matriz de Autovalores(eigen)

diagaut=matrix(diag(eigenfact$values),ncol=9,nrow=9)

##Matriz de Cargas Factoriales

alfa=eigenfact$vectors%*%sqrt(diagaut)

Loaded=alfa%*%t(alfa)

##comunalidad

Loaded2=diag(Loaded)
Loaded2
## 

especificidad

him=alfa[,1:2]%*%t(Loaded[,1:2])

##comunalidad

hi2m=diag(him)
hi2m
## 
##Matriz de Comunalidades

data.frame(varible=names(Thurstone[,1:9]), varible=inicial=round(Loaded2,3),extraccion=round(hi2m,3))
data.frame(inicial=round(Loaded2,3),extraccion=round(hi2m,3))
##Varianza Total

sum(eigenfact$values)

##Gráfico de sedimentación

barplot(eigenfact$values)

##matriz de carga factoriales

data.frame(comp1=round(-alfa[,1],3),comp1=round(alfa[,2],3))

Calculando la Matriz de Puntuaciones Factoriales

#F=bX(ojo para escribir la ecuación)

B=solve(Thurstone)%*%alfa

data.frame(varible=names(Thurstone[,1:9]),Fac1=round(-B[,1],3),Fac2=round(B[,2],3))

B1=data.frame(B)
ggplot(B1,aes(-B1[,1],B1[,2],label=rownames(B1)))+geom_point()+geom_text(vjust = 2)+
xlab("Fact 1")+ylab("Fact 2")+geom_hline(yintercept=0,size=1)+geom_vline(xintercept=0,size=1)

###analisis factorial con roación
##rotate=c("none", "varimax", "quartimax", "bentlerT", "equamax", "varimin", "geominT" and 
##"bifactor" are orthogonal rotations. "Promax", "promax", "oblimin", "simplimax", "bentlerQ, 
##"geominQ" and "biquartimin" and "cluster" )
 
##fm Factoring method fm="minres" will do a minimum residual as will fm="uls". Both of these use
## a first derivative. fm="ols" differs very slightly from "minres" in that it minimizes the 
##entire residual matrix using an OLS procedure but uses the empirical first derivative. 
##This will be slower. fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does
## a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, 
##fm="ml" will do a maximum likelihood factor analysis. fm="minchi" will minimize the sample 
##size weighted chi square when treating pairwise correlations with different number of subjects
## per pair. fm ="minrank" will do a minimum rank factor analysis. "old.min" will do minimal 
##residual the way it was done prior to April, 2017 (see discussion below). fm="alpha" will do
## alpha factor analysis as described in Kaiser and Coffey (1965)
 
##scores the default="regression" finds factor scores using regression. Alternatives for 
##estimating factor scores include simple regression ("Thurstone"), correlaton preserving
##("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms 
##( factor.scores). Although scores="tenBerge" is probably preferred for most solutions, 
##it will lead to problems with some improper correlation matrices. 
 

fac.varimax=fa(Thurstone,nfactors=2,fm="pa",rotate="varimax",scores="TRUE")

plot(fac.varimax,labels=row.names(Thurstone),cex=.7, ylim=c(-.8,.8))


####puntuaciones factoriales

punt.fac=factor.scores(Thurstone,fac.varimax,method="Thurstone")
punt.fac$weights

