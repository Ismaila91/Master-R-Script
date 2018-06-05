###---Ce programme calcule le PCS entre la loi EV1+ et (Normale, Logistique, Student, Gev, Gamma3)
###---en basant la comparaison sur la valeur de la statistique du test---------####



### Suppression des objets existants dans la m?moire ###
rm(list=ls())

###---Répertoire de travail---###
setwd("C:/Users/Ismaila/Documents")

###---Importer le package ADGofTest pour faire le test Anderson-Darling---###
require(ADGofTest)
require(FAdist)
require(data.table)

NumberIterations <- 500
seed <- 2

###---On peut changer shape.gev,shape.gamma3 et dl.student suivant nos besoins---###
shape.gev <- c(-0.2,-0.1,0.1,0.2)
shape.gamma3 <- c(2,3,10)
dl.student <- c(2,10)

###---Sample.size est utilis? pour l'affichage dans Excel, si on change Taille.Echantillon, il faut---###
###---changer aussi Sample.size, au fait les valeurs de Taille.Echantillon doivent etre les memes dans Sample.size---###
Taille.Echantillon <- c(10,20,40,60,80,100)
Sample.size <- rbind("","n",10,20,40,60,80,100)

###---Fonction de ?logvraisemblance pour la loi EV1+ ---###
#(Note: Y ~ EV1- --> x = -Y ~ EV1+)
Vraisemblance.EV<- function(X,parametres){
  n<-length(X)
  return(-(-n*log(parametres[2])-1/parametres[2]*sum(X-parametres[1])-sum(exp(-(X-parametres[1])/parametres[2]))))
}
###---Fonction de ?logvraisemblance pour la loi EV1+ ---###
log.EV <- function(X,parametres){
  -sum(dgumbel(X,scale=parametres[1],location=parametres[2],log=TRUE))
}

###---Fonction de ?logvraisemblance pour la loi GEV(le parametre de forme est fix? a sh)---###
log.GEV <- function(X,sh,parametres){
  -sum(dgev(X,shape=sh,scale=parametres[1],location=parametres[2],log=TRUE))
}

###---Fonction de ?logvraisemblance pour la loi GEV(le parametre de forme est fix? a sh)---###
Vraisemblance.GEV <- function(X,sh,parametres){
  n <- length(X)
  a <- 1-(sh*(X-parametres[2]))/parametres[1]
  return(-(-n*log(parametres[1])+(1/sh-1)*sum(log(a))-sum(a^(1/sh))))
}

###---Fonction de ?logvraisemblance pour la loi GAMMA3(le parametre de forme est fix? a sh) ---###
log.GAMMA3 <- function(X,sh,parametres){
  -sum(dgamma3(X,shape=sh,scale=parametres[1],thres=parametres[2],log=TRUE))
}

###---Fonction de ?logvraisemblance pour la loi NORM---###
Vraisemblance.NORM<- function(X,parametres){
  n<-length(X)
  return(-(-n*log(parametres[2])-(n/2)*log(2*pi)-(1/2)*sum(((X-parametres[1])/parametres[2])^2)))
}

###---Fonction de ?logvraisemblance pour la loi LOG---###
Vraisemblance.LOG<- function(X,parametres){
  n<-length(X)
  return(-(1/parametres[2]*sum(X-parametres[1])-n*log(parametres[2])-2*sum(log(1+exp((X-parametres[1])/parametres[2])))))
}

###---Fonction de logvraisemblance de la loi t de Student d?cal?---###
log.STU.DCL <- function(X,parametres){
  -sum(dt(X-parametres[2],df=parametres[1],ncp=0,log=TRUE))
}


###---Fonction qui cr?? un ?chantillon de loi EV1+ ---###
xEV<-function(n,parametres){
  U <- runif(n,0,1)
  y <- parametres[1]-parametres[2]*log(-log(U))
  return(y)
}

###---Fonction qui cr?? un ?chantillon de loi GEV---###
xGEV <- function(n,parametres){
  U <- runif(n,0,1)
  y <- parametres[3] + (parametres[2]/parametres[1])*(1-(-log(U))^parametres[1])
  return(y)
}

###---Fonction qui cr?? un ?chantillon de loi Normale standard---###
xNOR<-function(n){
  U <- runif(n,0,1)
  y <- qnorm(U)
  return(y)
}

###---Fonction qui cr?? un ?chantillon de loi LOG---###
xLOG<-function(n,parametres){
  U<-runif(n,0,1)
  y<-parametres[1]+parametres[2]*log(U/(1-U))
  return(y)
}


###---Fonction de r?partition de la loi EV1+ ---###
Repartition.EV<-function(X){
  parametres<-optim(par=c(0,1),fn=Vraisemblance.EV,X=X)$par
  return( exp(-exp(-(X-parametres[1])/parametres[2])) )
}

###---Fonction de r?partition de la loi Normale standard---###
Repartition.NORM<-function(X){
  parametres <- optim(par=c(0,1),fn=Vraisemblance.NORM,X=X)$par
  return(pnorm(X,parametres[1],parametres[2]))
}

###---Fonction de r?partition de la loi LOG---###
Repartition.LOG <-function(X){
  parametres<-optim(par=c(0,1),fn=Vraisemblance.LOG,X=X)$par
  return(exp((X-parametres[1])/parametres[2])/(1+exp((X-parametres[1])/parametres[2])))
}

########----MLE de GEV--------########
optim.gev <- function(shape){
  return(optim(c(1.12,0.07),log.GEV,sh=shape,X=X)$par) 
}

####-----MLE de GAMMA3---------#########
optim.gam3 <- function(shape){
  return(optim(c(1,-10),log.GAMMA3,sh=shape,X=X)$par)
}

########----MLE de GEV tri?--------########
optim.gev.sort <- function(shape){
  return(optim(c(1.12,0.07),log.GEV,sh=shape,X=sort(X))$par) 
}

####-----MLE de GAMMA3 tri?---------#########
optim.gam3.sort <- function(shape){
  return(optim(c(1,-10),log.GAMMA3,sh=shape,X=sort(X))$par)
}

#####------Transformation en normalit? de GEV-------#########
z.gev.sort <- function(parametres){
  return(qnorm(pgev(sort(X),parametres[1],parametres[2],parametres[3])))
}

#####------Transformation en normalit? de GAM3-------#########
z.gam3.sort <- function(parametres){
  return(qnorm(pgamma3(sort(X),parametres[1],parametres[2],parametres[3])))
}

#####------Calcul de la statistique de Anderson-Darling pour GEV-----######
ander.darling.gev <- function(parametres){
  return(ad.test(X,pgev,parametres[1],parametres[2],parametres[3])$statistic)
}

#####------Calcul de la statistique de Anderson-Darling pour GAMMA3-----######
ander.darling.gam3 <- function(parametres){
  return(ad.test(X,pgamma3,parametres[1],parametres[2],parametres[3])$statistic)
}

#####------Calcul de la statistique TNSW pour GEV-------#########
tnsw.gev <- function(parametres){
  return(shapiro.test(qnorm(pgev(X,parametres[1],parametres[2],parametres[3])))$statistic)
}

#####------Calcul de la statistique TNSW pour GAMMA3-------#########
tnsw.gam3 <- function(parametres){
  return(shapiro.test(qnorm(pgamma3(X,parametres[1],parametres[2],parametres[3])))$statistic)
}

####----Calcul de la log-vraisemblance pour GEV-----#####
ll.gev <- function(shape){
  return(-optim(c(1.12,0.07),log.GEV,sh=shape,X=X)$value) 
}

####----Calcul de la log-vraisemblance pour GAMMA3-----#####
ll.gam3 <- function(shape){
  return(-optim(c(1,-10),log.GAMMA3,sh=shape,X=X)$value)
}


############################################################################################################################

###-----Discrimination between the Gumbel and the Standard Normal Distributions-----------###

############-----------M?thode RML--------#################
###---PCS(?chantillon Gumbel, alternatif normale) avec la m?thode RML---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.NOR <- rep(0,Nombre.Iteration)
PCS.RML.GUM.NOR <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    LLGUM <- -optim(par=c(1,0),log.EV,X=X)$value
    LLNORM <- -optim(par=c(0,1),fn=Vraisemblance.NORM,X=X)$value
    if(LLGUM-LLNORM>0){res.GUM.NOR[i] <- 1}else{res.GUM.NOR[i] <- 0}
  }
  PCS.RML.GUM.NOR[j] <- sum(res.GUM.NOR)/Nombre.Iteration
  j <- j+1
}

PCS.RML.GUM.NOR <- round(PCS.RML.GUM.NOR,2)*100

###---PCS(?chantillon Normal, alternatif Gumbel) avec la m?thode RML---###
j <- 1
res.NOR.GUM <- rep(0,Nombre.Iteration)
PCS.RML.NOR.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xNOR(n)
    LLGUM <- -optim(par=c(0,1),fn=Vraisemblance.EV,X=X)$value
    LLNORM <- -optim(par=c(0,1),fn=Vraisemblance.NORM,X=X)$value
    if(LLGUM-LLNORM<0){res.NOR.GUM[i] <- 1}else{res.NOR.GUM[i] <- 0}
  }
  PCS.RML.NOR.GUM[j] <- sum(res.NOR.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.RML.NOR.GUM <- round(PCS.RML.NOR.GUM,2)*100

GUM.NOR.RML <- cbind(Taille.Echantillon,PCS.RML.GUM.NOR,PCS.RML.NOR.GUM)
GUM.NOR.RML

###########-----------M?thode TNSW-----------###############
###---PCS(?chantillon Gumbel, alternatif normale) avec la m?thode TN.SW---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.NOR <- rep(0,Nombre.Iteration)
PCS.TNSW.GUM.NOR <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    # Calcul de la statistique du test pour la loi Gumbel
    Z.GUM <- qnorm(Repartition.EV(X))
    S.GUM <- shapiro.test(Z.GUM)$statistic
    # Calcul de la statistique du test pour la loi Normale
    Z.NORM <- qnorm(Repartition.NORM(X))
    S.NORM <- shapiro.test(Z.NORM)$statistic
    if(S.NORM<S.GUM){res.GUM.NOR[i] <- 1}else{res.GUM.NOR[i] <- 0}
  }
  PCS.TNSW.GUM.NOR[j] <- sum(res.GUM.NOR)/Nombre.Iteration
  j <- j+1
}

PCS.TNSW.GUM.NOR <- round(PCS.TNSW.GUM.NOR,2)*100

###---PCS(?chantillon Normal, alternatif Gumbel) avec la m?thode TN.SW---###
j <- 1
res.NOR.GUM <- rep(0,Nombre.Iteration)
PCS.TNSW.NOR.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xNOR(n)
    # Calcul de la statistique du test pour la loi Normale
    S.NORM <- shapiro.test(X)$statistic
    # Calcul de la statistique pour la loi Gumbel
    Z.GUM <- qnorm(Repartition.EV(X))
    S.GUM <- shapiro.test(Z.GUM)$statistic
    if(S.GUM<S.NORM){res.NOR.GUM[i] <- 1}else{res.NOR.GUM[i] <- 0}
  }
  PCS.TNSW.NOR.GUM[j] <- sum(res.NOR.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.TNSW.NOR.GUM <- round(PCS.TNSW.NOR.GUM,2)*100

GUM.NOR.TNSW <- cbind(Taille.Echantillon,PCS.TNSW.GUM.NOR,PCS.TNSW.NOR.GUM)
GUM.NOR.TNSW

#########---------------M?thode Anderson-darling-----------##########
###---PCS(?chantillon Gumbel, alternatif normale) avec la m?thode AD---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.NOR <- rep(0,Nombre.Iteration)
PCS.AD.GUM.NOR <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    AGUM <- ad.test(X,Repartition.EV)$statistic
    ANOR <- ad.test(X,Repartition.NORM)$statistic
    if(ANOR>AGUM){res.GUM.NOR[i] <- 1}else{res.GUM.NOR[i] <- 0}
    
  }
  PCS.AD.GUM.NOR[j] <- sum(res.GUM.NOR)/Nombre.Iteration
  j <- j+1
}

PCS.AD.GUM.NOR <- round(PCS.AD.GUM.NOR,2)*100

###---PCS(?chantillon Normal, alternatif Gumbel) avec la m?thode AD---###
j <- 1
res.NOR.GUM <- rep(0,Nombre.Iteration)
PCS.AD.NOR.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xNOR(n)
    ANOR <- ad.test(X,Repartition.NORM)$statistic
    AGUM <- ad.test(X,Repartition.EV)$statistic
    if(AGUM>ANOR){res.NOR.GUM[i] <- 1}else{res.NOR.GUM[i] <- 0}
  }
  PCS.AD.NOR.GUM[j] <- sum(res.NOR.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.AD.NOR.GUM <- round(PCS.AD.NOR.GUM,2)*100

GUM.NOR.AD <- cbind(Taille.Echantillon,PCS.AD.GUM.NOR,PCS.AD.NOR.GUM)
GUM.NOR.AD

#########---------------M?thode PPCC Modifi? (TN.PPCC) Test-----------##########
###---PCS(?chantillon Gumbel, alternatif Normal) avec la m?thode tn.ppcc ---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.NOR <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.GUM.NOR <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  Quantiles <- (1:n)/(n+1)
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    Z.GUM.sort <- qnorm(Repartition.EV(sort(X)))
    Z.NORM.sort <- qnorm(Repartition.NORM(sort(X)))
    w <- qnorm(Quantiles)
    tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
    tn_ppcc.Normal <- sum((Z.NORM.sort-mean(Z.NORM.sort))*(w-mean(w)))/sqrt(sum(((sort(Z.NORM.sort)-mean(Z.NORM.sort))^2))*sum(((w-mean(w))^2)))
    if(tn_ppcc.Gumbel>tn_ppcc.Normal){res.GUM.NOR[i] <- 1}else{res.GUM.NOR[i] <- 0}
  }
  PCS.TN_PPCC.GUM.NOR[j] <- sum(res.GUM.NOR)/Nombre.Iteration
  j <- j+1
}

PCS.TN_PPCC.GUM.NOR <- round(PCS.TN_PPCC.GUM.NOR,2)*100

###---PCS(?chantillon Normal, alternatif Gumbel) avec la m?thode tn.ppcc---###
j <- 1
res.NOR.GUM <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.NOR.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  Quantiles <- (1:n)/(n+1)
  for(i in 1:Nombre.Iteration){
    X <- xNOR(n)
    Z.GUM.sort <- qnorm(Repartition.EV(sort(X)))
    w <- qnorm(Quantiles)
    tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
    tn_ppcc.Normal <- sum((sort(X)-mean(X))*(w-mean(w)))/sqrt(sum(((sort(X)-mean(X))^2))*sum(((w-mean(w))^2)))
    if(tn_ppcc.Normal>tn_ppcc.Gumbel){res.NOR.GUM[i] <- 1}else{res.NOR.GUM[i] <- 0}
  }
  PCS.TN_PPCC.NOR.GUM[j] <- sum(res.NOR.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.TN_PPCC.NOR.GUM <- round(PCS.TN_PPCC.NOR.GUM,2)*100

GUM.NOR.TN_PPCC <- cbind(Taille.Echantillon,PCS.TN_PPCC.GUM.NOR,PCS.TN_PPCC.NOR.GUM)
GUM.NOR.TN_PPCC


####################################################################################################################################################################

###-----Discrimination between the Gumbel and the Logistic Distributions-----------###

############-----------M?thode RML--------#################
###---PCS(?chantillon Gumbel, alternatif logistique) avec la m?thode RML---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.LOG <- rep(0,Nombre.Iteration)
PCS.RML.GUM.LOG <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    LLGUM <- -optim(par=c(1,0),log.EV,X=X)$value
    LLLOG <- -optim(par=c(0,1),fn=Vraisemblance.LOG,X=X)$value
    if(LLGUM-LLLOG>0){res.GUM.LOG[i] <- 1}else{res.GUM.LOG[i] <- 0}
  }
  PCS.RML.GUM.LOG[j] <- sum(res.GUM.LOG)/Nombre.Iteration
  j <- j+1
}

PCS.RML.GUM.LOG <- round(PCS.RML.GUM.LOG,2)*100

###---PCS(?chantillon Logistique, alternatif Gumbel) avec la m?thode RML---###
j <- 1
res.LOG.GUM <- rep(0,Nombre.Iteration)
PCS.RML.LOG.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xLOG(n,c(0,1))
    LLGUM <- -optim(par=c(0,1),fn=Vraisemblance.EV,X=X)$value
    LLLOG <- -optim(par=c(0,1),fn=Vraisemblance.LOG,X=X)$value
    if(LLGUM-LLLOG<0){res.LOG.GUM[i] <- 1}else{res.LOG.GUM[i] <- 0}
  }
  PCS.RML.LOG.GUM[j] <- sum(res.LOG.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.RML.LOG.GUM <- round(PCS.RML.LOG.GUM,2)*100

GUM.LOG.RML <- cbind(Taille.Echantillon,PCS.RML.GUM.LOG,PCS.RML.LOG.GUM)
GUM.LOG.RML

###########-----------M?thode TNSW-----------###############
###---PCS(?chantillon Gumbel, alternatif logistique) avec la m?thode TN.SW---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.LOG <- rep(0,Nombre.Iteration)
PCS.TNSW.GUM.LOG <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    # Calcul de la statistique du test pour la loi Gumbel
    Z.GUM <- qnorm(Repartition.EV(X))
    S.GUM <- shapiro.test(Z.GUM)$statistic
    # Calcul de la statistique du test pour la loi Logistique
    Z.LOG <- qnorm(Repartition.LOG(X))
    S.LOG <- shapiro.test(Z.LOG)$statistic
    if(S.LOG<S.GUM){res.GUM.LOG[i] <- 1}else{res.GUM.LOG[i] <- 0}
  }
  PCS.TNSW.GUM.LOG[j] <- sum(res.GUM.LOG)/Nombre.Iteration
  j <- j+1
}

PCS.TNSW.GUM.LOG <- round(PCS.TNSW.GUM.LOG,2)*100

###---PCS(?chantillon Logistique, alternatif Gumbel) avec la m?thode TN.SW---###
j <- 1
res.LOG.GUM <- rep(0,Nombre.Iteration)
PCS.TNSW.LOG.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xLOG(n,c(0,1))
    # Calcul de la statistique de test pour la loi Logistique
    Z.LOG <- qnorm(Repartition.LOG(X))
    S.LOG <- shapiro.test(Z.LOG)$statistic
    # Calcul de la statistique de test pour la loi Gumbel
    Z.GUM <- qnorm(Repartition.EV(X))
    S.GUM <- shapiro.test(Z.GUM)$statistic
    if(S.GUM<S.LOG){res.LOG.GUM[i] <- 1}else{res.LOG.GUM[i] <- 0}
  }
  PCS.TNSW.LOG.GUM[j] <- sum(res.LOG.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.TNSW.LOG.GUM <- round(PCS.TNSW.LOG.GUM,2)*100

GUM.LOG.TNSW <- cbind(Taille.Echantillon,PCS.TNSW.GUM.LOG,PCS.TNSW.LOG.GUM)
GUM.LOG.TNSW

#########---------------M?thode Anderson-darling-----------##########
###---PCS(?chantillon Gumbel, alternatif logistique) avec la m?thode AD---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.LOG <- rep(0,Nombre.Iteration)
PCS.AD.GUM.LOG <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    ALOG <-ad.test(X,Repartition.LOG)$statistic
    AGUM <- ad.test(X,Repartition.EV)$statistic
    if(ALOG>AGUM){res.GUM.LOG[i] <- 1}else{res.GUM.LOG[i] <- 0}	
  }
  PCS.AD.GUM.LOG[j] <- sum(res.GUM.LOG)/Nombre.Iteration
  j <- j+1
}

PCS.AD.GUM.LOG <- round(PCS.AD.GUM.LOG,2)*100

###---PCS(?chantillon Logistique, alternatif Gumbel) avec la m?thode AD---###
j <- 1
res.LOG.GUM <- rep(0,Nombre.Iteration)
PCS.AD.LOG.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xLOG(n,c(0,1))
    ALOG <-ad.test(X,Repartition.LOG)$statistic
    AGUM <- ad.test(X,Repartition.EV)$statistic
    if(AGUM>ALOG){res.LOG.GUM[i] <- 1}else{res.LOG.GUM[i] <- 0}
  }
  PCS.AD.LOG.GUM[j] <- sum(res.LOG.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.AD.LOG.GUM <- round(PCS.AD.LOG.GUM,2)*100

GUM.LOG.AD <- cbind(Taille.Echantillon,PCS.AD.GUM.LOG,PCS.AD.LOG.GUM)
GUM.LOG.AD

#########---------------M?thode TN.PPCC-----------##########
###---PCS(?chantillon Gumbel, alternatif logistique) avec la m?thode tn.ppcc---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.LOG <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.GUM.LOG <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  Quantiles <- (1:n)/(n+1)
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    Z.GUM.sort <- qnorm(Repartition.EV(sort(X)))
    Z.LOG.sort <- qnorm(Repartition.LOG(sort(X)))
    w <- qnorm(Quantiles)
    tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
    tn_ppcc.Logis <- sum((Z.LOG.sort-mean(Z.LOG.sort))*(w-mean(w)))/sqrt(sum(((Z.LOG.sort-mean(Z.LOG.sort))^2))*sum(((w-mean(w))^2)))
    if(tn_ppcc.Gumbel>tn_ppcc.Logis){res.GUM.LOG[i] <- 1}else{res.GUM.LOG[i] <- 0}
  }
  PCS.TN_PPCC.GUM.LOG[j] <- sum(res.GUM.LOG)/Nombre.Iteration
  j <- j+1
}

PCS.TN_PPCC.GUM.LOG <- round(PCS.TN_PPCC.GUM.LOG,2)*100

###---PCS(?chantillon Logistique, alternatif Gumbel) avec la m?thode tn.ppcc---###
j <- 1
res.LOG.GUM <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.LOG.GUM <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  Quantiles <- (1:n)/(n+1)
  for(i in 1:Nombre.Iteration){
    X <- xLOG(n,c(0,1))
    Z.LOG.sort <- qnorm(Repartition.LOG(sort(X)))
    Z.GUM.sort <- qnorm(Repartition.EV(sort(X)))
    w <- qnorm(Quantiles)
    tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
    tn_ppcc.Logis <- sum((Z.LOG.sort-mean(Z.LOG.sort))*(w-mean(w)))/sqrt(sum(((Z.LOG.sort-mean(Z.LOG.sort))^2))*sum(((w-mean(w))^2)))
    if(tn_ppcc.Logis>tn_ppcc.Gumbel){res.LOG.GUM[i] <- 1}else{res.LOG.GUM[i] <- 0}
  }
  PCS.TN_PPCC.LOG.GUM[j] <- sum(res.LOG.GUM)/Nombre.Iteration
  j <- j+1
}

PCS.TN_PPCC.LOG.GUM <- round(PCS.TN_PPCC.LOG.GUM,2)*100

GUM.LOG.TN_PPCC <- cbind(Taille.Echantillon,PCS.TN_PPCC.GUM.LOG,PCS.TN_PPCC.LOG.GUM)
GUM.LOG.TN_PPCC



#############################################################################################################################################

###-------Discrimination between Gumbel and Pearson type 3 distributions--------###

############-----------M?thode RML--------#################
###---PCS(?chantillon Gumbel, alternatif  pearson type 3) avec la m?thode RML---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GAM3 <-  array(0,c(Nombre.Iteration,length(shape.gamma3)))
PCS.RML.GUM.GAM3 <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xEV(n,c(0,1))
    LLGAM3 <- sapply(shape.gamma3,ll.gam3)
    LLGUM <- -optim(par=c(1,0),log.EV,X=X)$value
    for(l in 1:length(shape.gamma3)) if(LLGUM-LLGAM3[l]>0){res.GUM.GAM3[i,l] <- 1}else{res.GUM.GAM3[i,l] <- 0}
  }
  for(l in 1:length(shape.gamma3)) PCS.RML.GUM.GAM3[l,j] <- sum(res.GUM.GAM3[,l])/Nombre.Iteration
  j <- j+1
}

PCS.RML.GUM.GAM3 <- round(PCS.RML.GUM.GAM3,2)*100

###---PCS(?chantillon Pearson type 3, alternatif Gumbel) avec la m?thode RML---###
m_gam3_rml <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
res.GAM3.GUM <- rep(0,Nombre.Iteration)
PCS.RML.GAM3.GUM <- rep(0,length(shape.gamma3))
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  for(a in shape.gamma3){
    for(i in 1:Nombre.Iteration){
      X <- rgamma3(n,a,1,0)
      LLGAM3 <- -optim(c(1,-10),log.GAMMA3,sh=a,X=X)$value
      LLGUM <- -optim(c(1,0),log.EV,X=X)$value
      if(LLGAM3-LLGUM>0){res.GAM3.GUM[i] <- 1}else{res.GAM3.GUM[i] <- 0}
    }
    PCS.RML.GAM3.GUM[j] <- sum(res.GAM3.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gam3_rml[,r]<-PCS.RML.GAM3.GUM
  r<-r+1
}

m_gam3_rml <- round(m_gam3_rml,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gamma3)){
  print(cbind(PCS.RML.GUM.GAM3[l,],m_gam3_rml[l,]))
}

###########-----------M?thode TNSW-----------###############
###---PCS(?chantillon Gumbel, alternatif pearson type 3) avec la m?thode TN.SW---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GAM3 <-  array(0,c(Nombre.Iteration,length(shape.gamma3)))
PCS.TNSW.GUM.GAM3 <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- xEV(n,c(0,1))
    # Calcul de la statistique du test pour la loi Pearson Type 3
    S.GAM3 <- apply(rbind(shape.gamma3,sapply(shape.gamma3,optim.gam3)),2,tnsw.gam3)
    # Calcul de la statistique du test pour la loi Gumbel
    S.GUM <- shapiro.test(qnorm(Repartition.EV(X)))$statistic	
    for(l in 1:length(shape.gamma3)) if(S.GUM>S.GAM3[l]){res.GUM.GAM3[i,l] <- 1}else{res.GUM.GAM3[i,l] <- 0}
  }
  for(l in 1:length(shape.gamma3)) PCS.TNSW.GUM.GAM3[l,j] <- sum(res.GUM.GAM3[,l])/Nombre.Iteration
  j <- j+1
}

PCS.TNSW.GUM.GAM3 <- round(PCS.TNSW.GUM.GAM3,2)*100

###---PCS(?chantillon Pearson type 3, alternatif Gumbel) avec la m?thode TNSW---###
m_gam3_tnsw <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
res.GAM3.GUM <- rep(0,Nombre.Iteration)
PCS.TNSW.GAM3.GUM <- rep(0,length(shape.gamma3))
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  for(a in shape.gamma3){
    for(i in 1:Nombre.Iteration){
      X <- rgamma3(n,a,1,0)
      param.GAM3 <- optim(c(1,-10),log.GAMMA3,sh=a,X=X)$par
      z.hat.gam3 <- qnorm(pgamma3(X,a,param.GAM3[1],param.GAM3[2]))
      S.GAM3 <- shapiro.test(z.hat.gam3)$statistic
      param.GUM <- optim(c(1,0),log.EV,X=X)$par
      z.hat.gumbel <- qnorm(pgumbel(X,param.GUM[1],param.GUM[2]))
      S.GUM <- shapiro.test(z.hat.gumbel)$statistic
      if(S.GAM3>S.GUM){res.GAM3.GUM[i] <- 1}else{res.GAM3.GUM[i] <- 0}
    }
    PCS.TNSW.GAM3.GUM[j] <- sum(res.GAM3.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gam3_tnsw[,r]<-PCS.TNSW.GAM3.GUM
  r<-r+1
}

m_gam3_tnsw <- round(m_gam3_tnsw,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gamma3)){
  print(cbind(PCS.TNSW.GUM.GAM3[l,],m_gam3_tnsw[l,]))
}

#########---------------M?thode Anderson-darling-----------##########
###---PCS(?chantillon Gumbel, alternatif Pearson type 3 ) avec la m?thode AD---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GAM3 <-  array(0,c(Nombre.Iteration,length(shape.gamma3)))
PCS.AD.GUM.GAM3 <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    AGAM3 <- apply(rbind(shape.gamma3,sapply(shape.gamma3,optim.gam3)),2,ander.darling.gam3)		
    AGUM <- ad.test(X,Repartition.EV)$statistic		
    for(l in 1:length(shape.gamma3)) if(AGUM<AGAM3[l]){res.GUM.GAM3[i,l] <- 1}else{res.GUM.GAM3[i,l] <- 0}	
  }
  for(l in 1:length(shape.gamma3)) PCS.AD.GUM.GAM3[l,j] <- sum(res.GUM.GAM3[,l])/Nombre.Iteration
  j <- j+1
}

PCS.AD.GUM.GAM3 <- round(PCS.AD.GUM.GAM3,2)*100

###---PCS(?chantillon Pearson type 3, alternatif Gumbel) avec la m?thode AD---###
m_gam3_ad <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
res.GAM3.GUM <- rep(0,Nombre.Iteration)
PCS.AD.GAM3.GUM <- rep(0,length(shape.gamma3))
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  for(a in shape.gamma3){
    for(i in 1:Nombre.Iteration){
      X <- rgamma3(n,a,1,0)
      param.GAM3 <- optim(c(1,-10),log.GAMMA3,sh=a,X=X)$par
      AGAM3 <- ad.test(X,pgamma3,a,param.GAM3[1],param.GAM3[2])$statistic
      param.GUM <- optim(c(1,0),log.EV,X=X)$par
      AGUM <- ad.test(X,pgumbel,param.GUM[1],param.GUM[2])$statistic
      if(AGAM3<AGUM){res.GAM3.GUM[i] <- 1}else{res.GAM3.GUM[i] <- 0}
    }
    PCS.AD.GAM3.GUM[j] <- sum(res.GAM3.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gam3_ad[,r]<-PCS.AD.GAM3.GUM
  r<-r+1
}

m_gam3_ad <- round(m_gam3_ad,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gamma3)){
  print(cbind(PCS.AD.GUM.GAM3[l,],m_gam3_ad[l,]))
}


#########---------------M?thode TN.PPCC-----------##########
###---PCS(?chantillon Gumbel, alternatif Pearson type 3 ) avec la m?thode tn.ppcc---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GAM3 <-  array(0,c(Nombre.Iteration,length(shape.gamma3)))
PCS.TN_PPCC.GUM.GAM3 <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  Quantiles <- (1:n)/(n+1)
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    Z.GAM3.sort <- apply(rbind(shape.gamma3,sapply(shape.gamma3,optim.gam3.sort)),2,z.gam3.sort)
    Z.GUM.sort <- qnorm(Repartition.EV(sort(X)))
    w <- qnorm(Quantiles)
    tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
    tn_ppcc.Gamma3 <- rep(0,length(shape.gamma3))
    for(l in 1:length(shape.gamma3)) {
      numer <- sum((Z.GAM3.sort[,l]-mean(Z.GAM3.sort[,l]))*(w-mean(w)))
      denom <- sqrt( sum(((Z.GAM3.sort[,l]-mean(Z.GAM3.sort[,l]))^2))* sum(((w-mean(w))^2))) 
      tn_ppcc.Gamma3[l] <- numer/denom
    }
    for(l in 1:length(shape.gamma3)) if(tn_ppcc.Gumbel>tn_ppcc.Gamma3[l]){res.GUM.GAM3[i,l] <- 1}else{res.GUM.GAM3[i,l] <- 0}
  }
  for(l in 1:length(shape.gamma3)) PCS.TN_PPCC.GUM.GAM3[l,j] <- sum(res.GUM.GAM3[,l])/Nombre.Iteration
  j <- j+1
}

PCS.TN_PPCC.GUM.GAM3 <- round(PCS.TN_PPCC.GUM.GAM3,2)*100

###---PCS(?chantillon Pearson type 3, alternatif Gumbel) avec la m?thode tn.ppcc---###
m_gam3_tn.ppcc <- array(0,c(length(shape.gamma3),length(Taille.Echantillon)),dimnames=list(shape.gamma3,Taille.Echantillon))
res.GAM3.GUM <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.GAM3.GUM <- rep(0,length(shape.gamma3))
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  Quantiles <- (1:n)/(n+1)
  for(a in shape.gamma3){
    for(i in 1:Nombre.Iteration){
      X <- rgamma3(n,a,1,0)
      param.GAM3 <- optim(c(1,-10),log.GAMMA3,sh=a,X=sort(X))$par
      Z.GAM3.sort <- qnorm(pgamma3(sort(X),a,param.GAM3[1],param.GAM3[2]))
      param.GUM <- optim(c(1,0),log.EV,X=sort(X))$par
      Z.GUM.sort <- qnorm(pgumbel(sort(X),param.GUM[1],param.GUM[2]))
      w <- qnorm(Quantiles)
      tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
      tn_ppcc.Gamma3 <- sum((Z.GAM3.sort-mean(Z.GAM3.sort))*(w-mean(w)))/sqrt(sum(((Z.GAM3.sort-mean(Z.GAM3.sort))^2))*sum(((w-mean(w))^2)))
      if(tn_ppcc.Gamma3>tn_ppcc.Gumbel){res.GAM3.GUM[i] <- 1}else{res.GAM3.GUM[i] <- 0}
    }
    PCS.TN_PPCC.GAM3.GUM[j] <- sum(res.GAM3.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gam3_tn.ppcc[,r]<-PCS.TN_PPCC.GAM3.GUM
  r<-r+1
}

m_gam3_tn.ppcc <- round(m_gam3_tn.ppcc,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gamma3)){
  print(cbind(PCS.TN_PPCC.GUM.GAM3[l,],m_gam3_tn.ppcc[l,]))
}


####################################################################################################################

###-------Discrimination between Gumbel and Gev distributions---------###

############-----------M?thode RML--------#################
###---PCS(?chantillon Gumbel, alternatif Gev ) avec la m?thode RML---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GEV <-  array(0,c(Nombre.Iteration,length(shape.gev)))
PCS.RML.GUM.GEV <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    LLGEV <- NULL
    while(is.null(LLGEV)){
      X <- xEV(n,c(0,1)) 
      try(LLGEV <- sapply(shape.gev,ll.gev),silent=TRUE)
    }
    LLGUM <- -optim(par=c(1,0),log.EV,X=X)$value
    for(l in 1:length(shape.gev)) if(LLGUM-LLGEV[l]>0){res.GUM.GEV[i,l] <- 1}else{res.GUM.GEV[i,l] <- 0}
  }
  for(l in 1:length(shape.gev)) PCS.RML.GUM.GEV[l,j] <- sum(res.GUM.GEV[,l])/Nombre.Iteration
  j <- j+1
}

PCS.RML.GUM.GEV <- round(PCS.RML.GUM.GEV,2)*100

###---PCS(?chantillon GEV, alternatif Gumbel) avec la m?thode RML---###
m_gev_rml <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon)) 
res.GEV.GUM <- rep(0,Nombre.Iteration)
PCS.RML.GEV.GUM <- rep(0,length(shape.gev))
set.seed(seed)
r<-1
for(n in Taille.Echantillon){
  j <- 1
  for(a in shape.gev){
    for(i in 1:Nombre.Iteration){
      X <- rgev(n,a,1,0)
      LLGEV <- -optim(c(1.12,0.07),log.GEV,sh=a,X=X)$value
      LLGUM <- -optim(c(1,0),log.EV,X=X)$value
      if(LLGEV-LLGUM>0){res.GEV.GUM[i] <- 1}else{res.GEV.GUM[i] <- 0}
    }
    PCS.RML.GEV.GUM[j] <- sum(res.GEV.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gev_rml[,r]<-PCS.RML.GEV.GUM
  r<-r+1
}

m_gev_rml <- round(m_gev_rml,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gev)){
  print(cbind(PCS.RML.GUM.GEV[l,],m_gev_rml[l,]))
}


###########-----------M?thode TNSW-----------###############
###---PCS(?chantillon Gumbel, alternatif  Gev ) avec la m?thode TN.SW---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GEV <-  array(0,c(Nombre.Iteration,length(shape.gev)))
PCS.TNSW.GUM.GEV <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    S.GEV <- NULL
    while(is.null(S.GEV)){
      X <- xEV(n,c(0,1))
      try(S.GEV <- apply(rbind(shape.gev,sapply(shape.gev,optim.gev)),2,tnsw.gev),silent=TRUE)		
    }
    # Calcul de la statistique du test pour la loi Gumbel
    S.GUM <- shapiro.test(qnorm(Repartition.EV(X)))$statistic		
    for(l in 1:length(shape.gev)) if(S.GUM>S.GEV[l]){res.GUM.GEV[i,l] <- 1}else{res.GUM.GEV[i,l] <- 0}
  }
  for(l in 1:length(shape.gev)) PCS.TNSW.GUM.GEV[l,j] <- sum(res.GUM.GEV[,l])/Nombre.Iteration
  j <- j+1
}

PCS.TNSW.GUM.GEV <- round(PCS.TNSW.GUM.GEV,2)*100

###---PCS(?chantillon GEV, alternatif Gumbel) avec la m?thode TNSW---###
m_gev_tnsw <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon))
res.GEV.GUM <- rep(0,Nombre.Iteration)
PCS.TNSW.GEV.GUM <- rep(0,length(shape.gev))
set.seed(seed)
Val.Stat.Gum <- rep(0,Nombre.Iteration)
Val.Stat.Gev <- rep(0,Nombre.Iteration)
r<-1
for(n in Taille.Echantillon){
  j <- 1
  for(a in shape.gev){
    for(i in 1:Nombre.Iteration){
      S.GEV <- NULL
      while(is.null(S.GEV)){ 
        X <- rgev(n,a,1,0)
        param.GEV <- optim(c(1.12,0.07),log.GEV,sh=a,X=X)$par
        z.hat.gev <- qnorm(pgev(X,a,param.GEV[1],param.GEV[2]))
        try(S.GEV <- as.numeric(shapiro.test(z.hat.gev)$statistic),silent=TRUE) 
      }
      param.GUM <- optim(c(1,0),log.EV,X=X)$par
      z.hat.gumbel <- qnorm(pgumbel(X,param.GUM[1],param.GUM[2]))
      S.GUM <- as.numeric(shapiro.test(z.hat.gumbel)$statistic)
      Val.Stat.Gev[i] <- S.GEV
      Val.Stat.Gum[i] <- S.GUM
    }
    Val.Stat.Gum <- ifelse(is.na(Val.Stat.Gum),0,Val.Stat.Gum)
    res.GEV.GUM <- ifelse(Val.Stat.Gum<Val.Stat.Gev,1,0)
    PCS.TNSW.GEV.GUM[j] <- sum(res.GEV.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gev_tnsw[,r]<-PCS.TNSW.GEV.GUM
  r<-r+1
}

m_gev_tnsw <- round(m_gev_tnsw,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gev)){
  print(cbind(PCS.TNSW.GUM.GEV[l,],m_gev_tnsw[l,]))
}

#########---------------M?thode Anderson-darling-----------##########
###---PCS(?chantillon Gumbel, alternatif Gev ) avec la m?thode AD---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GEV <-  array(0,c(Nombre.Iteration,length(shape.gev)))
PCS.AD.GUM.GEV <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    AGEV <- NULL
    while(is.null(AGEV)){
      X <- rgumbel(n,1,0)
      try(AGEV <- apply(rbind(shape.gev,sapply(shape.gev,optim.gev)),2,ander.darling.gev),silent=TRUE)
    }
    AGUM <- ad.test(X,Repartition.EV)$statistic	
    for(l in 1:length(shape.gev)) if(AGUM<AGEV[l]){res.GUM.GEV[i,l] <- 1}else{res.GUM.GEV[i,l] <- 0}
    
  }
  for(l in 1:length(shape.gev)) PCS.AD.GUM.GEV[l,j] <- sum(res.GUM.GEV[,l])/Nombre.Iteration
  j <- j+1
}

PCS.AD.GUM.GEV <- round(PCS.AD.GUM.GEV,2)*100

###---PCS(?chantillon GEV, alternatif Gumbel) avec la m?thode AD---###
m_gev_ad <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon))
res.GEV.GUM <- rep(0,Nombre.Iteration)
PCS.AD.GEV.GUM <- rep(0,length(shape.gev))
set.seed(seed)
r<-1
for(n in Taille.Echantillon){
  j <- 1
  for(a in shape.gev){
    for(i in 1:Nombre.Iteration){
      AGEV <- NULL
      while(is.null(AGEV)){
        X <- rgev(n,a,1,0)
        param.GEV <- optim(c(1.12,0.07),log.GEV,sh=a,X=X)$par
        try(AGEV <- ad.test(X,pgev,a,param.GEV[1],param.GEV[2])$statistic,silent=TRUE)
      }
      param.GUM <- optim(c(1,0),log.EV,X=X)$par
      AGUM <- ad.test(X,pgumbel,param.GUM[1],param.GUM[2])$statistic
      if(AGUM>AGEV){res.GEV.GUM[i] <- 1}else{res.GEV.GUM[i] <- 0}
    }
    PCS.AD.GEV.GUM[j] <- sum(res.GEV.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gev_ad[,r]<-PCS.AD.GEV.GUM
  r<-r+1
}

m_gev_ad <- round(m_gev_ad,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gev)){
  print(cbind(PCS.AD.GUM.GEV[l,],m_gev_ad[l,]))
}



#########---------------M?thode TN.PPCC-----------##########
###---PCS(?chantillon Gumbel, alternatif Gev ) avec la m?thode tn.ppcc---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.GEV <-  array(0,c(Nombre.Iteration,length(shape.gev)))
PCS.TN_PPCC.GUM.GEV <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  Quantiles <- (1:n)/(n+1)
  for(i in 1:Nombre.Iteration){
    Z.GEV.sort <- NULL
    while(is.null(Z.GEV.sort)){
      X <- xEV(n,c(0,1))
      try(Z.GEV.sort <- apply(rbind(shape.gev,sapply(shape.gev,optim.gev.sort)),2,z.gev.sort),silent=TRUE)
    }
    Z.GUM.sort <- qnorm(Repartition.EV(sort(X)))
    w <- qnorm(Quantiles)
    tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
    tn_ppcc.Gev <- rep(0,length(shape.gev))
    for(l in 1:length(shape.gev)) {
      numer <- sum((Z.GEV.sort[,l]-mean(Z.GEV.sort[,l]))*(w-mean(w)))
      denom <- sqrt( sum(((Z.GEV.sort[,l]-mean(Z.GEV.sort[,l]))^2))* sum(((w-mean(w))^2))) 
      tn_ppcc.Gev[l] <- numer/denom
    }
    tn_ppcc.Gev <- ifelse(is.nan(tn_ppcc.Gev),0,tn_ppcc.Gev)
    for(l in 1:length(shape.gev)) if(tn_ppcc.Gumbel>tn_ppcc.Gev[l]){res.GUM.GEV[i,l] <- 1}else{res.GUM.GEV[i,l] <- 0}
  }
  for(l in 1:length(shape.gev)) PCS.TN_PPCC.GUM.GEV[l,j] <- sum(res.GUM.GEV[,l])/Nombre.Iteration
  j <- j+1
}

PCS.TN_PPCC.GUM.GEV <- round(PCS.TN_PPCC.GUM.GEV,2)*100

###---PCS(?chantillon GEV, alternatif Gumbel) avec la m?thode tn.ppcc---###
m_gev_tn.ppcc <- array(0,c(length(shape.gev),length(Taille.Echantillon)),dimnames=list(shape.gev,Taille.Echantillon))
res.GEV.GUM <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.GEV.GUM <- rep(0,length(shape.gev))
tn_ppcc.Gev.tab <- rep(0,Nombre.Iteration)
tn_ppcc.Gumbel.tab <- rep(0,Nombre.Iteration)
set.seed(seed)
r<-1
for(n in Taille.Echantillon){
  j <- 1
  Quantiles <- (1:n)/(n+1)
  for(a in shape.gev){
    for(i in 1:Nombre.Iteration){
      Z.GEV.sort <- NULL
      while(is.null(Z.GEV.sort)){ 
        X <- rgev(n,a,1,0)
        param.GEV <- optim(c(1.12,0.07),log.GEV,sh=a,X=sort(X))$par
        try(Z.GEV.sort <- qnorm(pgev(sort(X),a,param.GEV[1],param.GEV[2])),silent=T) 
      }
      param.GUM <- optim(c(1,0),log.EV,X=sort(X))$par
      Z.GUM.sort <- qnorm(pgumbel(sort(X),param.GUM[1],param.GUM[2]))
      w <- qnorm(Quantiles)
      tn_ppcc.Gumbel.tab[i] <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
      tn_ppcc.Gev.tab[i] <- sum((Z.GEV.sort-mean(Z.GEV.sort))*(w-mean(w)))/sqrt(sum(((Z.GEV.sort-mean(Z.GEV.sort))^2))*sum(((w-mean(w))^2)))
    }
    tn_ppcc.Gev.tab <- ifelse(is.nan(tn_ppcc.Gev.tab),0,tn_ppcc.Gev.tab)
    res.GEV.GUM <- ifelse(tn_ppcc.Gev.tab>tn_ppcc.Gumbel.tab,1,0)
    PCS.TN_PPCC.GEV.GUM[j] <- sum(res.GEV.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_gev_tn.ppcc[,r]<-PCS.TN_PPCC.GEV.GUM
  r<-r+1
}

m_gev_tn.ppcc <- round(m_gev_tn.ppcc,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(shape.gev)){
  print(cbind(PCS.TN_PPCC.GUM.GEV[l,],m_gev_tn.ppcc[l,]))
}


###############################################################################################################################################

###------Discrimination between the Gumbel and the Shifted Student t (T+m) distributions----------###

############-----------M?thode RML--------#################
###---PCS(?chantillon Gumbel, alternatif student) avec la m?thode RML---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.STU <- rep(0,Nombre.Iteration)
PCS.RML.GUM.STU <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    LLGUM <- -optim(par=c(1,0),log.EV,X=X)$value
    LLSTU<- -optim(par=c(2,1),log.STU.DCL,X=X)$value
    if(LLGUM-LLSTU>0){res.GUM.STU[i] <- 1}else{res.GUM.STU[i] <- 0}
  }
  PCS.RML.GUM.STU[j] <- sum(res.GUM.STU)/Nombre.Iteration
  j <- j+1
}

PCS.RML.GUM.STU <- round(PCS.RML.GUM.STU,2)*100

###---PCS(?chantillon student, alternatif Gumbel) avec la m?thode RML---###
m_student_rml <- array(0,c(length(dl.student),length(Taille.Echantillon)),dimnames=list(dl.student,Taille.Echantillon))
res.STU.GUM <- rep(0,Nombre.Iteration)
PCS.RML.STU.GUM <- rep(0,length(dl.student))
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  for(k in dl.student){
    for(i in 1:Nombre.Iteration){
      LLGUM <- NULL
      while(is.null(LLGUM)){
        X = rnorm(n)/sqrt(rchisq(n,k)/k)
        try(LLGUM <- -optim(par=c(0,1),fn=Vraisemblance.EV,X=X)$value,silent=TRUE)
      }
      LLSTU <- -optim(par=c(2,1),log.STU.DCL,X=X)$value
      if(LLGUM-LLSTU<0){res.STU.GUM[i] <- 1}else{res.STU.GUM[i] <- 0}
    }
    PCS.RML.STU.GUM[j] <- sum(res.STU.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_student_rml[,r] <- PCS.RML.STU.GUM
  r <- r+1
}

m_student_rml <- round(m_student_rml,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(dl.student)){
  print(cbind(PCS.RML.GUM.STU,m_student_rml[l,]))
}

###########-----------M?thode TNSW-----------###############
###---PCS(?chantillon Gumbel, alternatif student) avec la m?thode TN.SW---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.STU <- rep(0,Nombre.Iteration)
PCS.TNSW.GUM.STU <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    # Calcul de la statistique du test pour la loi Gumbel
    Z.GUM <- qnorm(Repartition.EV(X))
    S.GUM <- shapiro.test(Z.GUM)$statistic
    # Calcul de la statistique du test pour la loi student
    param.stu <- optim(c(2,1),f=log.STU.DCL,X=X)$par
    Z.STU <- qnorm(pt(X-param.stu[2],df=param.stu[1],ncp=0))
    S.STU <- shapiro.test(Z.STU)$statistic
    if(S.STU<S.GUM){res.GUM.STU[i] <- 1}else{res.GUM.STU[i] <- 0}
  }
  PCS.TNSW.GUM.STU[j] <- sum(res.GUM.STU)/Nombre.Iteration
  j <- j+1
}

PCS.TNSW.GUM.STU <- round(PCS.TNSW.GUM.STU,2)*100

###---PCS(?chantillon student, alternatif Gumbel) avec la m?thode TN.SW---###
m_student_tnsw <- array(0,c(length(dl.student),length(Taille.Echantillon)),dimnames=list(dl.student,Taille.Echantillon))
res.STU.GUM <- rep(0,Nombre.Iteration)
PCS.TNSW.STU.GUM <- rep(0,length(dl.student))
Val.Stat.Gum <- rep(0,Nombre.Iteration)
Val.Stat.Stu <- rep(0,Nombre.Iteration)
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  for(k in dl.student){
    for(i in 1:Nombre.Iteration){
      # Calcul de la statistique de test pour la loi Gumbel
      S.GUM <- NULL
      while(is.null(S.GUM)){
        X = rnorm(n)/sqrt(rchisq(n,k)/k)
        try(S.GUM <- shapiro.test(qnorm(Repartition.EV(X)))$statistic,silent=TRUE)
      }
      # Calcul de la statistique de test pour la loi student
      param.stu <- optim(c(2,1),log.STU.DCL,X=X)$par
      Z.STU <- qnorm(pt(X-param.stu[2],df=param.stu[1],ncp=0))
      S.STU <- shapiro.test(Z.STU)$statistic
      Val.Stat.Gum[i] <- S.GUM
      Val.Stat.Stu[i] <- S.STU
    }
    Val.Stat.Gum <- ifelse(is.na(Val.Stat.Gum),0,Val.Stat.Gum)
    res.STU.GUM <- ifelse(Val.Stat.Gum<Val.Stat.Stu,1,0)
    PCS.TNSW.STU.GUM[j] <- sum(res.STU.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_student_tnsw[,r] <- PCS.TNSW.STU.GUM
  r <- r+1
}

m_student_tnsw <- round(m_student_tnsw,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(dl.student)){
  print(cbind(PCS.TNSW.GUM.STU,m_student_tnsw[l,]))
}

#########---------------M?thode Anderson-darling-----------##########
###---PCS(?chantillon Gumbel, alternatif student) avec la m?thode AD---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.STU <- rep(0,Nombre.Iteration)
PCS.AD.GUM.STU <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    AGUM <- ad.test(X,Repartition.EV)$statistic
    param.stu <- optim(c(2,1),log.STU.DCL,X=X)$par
    ASTU <- ad.test(X-param.stu[2],pt,df=param.stu[1],ncp=0)$statistic
    if(ASTU>AGUM){res.GUM.STU[i] <- 1}else{res.GUM.STU[i] <- 0}
  }
  PCS.AD.GUM.STU[j] <- sum(res.GUM.STU)/Nombre.Iteration
  j <- j+1
}

PCS.AD.GUM.STU <- round(PCS.AD.GUM.STU,2)*100

###---PCS(?chantillon student, alternatif Gumbel) avec la m?thode AD---###
m_student_ad <- array(0,c(length(dl.student),length(Taille.Echantillon)),dimnames=list(dl.student,Taille.Echantillon)) 
res.STU.GUM <- rep(0,Nombre.Iteration)
PCS.AD.STU.GUM <- rep(0,length(dl.student))
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  for(k in dl.student){
    for(i in 1:Nombre.Iteration){
      AGUM <- NULL
      while(is.null(AGUM)){
        X = rnorm(n)/sqrt(rchisq(n,k)/k)
        try(AGUM <- ad.test(X,Repartition.EV)$statistic,silent=TRUE)
      }
      param.stu <- optim(c(2,1),log.STU.DCL,X=X)$par
      ASTU <-ad.test(X-param.stu[2],pt,df=param.stu[1],ncp=0)$statistic
      if(AGUM>ASTU){res.STU.GUM[i] <- 1}else{res.STU.GUM[i] <- 0}
    }
    PCS.AD.STU.GUM[j] <- sum(res.STU.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_student_ad[,r] <- PCS.AD.STU.GUM
  r <- r+1
}

m_student_ad <- round(m_student_ad,2)*100

###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(dl.student)){
  print(cbind(PCS.AD.GUM.STU,m_student_ad[l,]))
}

#########---------------M?thode TN.PPCC-----------##########
###---PCS(?chantillon Gumbel, alternatif student) avec la m?thode tn.ppcc---###
j <- 1
Nombre.Iteration <- NumberIterations
res.GUM.STU <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.GUM.STU <- rep(0,length(Taille.Echantillon))
set.seed(seed)
for(n in Taille.Echantillon){
  Quantiles <- (1:n)/(n+1)
  for(i in 1:Nombre.Iteration){
    X <- rgumbel(n,1,0)
    Z.GUM.sort <- qnorm(Repartition.EV(sort(X)))
    param.stu <- optim(c(2,1),f=log.STU.DCL,X=sort(X))$par
    Z.STU.sort <- qnorm(pt(sort(X-param.stu[2]),df=param.stu[1],ncp=0))
    w <- qnorm(Quantiles)
    tn_ppcc.Gumbel <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
    tn_ppcc.Student <- sum((Z.STU.sort-mean(Z.STU.sort))*(w-mean(w)))/sqrt(sum(((Z.STU.sort-mean(Z.STU.sort))^2))*sum(((w-mean(w))^2)))
    if(tn_ppcc.Gumbel>tn_ppcc.Student){res.GUM.STU[i] <- 1}else{res.GUM.STU[i] <- 0}
  }
  PCS.TN_PPCC.GUM.STU[j] <- sum(res.GUM.STU)/Nombre.Iteration
  j <- j+1
}

PCS.TN_PPCC.GUM.STU <- round(PCS.TN_PPCC.GUM.STU,2)*100

###---PCS(?chantillon student, alternatif Gumbel) avec la m?thode tn.ppcc---###
m_student_tn.ppcc <- array(0,c(length(dl.student),length(Taille.Echantillon)),dimnames=list(dl.student,Taille.Echantillon)) 
res.STU.GUM <- rep(0,Nombre.Iteration)
PCS.TN_PPCC.STU.GUM <- rep(0,length(dl.student))
tn_ppcc.Gumbel.tab <- rep(0,Nombre.Iteration)
tn_ppcc.Student.tab <- rep(0,Nombre.Iteration)
set.seed(seed)
r <- 1
for(n in Taille.Echantillon){
  j <- 1
  Quantiles <- (1:n)/(n+1)
  for(k in dl.student){
    for(i in 1:Nombre.Iteration){
      Z.GUM.sort <- NULL
      while(is.null(Z.GUM.sort)){
        X = rnorm(n)/sqrt(rchisq(n,k)/k)
        try(Z.GUM.sort <- qnorm(Repartition.EV(sort(X))),silent=TRUE)
      }
      param.stu <- optim(c(2,1),log.STU.DCL,X=sort(X))$par
      Z.STU.sort <- qnorm(pt(sort(X-param.stu[2]),df=param.stu[1],ncp=0))
      w <- qnorm(Quantiles)
      tn_ppcc.Gumbel.tab[i] <- sum((Z.GUM.sort-mean(Z.GUM.sort))*(w-mean(w)))/sqrt(sum(((Z.GUM.sort-mean(Z.GUM.sort))^2))*sum(((w-mean(w))^2)))
      tn_ppcc.Student.tab[i] <- sum((Z.STU.sort-mean(Z.STU.sort))*(w-mean(w)))/sqrt(sum(((Z.STU.sort-mean(Z.STU.sort))^2))*sum(((w-mean(w))^2)))
    }
    tn_ppcc.Gumbel.tab <- ifelse(is.nan(tn_ppcc.Gumbel.tab),0,tn_ppcc.Gumbel.tab)
    res.STU.GUM <- ifelse(tn_ppcc.Student.tab>tn_ppcc.Gumbel.tab,1,0)
    PCS.TN_PPCC.STU.GUM[j] <- sum(res.STU.GUM)/Nombre.Iteration
    j <- j+1
  }
  m_student_tn.ppcc[,r] <- PCS.TN_PPCC.STU.GUM
  r <- r+1
}

m_student_tn.ppcc <- round(m_student_tn.ppcc,2)*100


###----Affiche les r?sultats par param?tre-----###
for(l in 1:length(dl.student)){
  print(cbind(PCS.TN_PPCC.GUM.STU,m_student_tn.ppcc[l,]))
}

##############################################################################################################################################################

#############--------------------Affichage r?sultat--------------------------###############

###---------------Gumbel et Normal-----------------###
AD.GUM.NOR <- matrix(paste(PCS.AD.GUM.NOR,paste0(PCS.AD.NOR.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
TNSW.GUM.NOR <- matrix(paste(PCS.TNSW.GUM.NOR,paste0(PCS.TNSW.NOR.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
RML.GUM.NOR <- matrix(paste(PCS.RML.GUM.NOR,paste0(PCS.RML.NOR.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
TN_PPCC.GUM.NOR <- matrix(paste(PCS.TN_PPCC.GUM.NOR,paste0(PCS.TN_PPCC.NOR.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
PCS.GUM.NOR <- cbind(RML.GUM.NOR,TNSW.GUM.NOR,AD.GUM.NOR,TN_PPCC.GUM.NOR)

###-------------Gumbel et Logistique----------###
AD.GUM.LOG <- matrix(paste(PCS.AD.GUM.LOG,paste0(PCS.AD.LOG.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
TNSW.GUM.LOG <- matrix(paste(PCS.TNSW.GUM.LOG,paste0(PCS.TNSW.LOG.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
RML.GUM.LOG <- matrix(paste(PCS.RML.GUM.LOG,paste0(PCS.RML.LOG.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
TN_PPCC.GUM.LOG <- matrix(paste(PCS.TN_PPCC.GUM.LOG,paste0(PCS.TN_PPCC.LOG.GUM,rep(")")),sep='('),length(Taille.Echantillon),1)
PCS.GUM.LOG <- cbind(RML.GUM.LOG,TNSW.GUM.LOG,AD.GUM.LOG,TN_PPCC.GUM.LOG)



###---------Gumbel et Pearson type 3-----------###
L_GUM_GAM3 = list()
for(i in 1:length(shape.gamma3)){
  L_GUM_GAM3[[i]] <- cbind(paste(PCS.RML.GUM.GAM3[i,],paste0(m_gam3_rml[i,],rep(")")),sep='('),paste(PCS.TNSW.GUM.GAM3[i,],paste0(m_gam3_tnsw[i,],rep(")")),sep='('),paste(PCS.AD.GUM.GAM3[i,],paste0(m_gam3_ad[i,],rep(")")),sep='('),paste(PCS.TN_PPCC.GUM.GAM3[i,],paste0(m_gam3_tn.ppcc[i,],rep(")")),sep='('))
}

###---------Gumbel et Gev-----------###
L_GUM_GEV = list()
for(i in 1:length(shape.gev)){
  L_GUM_GEV[[i]] <- cbind(paste(PCS.RML.GUM.GEV[i,],paste0(m_gev_rml[i,],rep(")")),sep='('),paste(PCS.TNSW.GUM.GEV[i,],paste0(m_gev_tnsw[i,],rep(")")),sep='('),paste(PCS.AD.GUM.GEV[i,],paste0(m_gev_ad[i,],rep(")")),sep='('),paste(PCS.TN_PPCC.GUM.GEV[i,],paste0(m_gev_tn.ppcc[i,],rep(")")),sep='('))
}

###--------Gumbel et student--------###
L_GUM_STU = list()
for(i in 1:length(dl.student)){
  L_GUM_STU[[i]] <- cbind(paste(PCS.RML.GUM.STU,paste0(m_student_rml[i,],rep(")")),sep='('),paste(PCS.TNSW.GUM.STU,paste0(m_student_tnsw[i,],rep(")")),sep='('),paste(PCS.AD.GUM.STU,paste0(m_student_ad[i,],rep(")")),sep='('),paste(PCS.TN_PPCC.GUM.STU,paste0(m_student_tn.ppcc[i,],rep(")")),sep='('))
}

###############################################################################################################################
##########------------------Exportation des PCS vers EXCEL-------------##################
Used.statistics <- c("RML","TNSW","AD","TN.PPCC")
Col.GUM.NOR <- rep("GUM(NOR)",length(Used.statistics))
Col.GUM.LOG <- rep("GUM(LOG)",length(Used.statistics))
PCS.GN <- rbind(Col.GUM.NOR,Used.statistics,PCS.GUM.NOR)
PCS.GL <- rbind(Col.GUM.LOG,Used.statistics,PCS.GUM.LOG)
PCS.GNL <- cbind(PCS.GN,PCS.GL)

###---Nom colonne pour Gumbel et Gamma3---###
Gam3_para = rep("a",length(shape.gamma3))
for(i in 1:length(shape.gamma3)){
  Gam3_para[i] <- paste("GAM",shape.gamma3[i],sep="_")
}

Gum_Gam3 = matrix(0,length(shape.gamma3),length(Used.statistics))
for(i in 1:length(shape.gamma3)){
  Gum_Gam3[i,] <- matrix(rep(paste("GUM",paste0(Gam3_para[i],rep(")")),sep='('),length(Used.statistics)),1,length(Used.statistics))
}
Gum_Gam3 = matrix(t(Gum_Gam3),1,length(shape.gamma3)*length(Used.statistics))
Used.statistics.gam3 <- rbind(Gum_Gam3,rep(Used.statistics,length(shape.gamma3)))
PCS.GUM.GAM3 <- rbindlist(lapply(list(Used.statistics.gam3,L_GUM_GAM3),as.data.frame))

###---Nom colonne pour Gumbel et Gev---###
Gev_para = rep("a",length(shape.gev))
for(i in 1:length(shape.gev)){
  Gev_para[i] <- paste("GEV",shape.gev[i],sep="_")
}

Gum_Gev = matrix(0,length(shape.gev),length(Used.statistics))
for(i in 1:length(shape.gev)){
  Gum_Gev[i,] <- matrix(rep(paste("GUM",paste0(Gev_para[i],rep(")")),sep='('),length(Used.statistics)),1,length(Used.statistics))
}
Gum_Gev = matrix(t(Gum_Gev),1,length(shape.gev)*length(Used.statistics))
Used.statistics.gev <- rbind(Gum_Gev,rep(Used.statistics,length(shape.gev)))
PCS.GUM.GEV <- rbindlist(lapply(list(Used.statistics.gev,L_GUM_GEV),as.data.frame))

###---Nom colonne pour Gumbel et Student---###
Stu_para = rep("a",length(dl.student))
for(i in 1:length(dl.student)){
  Stu_para[i] <- paste("STU",dl.student[i],sep="_")
}

Gum_Stu = matrix(0,length(dl.student),length(Used.statistics))
for(i in 1:length(dl.student)){
  Gum_Stu[i,] <- matrix(rep(paste("GUM",paste0(Stu_para[i],rep(")")),sep='('),length(Used.statistics)),1,length(Used.statistics))
}
Gum_Stu = matrix(t(Gum_Stu),1,length(dl.student)*length(Used.statistics))
Used.statistics.stu <- rbind(Gum_Stu,rep(Used.statistics,length(dl.student)))
PCS.GUM.STU <- rbindlist(lapply(list(Used.statistics.stu,L_GUM_STU),as.data.frame))

PCS.Result <- cbind(PCS.GNL,PCS.GUM.GAM3,PCS.GUM.GEV,PCS.GUM.STU)
PCS.Result <- data.frame(PCS.Result,row.names=Sample.size)
write.table(PCS.Result,"PCS.Result.csv",quote=FALSE,row.names=T,col.names=NA)

#######################################################################################################################

#######-------------Moyennes et Diff?rences-----------------------##########

###----------------------Gumbel et Normal------------------------###
######-----------------Anderson-Darling-------------------######
PCS.MEAN.AD.GUM.NOR <- round((GUM.NOR.AD[,2]+GUM.NOR.AD[,3])/2,2)
PCS.DIFF.AD.GUM.NOR <- round(abs(GUM.NOR.AD[,2]-GUM.NOR.AD[,3]),2)
MEAN.DIFF.AD.GUM.NOR <- matrix(paste(PCS.MEAN.AD.GUM.NOR,paste0(PCS.DIFF.AD.GUM.NOR,rep(")")),sep='('),length(Taille.Echantillon),1)
######-----------------TNSW-------------------######
PCS.MEAN.TNSW.GUM.NOR <- round((GUM.NOR.TNSW[,2]+GUM.NOR.TNSW[,3])/2,2)
PCS.DIFF.TNSW.GUM.NOR <- round(abs(GUM.NOR.TNSW[,2]-GUM.NOR.TNSW[,3]),2)
MEAN.DIFF.TNSW.GUM.NOR <- matrix(paste(PCS.MEAN.TNSW.GUM.NOR,paste0(PCS.DIFF.TNSW.GUM.NOR,rep(")")),sep='('),length(Taille.Echantillon),1)
######-----------------RML-------------------######
PCS.MEAN.RML.GUM.NOR <- round((GUM.NOR.RML[,2]+GUM.NOR.RML[,3])/2,2)
PCS.DIFF.RML.GUM.NOR <- round(abs(GUM.NOR.RML[,2]-GUM.NOR.RML[,3]),2)
MEAN.DIFF.RML.GUM.NOR <- matrix(paste(PCS.MEAN.RML.GUM.NOR,paste0(PCS.DIFF.RML.GUM.NOR,rep(")")),sep='('),length(Taille.Echantillon),1)
######-----------------TN.PPCC-------------------######
PCS.MEAN.TN_PPCC.GUM.NOR <- round((GUM.NOR.TN_PPCC[,2]+GUM.NOR.TN_PPCC[,3])/2,2)
PCS.DIFF.TN_PPCC.GUM.NOR <- round(abs(GUM.NOR.TN_PPCC[,2]-GUM.NOR.TN_PPCC[,3]),2)
MEAN.DIFF.TN_PPCC.GUM.NOR <- matrix(paste(PCS.MEAN.TN_PPCC.GUM.NOR,paste0(PCS.DIFF.TN_PPCC.GUM.NOR,rep(")")),sep='('),length(Taille.Echantillon),1)
########---------Concat?nation des moyennes et diff?rences de PCS------------###########
MEAN.DIFF.GUM.NOR <- cbind(MEAN.DIFF.RML.GUM.NOR,MEAN.DIFF.TNSW.GUM.NOR,MEAN.DIFF.AD.GUM.NOR,MEAN.DIFF.TN_PPCC.GUM.NOR)

###----------------------Gumbel et Logistique------------------------###
######-----------------Anderson-Darling-------------------######
PCS.MEAN.AD.GUM.LOG <- round((GUM.LOG.AD[,2]+GUM.LOG.AD[,3])/2,2)
PCS.DIFF.AD.GUM.LOG <- round(abs(GUM.LOG.AD[,2]-GUM.LOG.AD[,3]),2)
MEAN.DIFF.AD.GUM.LOG <- matrix(paste(PCS.MEAN.AD.GUM.LOG,paste0(PCS.DIFF.AD.GUM.LOG,rep(")")),sep='('),length(Taille.Echantillon),1)
######-----------------TNSW-------------------######
PCS.MEAN.TNSW.GUM.LOG <- round((GUM.LOG.TNSW[,2]+GUM.LOG.TNSW[,3])/2,2)
PCS.DIFF.TNSW.GUM.LOG <- round(abs(GUM.LOG.TNSW[,2]-GUM.LOG.TNSW[,3]),2)
MEAN.DIFF.TNSW.GUM.LOG <- matrix(paste(PCS.MEAN.TNSW.GUM.LOG,paste0(PCS.DIFF.TNSW.GUM.LOG,rep(")")),sep='('),length(Taille.Echantillon),1)
######-----------------RML-------------------######
PCS.MEAN.RML.GUM.LOG <- round((GUM.LOG.RML[,2]+GUM.LOG.RML[,3])/2,2)
PCS.DIFF.RML.GUM.LOG <- round(abs(GUM.LOG.RML[,2]-GUM.LOG.RML[,3]),2)
MEAN.DIFF.RML.GUM.LOG <- matrix(paste(PCS.MEAN.RML.GUM.LOG,paste0(PCS.DIFF.RML.GUM.LOG,rep(")")),sep='('),length(Taille.Echantillon),1)
######-----------------TN.PPCC-------------------######
PCS.MEAN.TN_PPCC.GUM.LOG <- round((GUM.LOG.TN_PPCC[,2]+GUM.LOG.TN_PPCC[,3])/2,2)
PCS.DIFF.TN_PPCC.GUM.LOG <- round(abs(GUM.LOG.TN_PPCC[,2]-GUM.LOG.TN_PPCC[,3]),2)
MEAN.DIFF.TN_PPCC.GUM.LOG <- matrix(paste(PCS.MEAN.TN_PPCC.GUM.LOG,paste0(PCS.DIFF.TN_PPCC.GUM.LOG,rep(")")),sep='('),length(Taille.Echantillon),1)
########---------Concat?nation des moyennes et diff?rences de PCS------------###########
MEAN.DIFF.GUM.LOG <- cbind(MEAN.DIFF.RML.GUM.LOG,MEAN.DIFF.TNSW.GUM.LOG,MEAN.DIFF.AD.GUM.LOG,MEAN.DIFF.TN_PPCC.GUM.LOG)


###----------------------Gumbel et Pearson type 3------------------------##########
#############-------------Anderson-Darling--------------################
PCS.MEAN.AD.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
PCS.DIFF.AD.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
for(l in 1:length(shape.gamma3)){ 
  PCS.MEAN.AD.GUM.GAM3[l,] <- round((PCS.AD.GUM.GAM3[l,]+m_gam3_ad[l,])/2,2)
  PCS.DIFF.AD.GUM.GAM3[l,] <- round(abs(PCS.AD.GUM.GAM3[l,]-m_gam3_ad[l,]),2)
}
#############-------------TNSW--------------################
PCS.MEAN.TNSW.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
PCS.DIFF.TNSW.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
for(l in 1:length(shape.gamma3)){ 
  PCS.MEAN.TNSW.GUM.GAM3[l,] <- round((PCS.TNSW.GUM.GAM3[l,]+m_gam3_tnsw[l,])/2,2)
  PCS.DIFF.TNSW.GUM.GAM3[l,] <- round(abs(PCS.TNSW.GUM.GAM3[l,]-m_gam3_tnsw[l,]),2)
}
#############-------------RML--------------################
PCS.MEAN.RML.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
PCS.DIFF.RML.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
for(l in 1:length(shape.gamma3)){ 
  PCS.MEAN.RML.GUM.GAM3[l,] <- round((PCS.RML.GUM.GAM3[l,]+m_gam3_rml[l,])/2,2)
  PCS.DIFF.RML.GUM.GAM3[l,] <- round(abs(PCS.RML.GUM.GAM3[l,]-m_gam3_rml[l,]),2)
}
#############-------------TN.PPCC--------------################
PCS.MEAN.TN_PPCC.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
PCS.DIFF.TN_PPCC.GUM.GAM3 <- matrix(0,length(shape.gamma3),length(Taille.Echantillon))
for(l in 1:length(shape.gamma3)){ 
  PCS.MEAN.TN_PPCC.GUM.GAM3[l,] <- round((PCS.TN_PPCC.GUM.GAM3[l,]+m_gam3_tn.ppcc[l,])/2,2)
  PCS.DIFF.TN_PPCC.GUM.GAM3[l,] <- round(abs(PCS.TN_PPCC.GUM.GAM3[l,]-m_gam3_tn.ppcc[l,]),2)
}

########-------Concat?nation des moyennes et diff?rences de PCS---------############
MEAN_DIFF_GUM_GAM3 = list()
for(i in 1:length(shape.gamma3)){
  MEAN_DIFF_GUM_GAM3[[i]] <- cbind(paste(PCS.MEAN.RML.GUM.GAM3[i,],paste0(PCS.DIFF.RML.GUM.GAM3[i,],rep(")")),sep='('),paste(PCS.MEAN.TNSW.GUM.GAM3[i,],paste0(PCS.DIFF.TNSW.GUM.GAM3[i,],rep(")")),sep='('),paste(PCS.MEAN.AD.GUM.GAM3[i,],paste0(PCS.DIFF.AD.GUM.GAM3[i,],rep(")")),sep='('),paste(PCS.MEAN.TN_PPCC.GUM.GAM3[i,],paste0(PCS.DIFF.TN_PPCC.GUM.GAM3[i,],rep(")")),sep='('))
}

###----------------------Gumbel et Gev------------------------##########
#############-------------Anderson-Darling--------------################
PCS.MEAN.AD.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
PCS.DIFF.AD.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
for(l in 1:length(shape.gev)){ 
  PCS.MEAN.AD.GUM.GEV[l,] <- round((PCS.AD.GUM.GEV[l,]+m_gev_ad[l,])/2,2)
  PCS.DIFF.AD.GUM.GEV[l,] <- round(abs(PCS.AD.GUM.GEV[l,]-m_gev_ad[l,]),2)
}
#############-------------TNSW--------------################
PCS.MEAN.TNSW.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
PCS.DIFF.TNSW.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
for(l in 1:length(shape.gev)){ 
  PCS.MEAN.TNSW.GUM.GEV[l,] <- round((PCS.TNSW.GUM.GEV[l,]+m_gev_tnsw[l,])/2,2)
  PCS.DIFF.TNSW.GUM.GEV[l,] <- round(abs(PCS.TNSW.GUM.GEV[l,]-m_gev_tnsw[l,]),2)
}
#############-------------RML--------------################
PCS.MEAN.RML.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
PCS.DIFF.RML.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
for(l in 1:length(shape.gev)){ 
  PCS.MEAN.RML.GUM.GEV[l,] <- round((PCS.RML.GUM.GEV[l,]+m_gev_rml[l,])/2,2)
  PCS.DIFF.RML.GUM.GEV[l,] <- round(abs(PCS.RML.GUM.GEV[l,]-m_gev_rml[l,]),2)
}
#############-------------TN.PPCC--------------################
PCS.MEAN.TN_PPCC.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
PCS.DIFF.TN_PPCC.GUM.GEV <- matrix(0,length(shape.gev),length(Taille.Echantillon))
for(l in 1:length(shape.gev)){ 
  PCS.MEAN.TN_PPCC.GUM.GEV[l,] <- round((PCS.TN_PPCC.GUM.GEV[l,]+m_gev_tn.ppcc[l,])/2,2)
  PCS.DIFF.TN_PPCC.GUM.GEV[l,] <- round(abs(PCS.TN_PPCC.GUM.GEV[l,]-m_gev_tn.ppcc[l,]),2)
}

########-------Concat?nation des moyennes et diff?rences de PCS---------############
MEAN_DIFF_GUM_GEV = list()
for(i in 1:length(shape.gev)){
  MEAN_DIFF_GUM_GEV[[i]] <- cbind(paste(PCS.MEAN.RML.GUM.GEV[i,],paste0(PCS.DIFF.RML.GUM.GEV[i,],rep(")")),sep='('),paste(PCS.MEAN.TNSW.GUM.GEV[i,],paste0(PCS.DIFF.TNSW.GUM.GEV[i,],rep(")")),sep='('),paste(PCS.MEAN.AD.GUM.GEV[i,],paste0(PCS.DIFF.AD.GUM.GEV[i,],rep(")")),sep='('),paste(PCS.MEAN.TN_PPCC.GUM.GEV[i,],paste0(PCS.DIFF.TN_PPCC.GUM.GEV[i,],rep(")")),sep='('))
}

###----------------------Gumbel et Student------------------------##########
#############-------------Anderson-Darling--------------################
PCS.MEAN.AD.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
PCS.DIFF.AD.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
for(l in 1:length(dl.student)){ 
  PCS.MEAN.AD.GUM.STU[l,] <- round((PCS.AD.GUM.STU+m_student_ad[l,])/2,2)
  PCS.DIFF.AD.GUM.STU[l,] <- round(abs(PCS.AD.GUM.STU-m_student_ad[l,]),2)
}
#############-------------TNSW--------------################
PCS.MEAN.TNSW.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
PCS.DIFF.TNSW.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
for(l in 1:length(dl.student)){ 
  PCS.MEAN.TNSW.GUM.STU[l,] <- round((PCS.TNSW.GUM.STU+m_student_tnsw[l,])/2,2)
  PCS.DIFF.TNSW.GUM.STU[l,] <- round(abs(PCS.TNSW.GUM.STU-m_student_tnsw[l,]),2)
}
#############-------------RML--------------################
PCS.MEAN.RML.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
PCS.DIFF.RML.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
for(l in 1:length(dl.student)){ 
  PCS.MEAN.RML.GUM.STU[l,] <- round((PCS.RML.GUM.STU+m_student_rml[l,])/2,2)
  PCS.DIFF.RML.GUM.STU[l,] <- round(abs(PCS.RML.GUM.STU-m_student_rml[l,]),2)
}
#############-------------TN.PPCC--------------################
PCS.MEAN.TN_PPCC.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
PCS.DIFF.TN_PPCC.GUM.STU <- matrix(0,length(dl.student),length(Taille.Echantillon))
for(l in 1:length(dl.student)){ 
  PCS.MEAN.TN_PPCC.GUM.STU[l,] <- round((PCS.TN_PPCC.GUM.STU+m_student_tn.ppcc[l,])/2,2)
  PCS.DIFF.TN_PPCC.GUM.STU[l,] <- round(abs(PCS.TN_PPCC.GUM.STU-m_student_tn.ppcc[l,]),2)
}

########-------Concat?nation des moyennes et diff?rences de PCS---------############
MEAN_DIFF_GUM_STU = list()
for(i in 1:length(dl.student)){
  MEAN_DIFF_GUM_STU[[i]] <- cbind(paste(PCS.MEAN.RML.GUM.STU[i,],paste0(PCS.DIFF.RML.GUM.STU[i,],rep(")")),sep='('),paste(PCS.MEAN.TNSW.GUM.STU[i,],paste0(PCS.DIFF.TNSW.GUM.STU[i,],rep(")")),sep='('),paste(PCS.MEAN.AD.GUM.STU[i,],paste0(PCS.DIFF.AD.GUM.STU[i,],rep(")")),sep='('),paste(PCS.MEAN.TN_PPCC.GUM.STU[i,],paste0(PCS.DIFF.TN_PPCC.GUM.STU[i,],rep(")")),sep='('))
}

########################---Exportation des moyennes et diff?rences de PCS vers Excel---#########################
MEAN.DIFF.PCS.GUM.NOR <- rbind(Col.GUM.NOR,Used.statistics,MEAN.DIFF.GUM.NOR)
MEAN.DIFF.PCS.GUM.LOG <- rbind(Col.GUM.LOG,Used.statistics,MEAN.DIFF.GUM.LOG)
MEAN.DIFF.PCS.GUM.GAM3 <- rbindlist(lapply(list(Used.statistics.gam3,MEAN_DIFF_GUM_GAM3),as.data.frame))
MEAN.DIFF.PCS.GUM.GEV <- rbindlist(lapply(list(Used.statistics.gev,MEAN_DIFF_GUM_GEV),as.data.frame))
MEAN.DIFF.PCS.GUM.STU <- rbindlist(lapply(list(Used.statistics.stu,MEAN_DIFF_GUM_STU),as.data.frame))
PCS.MEAN.DIFF.Result <- cbind(MEAN.DIFF.PCS.GUM.NOR,MEAN.DIFF.PCS.GUM.LOG,MEAN.DIFF.PCS.GUM.GAM3,MEAN.DIFF.PCS.GUM.GEV,MEAN.DIFF.PCS.GUM.STU)
PCS.MEAN.DIFF.Result <- data.frame(PCS.MEAN.DIFF.Result,row.names=Sample.size)
write.table(PCS.MEAN.DIFF.Result,"PCS.Mean.Diff.Result.csv",quote=FALSE,row.names=T,col.names=NA)
