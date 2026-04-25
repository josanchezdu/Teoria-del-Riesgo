

# Jorge Andres Sanchez Duarte
#Numero de reclamos, el 6 identifica los mayores a 6
library(polynom)
k<-c(0:6)
n_k<-c(25356,1521,282,58,16,4,1) # simplemente, cambie los datos por la cartera que desee
Tcasos<-sum(n_k)
# Frecuencia relativa
freqrel<-(n_k/Tcasos)
# Estimadores
n_bar <- sum(k*n_k)/Tcasos
s_cua <- (sum(k^2*n_k)/Tcasos)-(n_bar^2)
n_bar;s_cua

### Prueba Chicuadrado
Chi_cua<-function(obs, esp,alpha=0.05,m=7,gl=3){
  chi<-sum((obs-esp)^2/esp)
  vchi<-qchisq(1-alpha,df = gl )
  if(chi<vchi){
    return(list("No hay evidencia para rechazar H_0",
                "Estadistica"=chi,"Valor critico"=vchi))
  }
  if(chi>=vchi){
    return(list("Hay evidencia para rechazar H_0",
                "Estadistica"=chi,"Valor critico"=vchi))
  }
}


# Momentos ------------------------------------------
m_4=0
for(k in 4:(length(n_k)-1)){
  m_4<-((gamma(k+1)/(gamma(k-3)))*freqrel[k+1])+m_4
}
m_3=0
for(k in 3:(length(n_k)-1)){
  m_3<-((gamma(k+1)/(gamma(k-2)))*freqrel[k+1])+m_3
}
m_4
m_2=0
for(k in 2:(length(n_k)-1)){
  m_2<-((gamma(k+1)/(gamma(k-1)))*freqrel[k+1])+m_2
}
m_1=0
for(k in 1:(length(n_k)-1)){
  m_1<-((gamma(k+1)/(gamma(k)))*freqrel[k+1])+m_1
}


#ZIP --------------------
(lambZI <- m_2/m_1)
(omegaZI <- m_1/lambZI)
#ZINB --------------------
rZINB <- (m_2^2)/((n_bar*m_3)-m_2^2)-1
pqZINB <- ((n_bar*m_3)-m_2^2)/(n_bar*m_2)
pZINB <- (1+pqZINB)^(-1)
wZINB <- n_bar/(rZINB*pqZINB)

#Mixtura Distribucion Poisson---------------------------
lambda<-n_bar
#Calculo de los lambdas-------------------------
coef1<-((n_bar*m_3)-m_2^2)
coef2<-(m_3-(n_bar*m_2))
coef3<-(m_2-(n_bar^2))
# Construccion del polinomio
polLam <- polynomial(c(coef1, -coef2, coef3))
r1 <- solve(polLam) # raíces
# lambda1 -------------------------------------------
lam1 <- max(r1)
#Lambda 2--------------------------------------------
lam2 <- min(r1)
# Theta ----------------------------------------------
theta<- coef2/coef3
# Gamma----------------------------------------------
gamma<- coef1/coef3
#omega---------------------

omega_1 <- (n_bar-lam2)/(lam1-lam2)
lam1;lam2;omega_1
####### Mixtura de Binomiales Negativas--------------
# parametros ----------------------------------------
q_bin <- 1-(n_bar/s_cua)
r_bin <- n_bar^2/(s_cua-n_bar)
q_bin;r_bin
#Calculo del r-------------------------
a1 <- n_bar*m_2*m_3
coef1<-(m_2^3+m_3^2-m_2*m_4+n_bar^2*m_4-2*a1)
coef2<-(7*m_2^3+4*m_3^2-3*m_2*m_4+4*n_bar^2*m_4 -12*a1)
coef3<-(16*m_2^3+3*m_3^2-2*m_2*m_4+5*n_bar^2*m_4-22*a1)
coef4<-(2*(6*m_2^3+n_bar^2*m_4-6*a1))
# Construccion del polinomio
polr <- polynomial(c(coef4, coef3, coef2, coef1))
r2 <- solve(polr)
r=max(r2)
denpar <- ((r+2)*((r*m_2)-((r+1)*n_bar^2)))
# Theta--------------------
thetabn <- ((r*m_3)-((r+2)*n_bar*m_2))/denpar
# gamma ---------------------
gammabn <-((n_bar*m_3)-((r+1)^(-1)*(r+2)*m_2^2))/denpar
## Estimacion de P
p_1 <- 1/((thetabn/2)+sqrt((thetabn/2)^2-gammabn)+1)
p_2 <- 1/((thetabn/2)-sqrt((thetabn/2)^2-gammabn)+1)
## Estimacion del omega
omega_2 <- ((n_bar*p_2-r*(1-p_2))*p_1)/(r*(p_2-p_1))
p_1;p_2;r;omega_2
#Tabla de comparacion-----------------------------------
#Estimacion de las probabilidades ZIP
Est_poisZI1 = (1-omegaZI)*c(1,0,0,0,0,0,0)
for (i in 1: (k+1)){
  Est_poisZI1[i] = Est_poisZI1[i] + omegaZI*dpois(x = i-1,lambda = lambZI)
}
#Estimacion No. de casos con Poisson
est_pZI<-round(Est_poisZI1*Tcasos,3)
#Casos Totales con Poisson
sum(est_pZI)

#Se agrupan los datos
ajuspoisZI <- c(est_pZI[1:4],sum(est_pZI[5:7]))
# Se agrupa de la misma manera los observados
n_ka<-c(n_k[1:4],sum(n_k[5:7]))
Chi_pZI<-sum((n_ka-ajuspoisZI)^2/ajuspoisZI)
# Se hace la prueba Chi-cuadrado
Chi_cua(obs = n_ka,esp=ajuspoisZI,gl = 2)

#Estimacion de las probabilidades ZINB
Est_ZINB1 = (1-wZINB)*c(1,0,0,0,0,0,0)
for (i in 1:(k+1)){
  Est_ZINB1[i] = Est_ZINB1[i] + wZINB*dnbinom(x = i-1,size = rZINB,prob = (pZINB))
}

#Estimacion No. de casos con Poisson
est_dZINB<-round(Est_ZINB1*Tcasos,3)
#Casos Totales con Poisson
sum(est_dZINB)
#Se agrupan los datos
ajusZINB <- c(est_dZINB[1:5],sum(est_dZINB[6:7]))
# Se agrupa de la misma manera los observados
n_ka<-c(n_k[1:5],sum(n_k[6:7]))
Chi_ZINB<-sum((n_ka-ajusZINB)^2/ajusZINB)
# Se hace la prueba Chi-cuadrado
Chi_cua(obs = n_ka,esp=ajusZINB,gl = 2)

#Estimacion de las probabilidades
Est_mix_pois = c(0, 0, 0, 0, 0, 0, 0)
for(i in 1:(1+k)){
  Est_mix_pois[i] <-omega_1*dpois(x = i-1,lambda = lam1)+
    (1-omega_1)*dpois(x = i-1,lambda = lam2)
}
#Estimacion No. de casos con Poisson
est_mix_p<-round(Est_mix_pois*Tcasos,3)
#Casos Totales con Poisson
sum(est_mix_p)

#Se agrupan los datos
ajuspois1 <- c(est_mix_p[1:5],sum(est_mix_p[6:7]))
# Se agrupa de la misma manera los observados
n_ka<-c(n_k[1:5],sum(n_k[6:7]))
Chi_mix_p<-sum((n_ka-ajuspois1)^2/ajuspois1)
# Se hace la prueba Chi-cuadrado
Chi_cua(obs = n_ka,esp=ajuspois1,gl = 3)

# Estimacion de probabilidades con binomial
# negativa por el metodo de los momentos
Est_mix_binom = c(0, 0, 0, 0, 0, 0, 0)
for (i in 1:(k+1)){
  dbnp1 <- dnbinom(x = i-1,size = r,prob = (p_1));
  dbnp2 <- dnbinom(x = i-1,size = r,prob = (p_2));
  Est_mix_binom[i] <-dbnp1*omega_2+dbnp2*(1-omega_2)
}
est_mix_binn1 <- round(Est_mix_binom*Tcasos,3)
# Casos Totales con la estimacion binomial negativa
# por el metodo de los momentos
sum(est_mix_binn1)

ajusbinn11 <- c(est_mix_binn1[1:5],sum(est_mix_binn1[6:7]))
# Se agrupa de la misma manera los observados
n_ka<-c(n_k[1:5],sum(n_k[6:7]))
Chi_mix_binn1<-sum((n_ka-ajusbinn11)^2/ajusbinn11)

# Se hace la prueba Chi-cuadrado
Chi_cua(obs = n_ka,esp=ajusbinn11)
### Salida Final
test=c("test",".",round(Chi_pZI,4),round(Chi_ZINB,4),
       round(Chi_mix_p,4),round(Chi_mix_binn1,4))
salida=cbind(k,n_k,est_pZI,est_dZINB,est_mix_p,
             est_mix_binn1)
salida=rbind(salida,test)
salida
