library(MASS)
library(LaplacesDemon)#for rinvwishart()
library(tidyquant)
library(plotly)
library(tidyverse)
library(timetk)
#source("RejectionSamplingOBPrior.R")
# get list of names from a imported yahoo portfolio
source("getNamesFinance.R")
# Gathers stock data in nested list
#source("getData.R")

# Stocks from arbitrary "Sweden top 100" portfolio
stockNames <- getNames("Data/quotes.csv")
# All stocks on Swedish market
allstockNames <- getNames("Data/stocks.csv")

# Create data frame containing adjusted price for each stock and time
# Maybe later

# Doing analysis from post "https://www.codingfinance.com/post/2018-04-20-portfolio-stats/" for 
#Select stocks to analyze

tick <- stockNames[-1] %>% head(6) %>% str_sort()

fullNames <- get_stock_name(tick)

#####End date: 2021-03-29
###Weekly 
##52
#2020-03-26
##78
#2019-09-26-2021-03-29
##104
#2019-03-26-2021-03-29
##130
#2018-09-26-2021-03-29
###Monthly
##52
#2016-11-29
##78
#2014-09-29
##104
#2012-07-29
##130
#2010-05-29

#download price data
price_data <- tq_get(tick,
                     from = '2014-09-29',
                     to = '2021-03-29',
                     get = 'stock.prices')


# Calculate daily returns for assets
ret_data <- price_data %>%
  group_by(symbol) %>%
  tq_transmute(select = adjusted,
               mutate_fun = periodReturn,
               period = "monthly",
               col_rename = "ret")

# Wide format
ret_data_wide <- ret_data %>%
  spread(symbol, value = ret) %>%
  tk_xts()

# Approximate missing values 
# Should not affect daily data but in case of weekly/monthly data the NA approx 
# should be done before gathering to get weekly/monthly
ret_data_wide <- na.approx(ret_data_wide)
ret_data_wide <- ret_data_wide[-1,]
ret_data_long <- t(ret_data_wide)

ret_data_wide
ret_data_wide %>% nrow()



k_set<-c(5,10,25,40)
n_set<-c(50,75,100,130)
n_portfolio <- length(k_set)*length(n_set)

R_P<-matrix(0,n_portfolio,1)
V_P<-matrix(0,n_portfolio,1)
R_S<-matrix(0,n_portfolio,1)
V_S<-matrix(0,n_portfolio,1)
R_B<-matrix(0,n_portfolio,1)
V_B<-matrix(0,n_portfolio,1)
R_BI<-matrix(0,n_portfolio,1)
V_BI<-matrix(0,n_portfolio,1)


x_data<-ret_data_long

bar_x<-x_data%*%bi_n/n
bar_xI<-(n*bar_x+r_0*m0)/(n+r_0)
S<-x_data%*%G%*%t(x_data)
iS<-solve(S)
iSI<-solve(S+S0+n*r_0*(m0-bar_x)%*%t(m0-bar_x)/(n+r_0))

Q_S<-iS-iS%*%J_k%*%iS/sum(t(bi_k)%*%iS%*%bi_k)
s_S<-t(bar_x)%*%Q_S%*%bar_x
Q_SI<-iSI-iSI%*%J_k%*%iSI/sum(t(bi_k)%*%iSI%*%bi_k)
s_SI<-t(bar_xI)%*%Q_SI%*%bar_xI

R_S[b]<- (t(bar_x)%*%iS%*%bi_k)/(t(bi_k)%*%iS%*%bi_k)+igam*s_S/d_n
V_S[b]<- d_n/(t(bi_k)%*%iS%*%bi_k)+(igam^2)*s_S/d_n

R_B[b]<- (t(bar_x)%*%iS%*%bi_k)/(t(bi_k)%*%iS%*%bi_k)+igam*s_S/c_kn
V_B[b]<- c_kn/(t(bi_k)%*%iS%*%bi_k)+(igam^2)*s_S/c_kn

R_BI[b]<- (t(bar_xI)%*%iSI%*%bi_k)/(t(bi_k)%*%iSI%*%bi_k)+igam*s_SI/q_kn
V_BI[b]<- q_kn/(t(bi_k)%*%iSI%*%bi_k)+(igam^2)*s_SI/q_kn


######################################
B<-10^2 #Increase to 10^4 for table
r_0<-100
d_0<-100
gam<-50
igam<-1/gam

k_set<-c(5,10,25,40)
n_set<-c(50,75,100,130)

D_R_B<-matrix(0,4,4)
D_R_BI<-matrix(0,4,4)
D_R_S<-matrix(0,4,4)
D_R_H<-matrix(0,4,4)
D_R_O<-matrix(0,4,4)
D_V_B<-matrix(0,4,4)
D_V_BI<-matrix(0,4,4)
D_V_S<-matrix(0,4,4)
D_V_H<-matrix(0,4,4)
D_V_O<-matrix(0,4,4)

for (i in 1:4)
{
k<-k_set[i]
for (j in 1:4)
{n<-n_set[j]
d_n<-1/(n-1)
c_kn<-1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
q_kn<-1/(n+d_0-2*k-1)+(2*n+r_0+d_0-2*k-1)/((n+r_0)*(n+d_0-2*k-1)*(n+d_0-2*k-2))

bi_n<-matrix(1,n,1)
bi_k<-matrix(1,k,1)
G<-diag(n)-matrix(1,n,n)/n
J_k<-matrix(1,k,k)

R_P<-matrix(0,B,1)
V_P<-matrix(0,B,1)
R_S<-matrix(0,B,1)
V_S<-matrix(0,B,1)
R_B<-matrix(0,B,1)
V_B<-matrix(0,B,1)
R_BI<-matrix(0,B,1)
V_BI<-matrix(0,B,1)
R_H<-matrix(0,B,1)
V_H<-matrix(0,B,1)
R_O<-matrix(0,B,1)
V_O<-matrix(0,B,1)

for (b in 1:B)
{
##############mu and Sigma
mu<-runif(k,-0.01,0.01)
u<-runif(k,0.002,0.005)
rho<-0.6
R<-(1-rho)*diag(k)+rho*J_k

Sig<-diag(u)%*%R%*%diag(u)
iSig<-solve(Sig)
tSig<-t(chol(Sig))
Q_P<-iSig-iSig%*%J_k%*%iSig/sum(t(bi_k)%*%iSig%*%bi_k)
s_P<-t(mu)%*%Q_P%*%mu

R_P[b]<- (t(mu)%*%iSig%*%bi_k)/(t(bi_k)%*%iSig%*%bi_k)+igam*s_P
V_P[b]<- 1/(t(bi_k)%*%iSig%*%bi_k)+(igam^2)*s_P 

########################hyperparameters
eps<-0.01
m0<-mu+eps*runif(k,-0.01,0.01)
del<-0.5
S0<-Sig+del*diag(runif(k,0.002,0.005)^2)

###########data
x<-matrix(rep(mu,n),c(k,n))+tSig%*%matrix(rnorm(k*n),c(k,n))
x_data<-x

bar_x<-x_data%*%bi_n/n
bar_xI<-(n*bar_x+r_0*m0)/(n+r_0)
S<-x_data%*%G%*%t(x_data)
iS<-solve(S)
iSI<-solve(S+S0+n*r_0*(m0-bar_x)%*%t(m0-bar_x)/(n+r_0))

Q_S<-iS-iS%*%J_k%*%iS/sum(t(bi_k)%*%iS%*%bi_k)
s_S<-t(bar_x)%*%Q_S%*%bar_x
Q_SI<-iSI-iSI%*%J_k%*%iSI/sum(t(bi_k)%*%iSI%*%bi_k)
s_SI<-t(bar_xI)%*%Q_SI%*%bar_xI

R_S[b]<- (t(bar_x)%*%iS%*%bi_k)/(t(bi_k)%*%iS%*%bi_k)+igam*s_S/d_n
V_S[b]<- d_n/(t(bi_k)%*%iS%*%bi_k)+(igam^2)*s_S/d_n

R_B[b]<- (t(bar_x)%*%iS%*%bi_k)/(t(bi_k)%*%iS%*%bi_k)+igam*s_S/c_kn
V_B[b]<- c_kn/(t(bi_k)%*%iS%*%bi_k)+(igam^2)*s_S/c_kn

R_BI[b]<- (t(bar_xI)%*%iSI%*%bi_k)/(t(bi_k)%*%iSI%*%bi_k)+igam*s_SI/q_kn
V_BI[b]<- q_kn/(t(bi_k)%*%iSI%*%bi_k)+(igam^2)*s_SI/q_kn


###hierarchical prior
#hyperparameters
###hierarchical prior
d_h<-n #also tested 100
kappa_h<-n #also tested 100
eps_h1 <- 0.0001#Same eps values as in "Bayesian estimation of the global minimum variance portfolio"
eps_h2 <- 0.0001
xi<-runif(1,-0.001,0.001) #Tested values -100,100 to -0.0001,0.0001
xi_vec<-rep(xi,k)
eta<- 1/rgamma(1,eps_h1,eps_h2)

r_h<-(kappa_h*xi_vec+n*bar_x)/(kappa_h+n)

SH<-Sig+del*diag(runif(k,0.002,0.005)^2)

psi<-SH/eta+S+(kappa_h*n)/(kappa_h+n)*(xi_vec-bar_x)%*%t(xi_vec-bar_x)#S definition.... Remove (n-1) from S

#Number of samples drawn from posterior
B2<-100 #Increase to 1000 for table
iB2<-1/(B2-1)

x<-c()
for (h in 1:B2) {
Sig_h <- rinvwishart(d_h+n, psi)

mu_h<-matrix(rep(r_h,1),c(k,1))+t(chol(Sig_h/(kappa_h+n)))%*%matrix(rnorm(k),c(k,1))

x<-cbind(x,matrix(rep(mu_h,1),c(k,1))+t(chol(Sig_h))%*%matrix(rnorm(k),c(k,1)))
}
x_data_h<-x

bi_b2<-matrix(1,B2,1)
G_b2<-diag(B2)-matrix(1,B2,B2)/(B2-1)

bar_x_h<-x_data_h%*%bi_b2/B2
S_h<-x_data_h%*%G_b2%*%t(x_data_h)
iS_h<-solve(S_h)

Q_S_h<-iS_h-iS_h%*%J_k%*%iS_h/sum(t(bi_k)%*%iS_h%*%bi_k)
s_S_h<-t(bar_x_h)%*%Q_S_h%*%bar_x_h

R_H[b]<- (t(bar_x_h)%*%iS_h%*%bi_k)/(t(bi_k)%*%iS_h%*%bi_k)+igam*s_S_h/iB2 # divide by d_n??
V_H[b]<- iB2/(t(bi_k)%*%iS_h%*%bi_k)+(igam^2)*s_S_h/iB2

###objective-based prior

#Not complete yet, R_O and V_O contain the same values as the hierarchical prior
s2 <- mean(sum(diag(Sig)))
s_ob<-100 #sigma^2_ob in literature 
v_ob<-100
w_ob <- rep(1/k,k)
SOB<-Sig+del*diag(runif(k,0.002,0.005)^2)
m_ob <- s2/s_ob+n
r_ob <- ((s2/s_ob)*gam*Sig%*%w_ob+n*bar_x)/m_ob

#Sig_ob <- RejectionSample(acceptMax = 1, n=n, k = k, mean_ret = bar_x, Sig = Sig, S = S, w = w_ob, s2 = s2, s_ob = s_ob, SOB = SOB, v_ob = v_ob,gamma = gam)



R_O[b]<- (t(bar_x_h)%*%iS_h%*%bi_k)/(t(bi_k)%*%iS_h%*%bi_k)+igam*s_S_h/iB2# divide by d_n??
V_O[b]<- iB2/(t(bi_k)%*%iS_h%*%bi_k)+(igam^2)*s_S_h/iB2
}

D_R_S[i,j]<-mean(abs(R_S-R_P))
D_V_S[i,j]<-mean(abs(V_S-V_P))
D_R_B[i,j]<-mean(abs(R_B-R_P))
D_V_B[i,j]<-mean(abs(V_B-V_P))
D_R_BI[i,j]<-mean(abs(R_BI-R_P))
D_V_BI[i,j]<-mean(abs(V_BI-V_P))
D_R_H[i,j]<-mean(abs(R_H-R_P))
D_V_H[i,j]<-mean(abs(V_H-V_P))
D_R_O[i,j]<-mean(abs(R_O-R_P))
D_V_O[i,j]<-mean(abs(V_O-V_P))
}
}



R_tab<-NULL
for (i in 1:4)
{
  R_tab<-rbind(R_tab,D_R_S[i,],D_R_B[i,],D_R_BI[i,],D_R_H[i,],D_R_O[i,])  
}

V_tab<-NULL
for (i in 1:4)
{
  V_tab<-rbind(V_tab,D_V_S[i,],D_V_B[i,],D_V_BI[i,],D_V_H[i,],D_V_O[i,])  
}

RV_tab <- as.data.frame(cbind(R_tab,V_tab))


base <- ggplot() + xlim(0.0000005, 1)
base + geom_function(fun = function(x) c_kn/t(mean_ret)%*%Q%*%mean_ret*(x-R_gmv_diffuse)^2 + V_gmv_diffuse, aes(colour = "Diffuse")) +
  geom_function(fun = function(x) d_n/t(mean_ret)%*%Q%*%mean_ret*(x-R_gmv_freq)^2 + V_gmv_freq, aes(colour = "Freq")) +
  geom_function(fun = function(x) q_kn/t(mean_ret_c)%*%Q_c%*%mean_ret_c*(x-R_gmv_conjugate)^2 + V_gmv_conjugate, aes(colour = "Conjugate")) +
  #scale_x_continuous(labels = expectedAnnualReturn, limits = c(0, 0.08)) +
  ylim(0,0.8) +
  scale_color_manual(values = c('red','blue','green')) +
  labs(x="Return", y="Variance") +
  coord_flip()











