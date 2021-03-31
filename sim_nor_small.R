library(MASS)
set.seed(999)
B<-10^2
r_0<-100
d_0<-100
gam<-50
igam<-1/gam

k_set<-c(5,10,25,40)
n_set<-c(50,75,100,130)

D_R_B<-matrix(0,4,4)
D_R_BI<-matrix(0,4,4)
D_R_S<-matrix(0,4,4)
D_V_B<-matrix(0,4,4)
D_V_BI<-matrix(0,4,4)
D_V_S<-matrix(0,4,4)

D_R_S_mve<-matrix(0,4,4)
D_V_S_mve<-matrix(0,4,4)
D_R_S_mcd<-matrix(0,4,4)
D_V_S_mcd<-matrix(0,4,4)

for (i in 1:4)
{
k<-k_set[i]
i
for (j in 1:4)
{n<-n_set[j]
j
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

R_S_mve<-matrix(0,B,1)
V_S_mve<-matrix(0,B,1)

R_S_mcd<-matrix(0,B,1)
V_S_mcd<-matrix(0,B,1)

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
tx<-t(x)
S_mve_res<-cov.mve(tx)
bar_x_mve<-S_mve_res$center
S_mve<-S_mve_res$cov
S_mcd_res<-cov.mcd(tx)
bar_x_mcd<-S_mcd_res$center
S_mcd<-S_mcd_res$cov

iS_mve<-solve(S_mve)
Q_S_mve<-iS_mve-iS_mve%*%J_k%*%iS_mve/sum(t(bi_k)%*%iS_mve%*%bi_k)
s_S_mve<-t(bar_x_mve)%*%Q_S_mve%*%bar_x_mve
R_S_mve[b]<- (t(bar_x_mve)%*%iS_mve%*%bi_k)/(t(bi_k)%*%iS_mve%*%bi_k)+igam*s_S_mve
V_S_mve[b]<- 1/(t(bi_k)%*%iS_mve%*%bi_k)+(igam^2)*s_S_mve

iS_mcd<-solve(S_mcd)
Q_S_mcd<-iS_mcd-iS_mcd%*%J_k%*%iS_mcd/sum(t(bi_k)%*%iS_mcd%*%bi_k)
s_S_mcd<-t(bar_x_mcd)%*%Q_S_mcd%*%bar_x_mcd
R_S_mcd[b]<- (t(bar_x_mcd)%*%iS_mcd%*%bi_k)/(t(bi_k)%*%iS_mcd%*%bi_k)+igam*s_S_mcd
V_S_mcd[b]<- 1/(t(bi_k)%*%iS_mcd%*%bi_k)+(igam^2)*s_S_mcd

}

D_R_S_mve[i,j]<-mean(abs(R_S_mve-R_P))
D_V_S_mve[i,j]<-mean(abs(V_S_mve-V_P))
D_R_S_mcd[i,j]<-mean(abs(R_S_mcd-R_P))
D_V_S_mcd[i,j]<-mean(abs(V_S_mcd-V_P))

D_R_S[i,j]<-mean(abs(R_S-R_P))
D_V_S[i,j]<-mean(abs(V_S-V_P))
D_R_B[i,j]<-mean(abs(R_B-R_P))
D_V_B[i,j]<-mean(abs(V_B-V_P))
D_R_BI[i,j]<-mean(abs(R_BI-R_P))
D_V_BI[i,j]<-mean(abs(V_BI-V_P))
}
}

R_tab<-NULL
for (i in 1:4)
{
  R_tab<-rbind(R_tab,D_R_B[i,],D_R_BI[i,],D_R_S[i,],D_R_S_mve[i,],D_R_S_mcd[i,])  
}

V_tab<-NULL
for (i in 1:4)
{
  V_tab<-rbind(V_tab,D_V_B[i,],D_V_BI[i,],D_V_S[i,],D_V_S_mve[i,],D_V_S_mcd[i,])  
}

print(cbind(R_tab,V_tab),digits=4)

print(R_tab,digits=4)

print(V_tab,digits=4)
