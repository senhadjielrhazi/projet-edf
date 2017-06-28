///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Modele hyperbolqique/////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////
// Fonction Densite_NIG
////////////////////////////////////////////////////////
function densite_NIG = Densite_NIG(x,alpha,beta,delta,mu)
// Calcul de la densité NIG
// x         	   	vecteur position
// alpha beta delta mu 	parametres du NIG
///////////////////////////////////////////////////////
n = length(x);
densite_NIG = zeros(1,n);
gama = sqrt(alpha^2-beta^2);
y = sqrt(delta^2+(x-mu)^2);
densite_NIG = ((alpha*delta)/%pi)*exp(delta*gama+beta*(x-mu)).*(besselk(1,alpha*y)')./y;
endfunction


////////////////////////////////////////////////////////
// Fonction Densite_VG
////////////////////////////////////////////////////////
function densite_VG = Densite_VG(x,lambda,alpha,beta,mu)
// Calcul de la densité VG
// x         	   		vecteur position
// lambda alpha beta mu 	parametres du VG
///////////////////////////////////////////////////////
n = length(x);
densite_VG = zeros(1,n);
gama = sqrt(alpha^2-beta^2);
y = abs(x-mu);
densite_VG = ((gama^(2*lambda))/(sqrt(%pi)*((2*alpha)^(lambda-0.5))*gamma(lambda)))*exp(beta*(x-mu)).*(besselk(lambda-0.5,alpha*y)').*y^(lambda-0.5);
endfunction


////////////////////////////////////////////////////////
// Fonction Phi_NIG
////////////////////////////////////////////////////////
function phi_NIG = Phi_NIG(x,alpha,beta,delta,mu)
// Calcul de l exposant caracteristique laplace du NIG
// x         	   		vecteur position
// alpha beta delta mu 		parametres du NIG
///////////////////////////////////////////////////////
n = length(x);
phi_NIG = zeros(1,n);
gama1 = sqrt(alpha^2-beta^2);
gama2 = sqrt(alpha^2-(beta+x)^2);
phi_NIG = mu*x+delta*(gama1-gama2);
endfunction


////////////////////////////////////////////////////////
// Fonction Phi_VG
////////////////////////////////////////////////////////
function phi_VG = Phi_VG(x,lambda,alpha,beta,mu)
// Calcul de l exposant caracteristique laplace du VG
// x         	   		vecteur position
// lambda alpha beta mu 	parametres du VG
///////////////////////////////////////////////////////
n = length(x);
phi_VG = zeros(1,n);
gama1 = sqrt(alpha^2-beta^2);
gama2 = sqrt(alpha^2-(beta+x)^2);
phi_VG = mu*x+2*lambda*log(gama1./gama2);
endfunction


////////////////////////////////////////////////////////
// Fonction Var_NIG
////////////////////////////////////////////////////////
function var_NIG = Var_NIG(n,alpha,beta,delta,mu)
// Simulation dune variable aleatoire NIG
// n         	   		nombre de simulations
// alpha beta delta mu 		parametres du NIG
///////////////////////////////////////////////////////
gama = sqrt(alpha^2-beta^2);
V = rand(1,n,'normal')^2;
Z1 = (delta/gama)+(1/(2*gama^2))*(V-sqrt(V^2+4*gama*delta*V));
Z2 = (delta/gama)+(1/(2*gama^2))*(V+sqrt(V^2+4*gama*delta*V));
p1 = delta*ones(1,n)./(delta+gama*Z1);
U = rand(1,n,'uniform');
Z = Z1.*(U<p1)+Z2.*(1-(U<p1));
var_NIG = mu + beta*Z + sqrt(Z).*rand(1,n,'normal');
endfunction


////////////////////////////////////////////////////////
// Fonction Var_VG
////////////////////////////////////////////////////////
function var_VG = Var_VG(n,lambda,alpha,beta,mu)
// Simulation dune variable aleatoire VG
// n         	   		nombre de simulations
// lambda alpha beta mu 	parametres du VG
// lambda 			>= 0
///////////////////////////////////////////////////////
gama = sqrt(alpha^2-beta^2);
m = 0;
Z = zeros(1,n);
k1 = 0;
k2 = 0;
l  = 0;
lambda1 = 0;
c = exp(1)*sqrt(6/%pi);
if(lambda < 1)
lambda1 = lambda + 1;
else
lambda1 = lambda;
end
while(m < n)
U1 = rand(1,2*n,'uniform');
index = find((U1<>0)&(U1<>1));
U = U1(index);
V1 = (U-1/2)./sqrt(U.*(1-U))*sqrt(3*lambda1-3/4)+(lambda1-1);
index = find(V1>=0);
V = V1(index);
k1 = length(index);
if(k1<>0)
index = find(V^(lambda1-1).*exp(-V)/gamma(lambda1) > c*rand(1,k1,'uniform').*((1+(V-lambda1+1)^2/(3*lambda1-3/4))^(-3/2)/(2*sqrt(3*lambda1-3/4))));
k2 = length(index);
if(k2>(n-l))
Z((l+1):n) = V(index(1:(n-l)));
m = n;
else
Z((l+1):(l+k2)) = V(index);
m = l+k2;
l = l+k2;
end
end
end
if(lambda < 1)
Z = Z.*rand(1,n,'uniform')^(1/lambda);
end
var_VG = mu + beta*2*Z/gama^2 + sqrt(2*Z/gama^2).*rand(1,n,'normal');
endfunction


////////////////////////////////////////////////////////
// Fonction Simul_XT_NIG
////////////////////////////////////////////////////////
function simul_XT_NIG = Simul_XT_NIG(m,n,T0,Tf,t,alpha,beta,delta,mu,a,sigma)
// Simuler la variable XT pour option europeenne NIG
// n         	   		nombre de pas
// m         	   		nombre de tirages
// alpha beta delta mu		parametres du NIG 
// kk T0 Tf T			parametre du Call
// a sigma			parametre du processus
///////////////////////////////////////////////////////
h = (T0-t)/n;
Ti = t + (0:(n-1))*h;
simul_XT_NIG = zeros(1,m);
Val = -Phi_NIG(sigma*exp(-a*(Tf-Ti)),alpha,beta,delta,mu)*h;
for i=1:m
Y = sigma*exp(-a*(Tf-Ti)).*Var_NIG(n,alpha,beta,delta*h,mu*h);
simul_XT_NIG(i) = sum(Y+Val);
end
endfunction


////////////////////////////////////////////////////////
// Fonction Simul_XT_VG
////////////////////////////////////////////////////////
function simul_XT_VG = Simul_XT_VG(m,n,T0,Tf,t,lambda,alpha,beta,mu,a,sigma)
// Simuler la variable XT pour option europeenne VG
// n         	   		nombre de pas
// m         	   		nombre de tirages
// lambda alpha beta mu		parametres du VG 
// kk T0 Tf T			parametre du Call
// a sigma			parametre du processus
///////////////////////////////////////////////////////
h = (T0-t)/n;
Ti = t + (0:(n-1))*h;
simul_XT_VG = zeros(1,m);
Val = -Phi_VG(sigma*exp(-a*(Tf-Ti)),lambda,alpha,beta,mu)*h;
for i=1:m
Y = sigma*exp(-a*(Tf-Ti)).*Var_VG(n,lambda*h,alpha,beta,mu*h);
simul_XT_VG(i) = sum(Y+Val);
end
endfunction


////////////////////////////////////////////////////////
// Fonction Kernel_NIG
////////////////////////////////////////////////////////
function kernel_NIG = Kernel_NIG(s,k,eta,T0,Tf,t,alpha,beta,delta,mu,a,sigma)
// Calcul la fonction d integration du FTT
// h				pas de temps
// N         	   		nombre de pas
// alpha beta delta mu		parametres du NIG 
// k T0 Tf T			parametre du Call
// a sigma			parametre du processus
// A				borne d integration
// eta				constante
///////////////////////////////////////////////////////
I = %i;
reel  = integrate('real(Phi_NIG((I*s+1+eta)*sigma*exp(-a*(Tf-y)),alpha,beta,delta,mu)-(I*s+1+eta)*Phi_NIG(sigma*exp(-a*(Tf-y)),alpha,beta,delta,mu))','y',t,T0);
imagi  = integrate('imag(Phi_NIG((I*s+1+eta)*sigma*exp(-a*(Tf-y)),alpha,beta,delta,mu)-(I*s+1+eta)*Phi_NIG(sigma*exp(-a*(Tf-y)),alpha,beta,delta,mu))','y',t,T0);
kernel_NIG = (exp(reel)/((eta^2+eta-s^2)^2+s^2*(2*eta+1)^2))*((eta^2+eta-s^2)*cos(-s*k+imagi)+(s*(2*eta+1))*sin(-s*k+imagi));
endfunction


////////////////////////////////////////////////////////
// Fonction Kernel_VG
////////////////////////////////////////////////////////
function kernel_VG = Kernel_VG(s,k,eta,T0,Tf,t,lambda,alpha,beta,mu,a,sigma)
// Calcul la fonction d integration du FTT
// h				pas de temps
// N         	   		nombre de pas
// lambda alpha beta mu		parametres du VG 
// k T0 Tf T			parametre du Call
// a sigma			parametre du processus
// A				borne d integration
// eta				constante
///////////////////////////////////////////////////////
I = %i;
reel  = integrate('real(Phi_VG((I*s+1+eta)*sigma*exp(-a*(Tf-y)),lambda,alpha,beta,mu)-(I*s+1+eta)*Phi_VG(sigma*exp(-a*(Tf-y)),lambda,alpha,beta,mu))','y',t,T0);
imagi  = integrate('imag(Phi_VG((I*s+1+eta)*sigma*exp(-a*(Tf-y)),lambda,alpha,beta,mu)-(I*s+1+eta)*Phi_VG(sigma*exp(-a*(Tf-y)),lambda,alpha,beta,mu))','y',t,T0);
kernel_VG = (exp(reel)/((eta^2+eta-s^2)^2+s^2*(2*eta+1)^2))*((eta^2+eta-s^2)*cos(-s*k+imagi)+(s*(2*eta+1))*sin(-s*k+imagi));
endfunction


////////////////////////////////////////////////////////
// Fonction Call_NIG_Fermee
////////////////////////////////////////////////////////
function call_NIG_Fermee = Call_NIG_Fermee(r,K,x,eta,T0,Tf,t,alpha,beta,delta,mu,a,sigma, epsilon,h)
// Calcul du prix call NIG avec une formule fermee 
// h				pas de temps
// N         	   		nombre de pas
// alpha beta delta mu		parametres du NIG 
// x K T0 Tf T			parametre du Call
// a sigma			parametre du processus
// A				borne d integration
// epsilon			erreur sur le prix
// eta				constante
///////////////////////////////////////////////////////
f1 = 0;
A1 = 0;
M = 0;
Integrale = 0;
k =  log(K/x);
A = x*exp(-eta*k-r*(T0-t)+integrate('Phi_NIG((eta+1)*sigma*exp(-a*(Tf-q)),alpha,beta,delta,mu)-(eta+1)*Phi_NIG(sigma*exp(-a*(Tf-q)),alpha,beta,delta,mu)','q',t,T0))/(%pi*(epsilon/3));
for j=1:1000
A1 = j*A/1000;
f1 = Kernel_NIG(A1,k,eta,T0,Tf,t,alpha,beta,delta,mu,a,sigma);
M = abs(f1)*(A-A1);
if(M<(epsilon/3))
break
end
end
A = A1;
for j=1:1000
A1 = j*A/1000;
f1 = Kernel_NIG(A1,k,eta,T0,Tf,t,alpha,beta,delta,mu,a,sigma);
M = abs(f1)*(A-A1);
if(M<(epsilon/3))
break
end
end
N = floor(A1/h)+1;
Val = zeros(2,N);
for i=1:N
Val(1,i) = (i-1)*h;
Val(2,i) = Kernel_NIG(Val(1,i),k,eta,T0,Tf,t,alpha,beta,delta,mu,a,sigma)
end
Integrale = inttrap(Val(1,:),Val(2,:));
call_NIG_Fermee = (x*exp(-eta*k)*exp(-r*(T0-t))/(%pi))*Integrale;
endfunction


////////////////////////////////////////////////////////
// Fonction Call_VG_Fermee
////////////////////////////////////////////////////////
function call_VG_Fermee = Call_VG_Fermee(r,K,x,eta,T0,Tf,t,lambda,alpha,beta,mu,a,sigma,epsilon,h)
// Calcul du prix call VG avec une formule fermee 
// h				pas de temps
// N         	   		nombre de pas
// lambda alpha beta mu		parametres du VG 
// x K T0 Tf T			parametre du Call
// a sigma			parametre du processus
// A				borne d integration
// epsilon			erreur sur le prix
// eta				constante
///////////////////////////////////////////////////////
f1 = 0;
A1 = 0;
M = 0;
Integrale = 0;
k =  log(K/x);
A = x*exp(-eta*k-r*(T0-t)+integrate('Phi_VG((eta+1)*sigma*exp(-a*(Tf-q)),lambda,alpha,beta,mu)-(eta+1)*Phi_VG(sigma*exp(-a*(Tf-q)),lambda,alpha,beta,mu)','q',t,T0))/(%pi*(epsilon/3));
for j=1:1000
A1 = j*A/1000;
f1 = Kernel_VG(A1,k,eta,T0,Tf,t,lambda,alpha,beta,mu,a,sigma);
M = abs(f1)*(A-A1);
if(M<(epsilon/3))
break
end
end
A = A1;
for j=1:1000
A1 = j*A/1000;
f1 = Kernel_VG(A1,k,eta,T0,Tf,t,lambda,alpha,beta,mu,a,sigma);
M = abs(f1)*(A-A1);
if(M<(epsilon/3))
break
end
end
N = floor(A1/h)+1;
Val = zeros(2,N);
for i=1:N
Val(1,i) = (i-1)*h;
Val(2,i) = Kernel_VG(Val(1,i),k,eta,T0,Tf,t,lambda,alpha,beta,mu,a,sigma)
end
Integrale = inttrap(Val(1,:),Val(2,:));
call_VG_Ferme = (x*exp(-eta*k)*exp(-r*(T0-t))/(%pi))*Integrale;
endfunction









//////////////////////////////////////////////Calage////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////
// Fonction Phi0_NIG
////////////////////////////////////////////////////////
function phi0_NIG = Phi0_NIG(n,h,phi1,omega,alpha,beta,delta,mu)
// Fonction support pour le calcul de phi0(i)
// n         	   	nombre de simulations
// h  	  		le pas de temps des simulations 
// alpha beta delta mu		parametres du NIG
// omega			prime de risque
// a sigma			parametre du processus
///////////////////////////////////////////////////////
//A = zeros(1,n);
//disp([phi1 omega alpha beta delta mu], "phi1,omega,alpha,beta,delta,mu");
//for i = 1:30
//A(i) = integrate('Phi_NIG(omega+phi1^(x/h),alpha,beta,delta,mu)','x',(i-1)*h,i*h)/h-Phi_NIG(omega,alpha,beta,delta,mu);
//end
//phi0_NIG = -cumsum(A)+phi1*cumsum({0,A(1:(n-1))});

phi0_NIG = zeros(1,n);

for i=1:n
phi0_NIG(i) = -integrate('Phi_NIG(omega+phi1^(x/h),alpha,beta,delta,mu)-Phi_NIG(omega,alpha,beta,delta,mu)','x',0,i*h)/h+phi1*integrate('Phi_NIG(omega+phi1^(x/h),alpha,beta,delta,mu)-Phi_NIG(omega,alpha,beta,delta,mu)','x',0,(i-1)*h)/h;
end

endfunction


////////////////////////////////////////////////////////
// Fonction Simul_NIG_Process
////////////////////////////////////////////////////////
function  simul_NIG_Process =  Simul_NIG_Process(n,h,a,sigma,omega,alpha,beta,delta,mu)
// Simuler un processus OU NIG 
// h				pas de temps
// n         	   		nombre de tirages
// alpha beta delta mu		parametres du NIG
// omega			prime de risque
// a sigma			parametre du processus
///////////////////////////////////////////////////////
phi1 = exp(-a*h);
omega_bar = omega/sigma;
alpha_bar = alpha/sigma;
beta_bar = beta/sigma;
delta_bar = delta*h*sigma;
mu_bar = mu*h*sigma;
simul_NIG_Process = zeros(1,n);
phi0_NIG = Phi0_NIG(n,h,phi1,omega_bar,alpha_bar,beta_bar,delta_bar,mu_bar);
simul_NIG_Process(1) = phi0_NIG(1) + Var_NIG(1,alpha_bar,beta_bar,delta_bar,mu_bar);
for i=2:n

simul_NIG_Process(i) = phi1*simul_NIG_Process(i-1) + phi0_NIG(i) + Var_NIG(1,alpha_bar,beta_bar,delta_bar,mu_bar);
end
endfunction



////////////////////////////////////////////////////////
// Fonction Calage_NIG
////////////////////////////////////////////////////////
function [alpha_e, beta_e, delta_e, mu_e] = Calaga_Const_NIG(n,h,Et)
// Calage des residus  du modele OU un facteur NIG 
// h				pas de temps
// n         	   		nombre de pas
// Et				Residus
// alpha beta delta mu		parametres du NIG
///////////////////////////////////////////////////////
alpha_bar_init = 7;
beta_bar_init = 2;
delta_bar_init = 7;
mu_bar_init = -3;

//LB = [0.00001 -40 0.00001 -100 ];
//UB = [60 40 70 100];'b',LB,UB,


//disp([alpha_bar_init beta_bar_init delta_bar_init mu_bar_init],"phi1_init omega_bar_init alpha_bar_init beta_bar_init delta_bar_init mu_bar_init");
[Lv, Param, gradopt] = optim(list(Log_vraisemblance_NIG,Et,n,h),[alpha_bar_init beta_bar_init delta_bar_init mu_bar_init],'gc')

alpha_e = Param(1);
beta_e = Param(2);
delta_e = Param(3);
mu_e = Param(4);
endfunction



////////////////////////////////////////////////////////
// Fonction Log_vraisemblance_NIG
////////////////////////////////////////////////////////
function [Lv, grad, ind] = Log_vraisemblance_NIG(x,ind,X,n,h)
// Log-vraisemblance pour un processus d'OU NIG
// h				pas de temps
// X				Residus
// n				nombre d echantillon
// grad     	   		le gradient de la Log-vraisemblance 
// Lv  	   	   		la vraisemblance
// x		   		parametres a optimiser phi1, omega...
////////////////////////////////////////////////////////
alpha = x(1);
beta = x(2);
delta = x(3);
mu = x(4);

//disp([alpha beta delta mu],"alpha beta delta mu");


if(abs(beta)>=abs(alpha) | alpha<0 | delta<0)// | alpha_bar+delta_bar>50)
x(1) = 10;
x(2) = 0;
x(3) = 10;
Lv = -100;
grad = 100*ones(1,4);
disp([alpha beta delta],'probleme');
else
gama = sqrt(alpha^2-beta^2);

S = (X-mu*ones(1,n))/(delta*h);
P = sqrt(1+S^2);
K = besselk(1,alpha*h*delta*P)'
R = besselk(2,alpha*h*delta*P)'./K;
//disp(size(K./P),"R");
Lv = -n*log(%pi)+n*log(alpha)+n*delta*h*gama+n*mean(beta*delta*h*S+log(K./P));

grad(1) = n*(2/alpha + delta*h*alpha/gama) - n*mean(delta*h*P.*R);
grad(2) = -n*delta*h*beta/gama + n*mean(delta*h*S);
grad(3) = n*(1/delta + gama*h) - n*mean(alpha*h*R./P);
grad(4) = -n*beta*h + n*mean(alpha*h*S.*R./P);

//disp("ok");

disp([Lv grad(1) grad(2) grad(3) grad(4) alpha beta delta mu],"lv grad alpha beta delta mu");

//Lv = Log_vraisemblance_NIG_Value({phi1,omega_bar,alpha_bar,beta_bar,delta_bar,mu_bar},n,h,X,X0);
//grad = numdiff(list(Log_vraisemblance_NIG_Value,n,h,X,X0),{phi1,omega_bar,alpha_bar,beta_bar,delta_bar,mu_bar});
end


Lv = -Lv;
grad = -grad;

endfunction



////////////////////////////////////////////////////////
// Fonction Phi0_NIG_Func
////////////////////////////////////////////////////////
function phi0_NIG_Func = Phi0_NIG_Func(theta,n,h,i)
// Fonction support pour le calcul de phi0(i)
// n         	   	nombre de simulations
// h  	  		le pas de temps des simulations 
// alpha beta delta mu		parametres du NIG
// omega			prime de risque
// a sigma			parametre du processus
///////////////////////////////////////////////////////
phi1  = theta(1);
omega = theta(2);
alpha = theta(3);
beta  = theta(4);
delta = theta(5);
mu    = theta(6);

phi0_NIG_Func = -integrate('Phi_NIG(omega+phi1^(x/h),alpha,beta,delta,mu)-Phi_NIG(omega,alpha,beta,delta,mu)','x',0,i*h)/h+phi1*integrate('Phi_NIG(omega+phi1^(x/h),alpha,beta,delta,mu)-Phi_NIG(omega,alpha,beta,delta,mu)','x',0,(i-1)*h)/h;

endfunction









////////////////////////////////////////////////////////
// Fonction Phi0_NIG
////////////////////////////////////////////////////////
function phi0_VG = Phi0_VG(n,h,phi1,omega,lambda,alpha,beta,mu)
// Fonction support pour le calcul de phi0(i)
// n         	   		nombre de simulations
// h  	  			le pas de temps des simulations 
// lambda alpha beta mu		parametres du VG
// omega			prime de risque
// a sigma			parametre du processus
///////////////////////////////////////////////////////
phi0_VG = zeros(1,n);
for i=1:n
phi0_VG(i) = -integrate('Phi_VG(omega+phi1^(x/h),lambda,alpha,beta,mu)-Phi_VG(omega,lambda,alpha,beta,mu)','x',0,i*h)/h+phi1*integrate('Phi_VG(omega+phi1^(x/h),lambda,alpha,beta,mu)-Phi_VG(omega,lambda,alpha,beta,mu)','x',0,(i-1)*h)/h;
end
endfunction


////////////////////////////////////////////////////////
// Fonction Simul_NIG_Process
////////////////////////////////////////////////////////
function  simul_VG_Process =  Simul_VG_Process(n,h,a,sigma,omega,lambda,alpha,beta,mu)
// Simuler un processus OU VG 
// h				pas de temps
// n         	   		nombre de tirages
// lambda alpha beta mu		parametres du VG
// omega			prime de risque
// a sigma			parametre du processus
///////////////////////////////////////////////////////
phi1 = exp(-a*h);
omega_bar = omega/sigma;
alpha_bar = alpha/sigma;
beta_bar = beta/sigma;
lambda_bar = lambda*h;
mu_bar = mu*h*sigma;
simul_VG_Process = zeros(1,n);
phi0_VG = Phi0_VG(n,h,phi1,omega_bar,lambda_bar,alpha_bar,beta_bar,mu_bar);
cc = Var_VG(1,lambda_bar,alpha_bar,beta_bar,mu_bar);
simul_VG_Process(1) = phi0_VG(1) + Var_VG(1,lambda_bar,alpha_bar,beta_bar,mu_bar);
for i=2:n
simul_VG_Process(i) = phi1*simul_VG_Process(i-1) + phi0_VG(i) + Var_VG(1,lambda_bar,alpha_bar,beta_bar,mu_bar);
end
endfunction


////////////////////////////////////////////////////////
// Fonction Phi0_VG_Func
////////////////////////////////////////////////////////
function phi0_VG_Func = Phi0_VG_Func(theta,n,h,i)
// Fonction support pour le calcul de phi0(i)
// n         	   		nombre de simulations
// h  	  			le pas de temps des simulations 
// lambda alpha beta mu		parametres du VG
// omega			prime de risque
// a sigma			parametre du processus
///////////////////////////////////////////////////////
phi1   = theta(1);
omega  = theta(2);
lambda = theta(3);
alpha  = theta(4);
beta   = theta(5);
mu     = theta(6);

phi0_VG_Func = -integrate('Phi_VG(omega+phi1^(x/h),lambda,alpha,beta,mu)-Phi_VG(omega,lambda,alpha,beta,mu)','x',0,i*h)/h+phi1*integrate('Phi_VG(omega+phi1^(x/h),lambda,alpha,beta,mu)-Phi_VG(omega,lambda,alpha,beta,mu)','x',0,(i-1)*h)/h;

endfunction


