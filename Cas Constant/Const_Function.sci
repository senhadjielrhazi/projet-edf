///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tendance Constante///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////
// Fonction Simul_Xt_Const
////////////////////////////////////////////////////////
function Xt_const = Simul_Xt_Const(n,h,mu,lambda,sigma)
// Simulation d'un processus d'OU avec tendance constante
// X(i+1) - phi1*X(i)-phi0 = sigma_e*N(0,1)
// n      nombre de simulations
// mu     coeff de retour a la moyenne 
// sigma  la volatilite
// lambda la vitesse de retour a la moyenne	
// h  	  le pas de temps des simulations 
///////////////////////////////////////////////////////
Xt_const = zeros(1,n);
phi1 = exp(-lambda*h);
phi0 = mu*(1-exp(-lambda*h));
sigma_e = sqrt((1-exp(-2*h*lambda))/(2*lambda))*sigma;

Xt_const(1) = phi0+sigma_e*rand(1,1,"normal");

for i=1:n-1
Xt_const(i+1) = phi0+(phi1*Xt_const(i))+(sigma_e*rand(1,1,"normal"));
end



////////////////////////////////////////////////////////
// Fonction Calage_Xt_Const
////////////////////////////////////////////////////////
function Param_Theta = Calage_Xt_Const(n,h,X,X0)
// Calage d'un processus d'OU avec tendance constante
// n         	   nombre de simulations
// X               un echantillon de X(ih)
// h		   le pas de temps	
///////////////////////////////////////////////////////
Y = zeros(1,n);
Y(1:n) = {X0,X(1:n-1)};
mu_init = mean(X);
phi1_init = (((X(2)-mu_init)/(X(1)-mu_init))-((X(1)-mu_init)/mu_init))/2;
phi0_init = mu_init*(1-phi1_init);
sigma2_e_init = mean((X-phi1_init*Y-phi0_init)^2);
phi2_init = log(abs(sigma2_e_init));

[Lv, Param, gradopt] = optim(list(Log_vraisemblance_Const,X,X0),[phi0_init phi1_init phi2_init],'gc')

Param_Theta(1) = Param(1)/(1-Param(2));
Param_Theta(2) = -log(Param(2))/h;
Param_Theta(3) = sqrt(exp(Param(3))*2*(-log(Param(2))/h)/(1-exp(2*log(Param(2)))));



////////////////////////////////////////////////////////
// Fonction Log_vraisemblance_Const
////////////////////////////////////////////////////////
function [Lv, grad, ind] = Log_vraisemblance_Const(x,ind,X,X0)
// Log-vraisemblance pour un processus d'OU avec tendance constante
// X(i+1) - phi1*X(i)-phi0= sigma_e*N(0,1)
// n         	   nombre de simulations
// X               un echantillon de (Xt)
// grad     	   le gradient de la Log-vraisemblance 
// Lv  	   	   la vraisemblance
// x		   parametres a optimiser phi0, phi1 et sigma2_e
////////////////////////////////////////////////////////
phi0 = x(1);
phi1 = x(2);
phi2 = x(3);
cst = %pi;
n = length(X);
Y = zeros(1,n);
Y(1:n) = {X0,X(1:n-1)}; 


Lv = -(n/2)*(log(2*cst)+phi2+mean((X-(phi1*Y)-phi0)^2*exp(-phi2)));

grad(1) = (n)*(mean(X-(phi1*Y)-phi0)*exp(-phi2));
grad(2) = (n)*(mean(Y.*(X-(phi1*Y)-phi0))*exp(-phi2));
grad(3) = (n/2)*(mean((X-(phi1*Y)-phi0)^2*exp(-phi2))-1);

Lv = -Lv;
grad = -grad;



////////////////////////////////////////////////////////
// Fonction MCO_Const
////////////////////////////////////////////////////////
function Param_MCO = MCO_Const(n,h,X,X0)
// Moindres carres sur le modele a tendance constante
// X(i+1) - phi1*X(i)-phi0 = sigma_e*N(0,1)
// X		un echantillon de (Xt)
// Param	parametres a optimiser mu,lambda et sigma
///////////////////////////////////////////////////////
A = ones(n,2);
Y = zeros(n,1);
Y(1:n) = {X0,X(1:n-1)}';
A(:,2) = Y;
B = X';
C = (inv(A'*A))*A'*B;
D = sqrt(mean((B-C(2)*Y-C(1))^2));

Param_MCO(1) = C(1)/(1-C(2));
Param_MCO(2) = -log(C(2))/h;
Param_MCO(3) = D/(sqrt(-h*(1-C(2)^2)/(2*log(C(2)))));



////////////////////////////////////////////////////////
// Fonction Residu_Const
////////////////////////////////////////////////////////
function [Resid, Autocorr, BP_Test, Probac, JB_Test, S, K, Proban] = Residu_Const(n,h,theta,X,X0,alpha)
// Calcul des residus sur le modele a tendance constante 
// Test de Box-Pierce pour les residus
// et test d'abscence d'aurocorrelation
// X(i+1) - phi1*X(i)-phi0= sigma_e*N(0,1)
// X               un echantillon de (Xt)
// theta	   parametres a optimiser mu,lambda et sigma
// n         	   nombre de simulations
// h		   le pas de temps	
// alpha	   le quantil de probabilite de rejet
///////////////////////////////////////////////////////
phi1 = exp(-theta(2)*h);
phi0 = theta(1)*(1-exp(-theta(2)*h));
sigma_e = sqrt((1-exp(-2*h*theta(2)))/(2*theta(2)))*theta(3);
Resid = zeros(1,n);
Y = zeros(1,n);
Y(1:n) = {X0,X(1:n-1)};
Resid = (X-phi1*Y-phi0)/sigma_e;
k = floor(n/3);
moyenne = mean(Resid);
var = mean((Resid-moyenne)^2);
Autocorr = zeros(1,k);
Zi = zeros(1,n);
for i=1:k
Zi = 0*Zi;
Zi(1:n-i) = Resid((i+1):n)-moyenne;
Autocorr(1,i) = mean((Resid-moyenne).*(Zi))/var;
end

//Khi-deux a k-p-q degres de liberte.
BP_Stat = n*k*mean(Autocorr^2);

[p,Probac] = cdfchi("PQ",BP_Stat,k-1);
BP_Test = 1*(Probac >= alpha);//1:accepter hypothese absence autocorrelation 0:rejet

u1 = mean(Resid);
u2 = mean((Resid-u1)^2);
u3 = mean((Resid-u1)^3);
u4 = mean((Resid-u1)^4);

S = u3^2/u2^3;
K = u4/u2^2;


JB_Stat = (n/6)*(S+(1/4)*(K-3)^2);
[p,Proban] = cdfchi("PQ",JB_Stat,2);

JB_Test = 1*(Proban >= alpha);//1:accepter hypothese de normalite 0:rejet
