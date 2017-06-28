///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tendance Variable///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////
// Fonction phi0_i
////////////////////////////////////////////////////////
function phi0_vect = Phi0(n,h,phi1,phi2,lnf,lnf0)
// Fonction support pour le calcul de phi0(i)
// lnforward 	le logatithme du forward
// sigma  	la volatilite
// lambda 	la vitesse de retour a la moyenne	
// h  	  	le pas de temps des simulations 
///////////////////////////////////////////////////////
phi0_vect = zeros(1,n);

phi0_vect(1) = lnf(1)-phi1*lnf0-(exp(phi2/2)/2);
phi0_vect(2:n) = lnf(2:n)-phi1*lnf(1:n-1)-((exp(phi2/2)/2)*(1+phi1^(2*(2:n)-1))/(1+phi1));


////////////////////////////////////////////////////////
// Fonction Simul_Xt_Var
////////////////////////////////////////////////////////
function Xt_var = Simul_Xt_Var(n,h,lambda,sigma,lnforward,lnf0)
// Simulation d'un processus d'OU avec tendance variable
// X(i+1) - phi1*X(i)-phi0_i = sigma_e*N(0,1)
// n     	nombre de simulations
// lnforward 	le logatithme du forward
// sigma  	la volatilite
// lambda 	la vitesse de retour a la moyenne	
// h  	  	le pas de temps des simulations 
///////////////////////////////////////////////////////
phi0_vect = zeros(1,n);
Xt_var = zeros(1,n);

phi1 = exp(-lambda*h);
sigma_e = sqrt((1-phi1^2)/(2*lambda))*sigma;
phi2 = log(sigma_e^2);

phi0_vect = Phi0(n,h,phi1,phi2,lnforward,lnf0);

Xt_var(1) = phi0_vect(1)+(sigma_e*rand(1,1,"normal"));

for i=1:n-1
Xt_var(i+1) = phi0_vect(i+1)+(phi1*Xt_var(i))+(sigma_e*rand(1,1,"normal"));
end




////////////////////////////////////////////////////////
// Fonction Calage_Xt_Var
////////////////////////////////////////////////////////
function Param_Theta = Calage_Xt_Var(n,h,X,X0,lnf,lnf0)
// Calage d'un processus d'OU avec tendance variable
// X(i) = phi0(i) + phi1*X(i-1) + sigma_e*normal(0,1)
// n         	   nombre de simulations
// Xt_Var          un echantillon de X(ih)
// h		   le pas de temps
// lnforward 	   le logatithme du forward
// sigma  	   la volatilite
// lambda	   la vitesse de retour a la moyenne	
///////////////////////////////////////////////////////
Y = zeros(1,n);
Y(1:n) = {X0,X(1:n-1)};

Param = zeros(1,3);
Param = MCO_Const(n,h,X,X0);
lambda_init = Param(2);
sigma_init = Param(3);

phi1_init = exp(-lambda_init*h);
sigma2_e_init = ((1-phi1_init^2)/(2*lambda_init))*sigma_init^2;
phi2_init = log(abs(sigma2_e_init));

//disp([lambda_init sigma_init],'lambda_init sigma_init');

[Lv, Param, gradopt] = optim(list(Log_vraisemblance_Var,n,h,X,X0,lnf,lnf0),[phi1_init phi2_init],'gc')

Param_Theta(1) = -log(abs(Param(1)))/h;//lambda
Param_Theta(2) = exp(Param(2)/2)*sqrt(abs(2*Param_Theta(1)/(1-exp(-2*Param_Theta(1)*h))));//sigma


////////////////////////////////////////////////////////
// Fonction Log_vraisemblance_Var
////////////////////////////////////////////////////////
function [Lv, grad, ind] = Log_vraisemblance_Var(x,ind,n,h,X,X0,lnf,lnf0)
// Log-vraisemblance pour un processus d'OU avec tendance variable
// X(i) = phi0(i) + phi1*X(i-1) + sigma_e*normal(0,1)
// n         	   nombre de simulations
// X               un echantillon de (Xt)
// lnforward 	   le logatithme du forward
// grad     	   le gradient de la Log-vraisemblance 
// Lv  	   	   la vraisemblance
// x		   parametres a optimiser lambda et sigma
////////////////////////////////////////////////////////
phi1 = x(1);
phi2 = x(2);
cst = %pi;
phi0_vect = zeros(1,n);

Y = zeros(1,n);
Y(1:n) = {X0,X(1:n-1)};
F = zeros(1,n);
F = {lnF0,lnf(1:n-1)};


H = (1+phi1^(2*(1:n)-1))/(1+phi1);
G = ((1+phi1)*(2*(1:n)-1).*phi1^(2*(1:n)-2)-(1+phi1^(2*(1:n)-1)))/((1+phi1)^2);

phi0_vect = Phi0(n,h,phi1,phi2,lnf,lnf0);

Lv = -(n/2)*(log(2*cst)+phi2+mean((X-phi1*Y-phi0_vect)^2*exp(-phi2)));

grad(1) = n*mean((X-phi1*Y-phi0_vect).*(Y-F-(G*exp(phi2)/2))*exp(-phi2));//(-phi1*h)*
grad(2) = (n/2)*(mean((X-phi1*Y-phi0_vect)^2*exp(-phi2)-(X-phi1*Y-phi0_vect).*H)-1);

//disp([lambda sigma phi1 phi2 Lv grad(1) grad(2)],'lambda sigma phi1 phi2 Lv grad12')

Lv = -Lv;
grad = -grad;



////////////////////////////////////////////////////////
// Fonction Residu_Var
////////////////////////////////////////////////////////
function [Resid, Autocorr, BP_Test, Probac, JB_Test, S, K, Proban] = Residu_Var(n,h,theta,X,X0,lnf,lnf0,alpha)
// Calcul des residus sur le modele a tendance variable 
// Test de Box-Pierce pour les residus
// et test d'abscence d'aurocorrelation
// X(i+1) - phi1*X(i)-phi0_i= sigma_e*N(0,1)
// X               un echantillon de (Xt)
// theta	   parametres a optimiser mu,lambda et sigma
// n         	   nombre de simulations
// lnforward 	   le logatithme du forward
// h		   le pas de temps	
// alpha	   le quantil de probabilite de rejet
///////////////////////////////////////////////////////
phi0_vect = zeros(1,n);
Resid = zeros(1,n);
Y = zeros(1,n);
Zi = zeros(1,n);

phi1 = exp(-theta(1)*h);
sigma_e = sqrt((1-exp(-2*h*theta(1)))/(2*theta(1)))*theta(2);
phi2 = log(sigma_e^2)
phi0_vect = Phi0(n,h,phi1,phi2,lnf,lnf0);


Y(1:n) = {X0,X(1:n-1)};
Resid = (X-phi1*Y-phi0_vect)/sigma_e;
k = floor(n/3);
moyenne = mean(Resid);
var = mean((Resid-moyenne)^2);
Autocorr = zeros(1,k);

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
