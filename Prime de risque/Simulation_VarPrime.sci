clear;

chdir('J:\Stage\Simulations')
getf('Modele un facteur gaussien\Prime de risque\VarPrime_Function.sci');disp('getf done');
getf('Modele un facteur gaussien\Cas Constant\Const_Function.sci');disp('getf done');



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tendance Variable////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
// donnees
////////////////////////////////////////////////////////
// f(0,t) la courbe forward en 0 ou plutot log(f(0,t))
// sigma = [100%,400%] annuelle
// h =1  heure, 1 jour
// lambda = [1,25] jours ou [10,200] annee
///////////////////////////////////////////////////////

//simulation heure par heure
lambda = 100;
h = 1/(24*365);
sigma = 4;
prime = 20;
lnF0 = 0;


// 1: simulation et prime de risque constante

disp([lambda sigma prime],'vrais valeurs de lambda, sigma et prime pure');

ForwardFile = "Data\EDF_Forward_01012002_31122002.txt" ;

lnForward = log((read(ForwardFile,-1,1))');
n = length(lnForward)-1;
lnForwardEdf = lnForward(2:(n+1));
Xt_VarEdf = zeros(1,n);
X0 = 0;
lnF0 = lnForward(1);

disp(h,'le pas de temps h');
disp(n,'le nombre d observation');

lnForwardPrime = lnForwardEdf + (sigma^2*prime/lambda);
lnF0Prime = lnF0 + (sigma^2*prime/lambda);

Xt_VarEdf = Simul_Xt_Var(n,h,lambda,sigma,lnForwardPrime,lnF0Prime);

//2: Calage par EMV
Theta_EMV = Calage_Xt_VarPrime(n,h,Xt_VarEdf,X0,lnForwardEdf,lnF0);
disp(Theta_EMV,'Les parametres : lambda, sigma prime pure estime par EMV estime');


xbasc(0);
xset("window",0);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(1:n,Xt_VarEdf,1:3,"061","",[0,0,20,20]);
xtitle('Evolution du logarithme du prix Xt');

xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(1:n,exp(Xt_VarEdf),1:3,"061","",[0,0,20,20]);
xtitle('Evolution du prix spot St');

xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(1:n,lnForwardEdf,1:3,"061","",[0,0,20,20]);
xtitle('Evolution du prix forward F(t,T)');

//3: autocorrelation des residus et test de box-pierce
alpha = 0.05;
[Resid, Autocorr, BP_Test, Probac, JB_Test, S, K, Proban] = Residu_VarPrime(n,h,Theta_EMV,Xt_VarEdf,X0,lnForwardEdf,lnF0,alpha);

xbasc(1);
xset("window",1);
xsetech([0,1/4,0.5,1/2],[-1,1,-1,1]);
histplot(100,Resid,1:3,"061","",[0,0,20,20]);
xtitle('Histogramme des residus');

xsetech([0.5,1/4,0.5,1/2],[-1,1,-1,1]);
plot2d(1:length(Autocorr),Autocorr,1:3,"061","",[0,0,20,20]);
xtitle('La fonction d autocorrelation en fonction du retard');

disp(Probac,'La probabilite d absence de autocorrelation');
disp(BP_Test,'Test Box-Pierce 1:absence d autocorrelation, 0:autocorrelation');

disp(Proban,'La probabilite de normalite');
disp(JB_Test,'Test Jarque Bera 1:normalite des residus, 0:rejet');

disp([S K],'Le skewness et le kurtosis');


