clear;

chdir('J:\Stage\Simulations\Modele un facteur gaussien')

getf('Cas Constant\Const_Function.sci');disp('getf done');

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tendance Constante///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
// donnees
////////////////////////////////////////////////////////
// mu = 25
// sigma = [100%,400%] annuelle
// h =1  heure, 1 jour
// lambda = [1,25] jours ou [10,200] annee
///////////////////////////////////////////////////////

//simulation heure par heure
mu = 3.35;
lambda = 100;
h = 1/(24*365);
sigma = 4;
n = 365*24;

//simulation jour par jour
//mu = 25;
//lambda = 100;
//h = 1/365;
//sigma = 4;
//n = 365;


//Type de donnees a trouver
Methode = 1; // 0: simulation du log spot et du prix avec saisonnalite, 
	     						// 1: autocorrelation des residus, calage et test de box-pierce, 
	     							     							

disp(h,'le pas de temps h');
disp(n,'le nombre d observation');

disp([mu lambda sigma],'mu_vrai lambda_vrai sigma_vrai');


// 0: simulation du log spot et du prix avec saisonnalite
if Methode == 0
Xt = zeros(1,n); 
Dt = ones(1,n);//saisonnalite

xbasc(0);
xset("window",0);

xsetech([0,1/4,1/2,1/2],[-1,1,-1,1]);
Xt_Const = Simul_Xt_Const(n,h,mu,lambda,sigma);

plot2d(1:n,Xt_Const,1:3,"061","",[0,0,20,20]);
plot2d(1:n,mu*ones(1,n),[2,2,4]);
xtitle('Evolution du log-prix spot');

xsetech([1/2,1/4,1/2,1/2],[-1,1,-1,1]);
plot2d(1:n,Dt.*exp(Xt_Const),1:3,"061","",[0,0,20,20]);
xtitle('Evolution du prix sopt');
end
  

// 1: autocorrelation des residus, calage et test de box-pierce
if Methode == 1
alpha = 0.05;
X0 = 0;

Xt_Const = Simul_Xt_Const(n,h,mu,lambda,sigma);
Param_EMV = Calage_Xt_Const(n,h,Xt_Const,X0);
disp(Param_EMV,"Valeurs estime par EMV");


//3: autocorrelation des residus et test de box-pierce
alpha = 0.05;
[Resid, Autocorr, BP_Test, Probac, JB_Test, S, K, Proban] = Residu_Const(n,h,Param_EMV,Xt_Const,X0,alpha);


xbasc(1);
xset("window",1);
xsetech([0,1/4,0.5,1/2],[-1,1,-1,1]);
histplot(100,Resid,1:3,"061","",[0,0,20,20]);
xtitle('Histogramme des residus');

xsetech([0.5,1/4,0.5,1/2],[-1,1,-1,1]);
plot2d(1:length(Autocorr),Autocorr,1:3,"061","",[0,0,20,20]);
xtitle('La fonction d autocorrelation en fonction de retard');


disp(Probac,'La probabilite d absence de autocorrelation');
disp(BP_Test,'Test Box-Pierce 1:absence d autocorrelation, 0:autocorrelation');


disp(Proban,'La probabilite de normalite');
disp(JB_Test,'Test Jarque Bera :normalite des residus, 0:rejet');

disp([S K],'Le skewness et le kurtosis');
end
