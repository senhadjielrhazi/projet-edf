clear;

chdir('J:\Stage\Simulations')
getf('Modele un facteur gaussien\Cas Variable\Var_Function.sci');disp('getf done');
getf('Modele un facteur gaussien\Cas Constant\Const_Function.sci');disp('getf done');



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tendance Variable////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
// donnees
////////////////////////////////////////////////////////
// lnf(0,t) le logarithme de la courbe forward en (0,t)
// sigma = [100%,400%] annuelle
// h =1/(24*365)  heure, 1/365 jour
// lambda = [1,25] jours ou [10,200] annee
///////////////////////////////////////////////////////

//simulation heure par heure
lambda = 100;
h = 1/(24*365);
sigma = 4;

//simulation jour par jour
//lambda = 10;
//h = 1/365;
//sigma = 1*sqrt(h); sigma_annuelle=1;
//n = 365*3;

		
//Type de donnees a trouver :
Methode = 1;  // 0: simulation du logarithme et du spot  
	            // 1: autocorrelation des residus, calage et test de box-pierce, 
	            // 4: Une simulations et un calage du Log-Prix par MCO-EMV,
	            // 5: Applications au cas de EDF 2004-2007


ForwardFile = "Data\EDF_Forward_01012002_31122002.txt";
MF = log((read(ForwardFile,-1,1))');
n  = length(MF)-1;
lnForward = MF(2:(n+1));
lnF0 = MF(1);


disp(h,'le pas de temps h');
disp(n,'le nombre d observation');
disp([lambda sigma],'vrais valeurs de lambda et sigma');


// 0: simulation et dessin 3 graphes pour 3 prix forward differents
if Methode == 0
Xt_Var = zeros(1,n); 

Xt_Var = Simul_Xt_Var(n,h,lambda,sigma,lnForward,lnF0);

xbasc(0);
xset("window",0);

xsetech([0,1/4,1/2,1/2],[-1,1,-1,1]);
plot2d(1:n,Xt_Var,1:3,"061","",[0,0,20,20]);
xtitle('Xt le logarithme du prix spot');

xsetech([1/2,1/4,1/2,1/2],[-1,1,-1,1]);
plot2d(1:n,exp(Xt_Var),1:3,"061","",[0,0,20,20]);
xtitle('St pour une tendance EDF 2004-2007');

end


// 1: autocorrelation des residus, calage et test de box-pierce
if Methode == 1
alpha = 0.05;
Xt_Var = Simul_Xt_Var(n,h,lambda,sigma,lnForward,lnF0);
X0 = 0;

Param_EMV = Calage_Xt_Var(n,h,Xt_Var,X0,lnForward,lnF0);
disp(Param_EMV,"valeurs estime des parametres");

[Resid, Autocorr, BP_Test, Probac, JB_Test, S, K, Proban] = Residu_Var(n,h,Param_EMV,Xt_Var,X0,lnForward,lnF0,alpha);

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
disp(JB_Test,'Test Jarque Bera 1:normalite des residus, 0:rejet');

disp([S K],'Le skewness et le kurtosis');
end




// 4: Une simulations et un calage du Log-Prix par EMV
if Methode == 4
Xt = zeros(1,n+1); 
F = zeros(1,n);
F(2:n) = lnForward(1:n-1);
I = cumsum(ones(1,n));
Tend = zeros(1,n);

Tend = lnForward + (1/(lambda*h))*(lnForward-F)-(sigma^2/(4*lambda))*(1+exp(-2*lambda*I*h));
Xt_Var = Simul_Xt_Var(n,h,lambda,sigma,lnForward);
Xt(2:n+1) = Xt_Var(1:n);

xbasc(0);
xset("window",0);
plot2d(1:n+1,Xt,1:3,"061","",[0,0,20,20]);
plot2d(1:n,Tend,[2,2,4]);
xtitle('Evolution du log-prix spot en fonction du temps en bleu la tendance');

Param_EMV = Calage_Xt_Var(n,h,Xt_Var,lnForward);

disp(Param_EMV,'Valeurs des parametres avec methode EMV');
end


// 5: Applications au cas de EDF 2004-2007
if Methode == 6

xbasc(0);
xset("window",0);
plot2d(exp(lnForward));
xtitle('La courbe forward initiale');

Xt_Var = Simul_Xt_Var(n,h,lambda,sigma,lnForward);

xbasc(1);
xset("window",1);
plot2d(Xt_Var);
xtitle('Le log-prix spot 2004-2007');

xbasc(2);
xset("window",2);
plot2d(exp(Xt_Var));
xtitle('Le prix spot 2004-2007');

Param = Calage_Xt_Var(n,h,Xt_Var,lnForward);

disp(n,'nombre de simulations n');

disp([Param(1) Param(2)],'vrais calage de lambda et sigma')
end

