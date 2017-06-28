clear;

chdir('J:\Stage\Simulations')

getf('Modele un facteur gaussien\Cas Constant\Const_Function.sci');disp('getf done');

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Tendance Constante///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

Calage_Day = %F;
Calage_Heure = %T;


//Calage pour la moyenne horaire -DayAhead-
if(Calage_Day)

h = 1/(365);
n = 0;

//Donnees d'entree
SpotFile = "Data\EDF_Spot_ 01082003_15062004_Des.txt";//EDF_Spot_01012002_31122002.txt//

//Matrice du prix spot et forward
MSHT = read(SpotFile,-1,1)';

nTHT = length(MSHT);

nTH = nTHT - (0)*24;

MSH = zeros(1,nTH);

MSH = MSHT(((0*24)+1):nTHT);

nT = nTH/24;
MS =zeros(1,nT);

for i = 1:nT
MS(1,i) = mean(MSH((1+(i-1)*24):(i*24)));
end


//Nombre d observations
nT = length(MS);
n = nT-1;

St_obs = zeros(1,n);
Xt_obs = zeros(1,n);

S0 = MS(1);
X0 = 0;

//Le processus log-prix et forward
St_obs = MS(2:nT);

Xt_obs = log(St_obs/S0);



//Informations
disp(h,'le pas de temps h');
disp(nT,'le nombre d observation Jours');



//0: Dessin de trois graphes ft, st, xt
xbasc(0);
xset("window",0);
xsetech([0,1/4,1/2,1/2],[-1,1,-1,1]);
plot2d(1:nT,{S0,St_obs},1:3,"061","",[0,0,20,20]);
xtitle('Evolution du prix spot St');
xsetech([1/2,1/4,1/2,1/2],[-1,1,-1,1]);
plot2d(1:nT,{X0,Xt_obs},1:3,"061","",[0,0,20,20]);
xtitle('Evolution du log-prix spot Xt=ln(St/S0)')


//1: Calage par MCO
Theta_MCO = MCO_Const(n,h,Xt_obs,X0);
disp(Theta_MCO,'Les parametres : mu, lambda, sigma par MCO sur Day-Ahead');


//2: Calage par EMV
Theta_EMV = Calage_Xt_Const(n,h,Xt_obs,X0);
disp(Theta_EMV,'Les parametres : mu, lambda, sigma par EMV constant sur Day-Ahead');



//3: autocorrelation des residus et test de box-pierce
alpha = 0.05;
[Resid, Autocorr, BP_Test, Probac, JB_Test, S, K, Proban] = Residu_Const(n,h,Theta_EMV,Xt_obs,X0,alpha);


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




//Calage heure par heure pour une heure journaliere donnee
if(Calage_Heure)

h = 1/(365);
n = 0;

//Donnees d'entree
SpotFile = "Data\EDF_Spot_01042004_08052004.txt";
ForwardFile = "Data\EDF_Forward_01042004_08052004.txt";

//Matrice du prix spot et forward
MSH = read(SpotFile,-1,1)';


nTH = length(MSH);

nT = nTH/24;
MSJ =zeros(24,nT);

for i = 1:nT
for j = 1:24
MSJ(j,i) = MSH(j+24*(i-1));
end
end

//Informations
disp(h,'le pas de temps h');
disp(nT,'le nombre d observation Heures');

EMV_Vect = zeros(3,24);

for k=1:24

MS = MSJ(k,:);

//Nombre d observations
nT = length(MS);
n = nT-1;

St_obs = zeros(1,n);
Xt_obs = zeros(1,n);

S0 = MS(1);
X0 = 0;

//Le processus log-prix et forward
St_obs = MS(2:nT);

Xt_obs = log(St_obs/S0);


//0: Dessin de trois graphes ft, st, xt
xbasc(0);
xset("window",0);
xsetech([0,1/4,1/2,1/2],[-1,1,-1,1]);
plot2d(1:nT,{S0,St_obs},1:3,"061","",[0,0,20,20]);
xtitle('Evolution du prix spot St');
xsetech([1/2,1/4,1/2,1/2],[-1,1,-1,1]);
plot2d(1:nT,{X0,Xt_obs},1:3,"061","",[0,0,20,20]);
xtitle('Evolution du log-prix spot Xt=ln(St/S0)')


//2: Calage par EMV
Theta_EMV_Const = Calage_Xt_Const(n,h,Xt_obs,X0);
disp(Theta_EMV_Const,'Les parametres : lambda, sigma par EMV variable Heure');


//3: autocorrelation des residus et test de box-pierce
alpha = 0.05;
[Resid, Autocorr, BP_Test, Probac, JB_Test, S, K, Proban] = Residu_Const(n,h,Theta_EMV_Const,Xt_obs,X0,alpha);

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
disp(JB_Test,'Test Jarque Bera :normalite des residus, 0:rejet');

disp([S K],'Le skewness et le kurtosis');

EMV_Vect(:,k) = Theta_EMV_Const;

end

xbasc(2);
xset("window",2);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(1:24,EMV_Vect(1,:),1:3,"061","",[0,0,20,20]);
xtitle('Evolution du mu dans un jour');
xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(1:24,EMV_Vect(2,:),1:3,"061","",[0,0,20,20]);
xtitle('Evolution du lambda dans un jour');
xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(1:24,EMV_Vect(3,:),1:3,"061","",[0,0,20,20]);
xtitle('Evolution du sigma dans un jour')

end


