clear;

chdir('F:\My Tempo\Simulations')

getf('Modele facteur HG\HG_Function.sci');disp('getf done');

NIG_DENSITE =%F;
VG_DENSITE = %F;
NIG_RESID_SIMULATION = %F;
VG_RESID_SIMULATION = %F;
NIG_OPTION_MC = %F;
VG_OPTION_MC = %F;
NIG_OPTION_FERMEE = %F;
VG_OPTION_FERMEE = %F;
SIMULATION_PRIXSPOT_NIG = %F;
SIMULATION_PRIXSPOT_VG = %F;



if(NIG_DENSITE)
//Densite du NIG
alpha = 0.33;
beta = 0;
delta = 0.8;
mu = 0 ;

n = 200;
m = 10;

X1 = zeros(1,n+1);
X2 = zeros(1,n+1);
X3 = zeros(1,n+1);
Y = zeros(1,n+1);

Y = -m+((0:n)*2*m)/n;
X1 = Densite_NIG(Y,alpha,beta,delta,mu);
X2 = Densite_NIG(Y,alpha,beta+0.32,delta,mu);
X3 = Densite_NIG(Y,alpha,beta,delta,2);


//0: Dessin de trois graphes ft, st, xt
xbasc(0);
xset("window",0);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,exp(-Y^2/2)/(sqrt(2*%pi)),1:3,"061","",[0,0,20,20]);
plot2d(Y,X1,[-2,0,0]);
//xtitle('Evolution du forward F(0,t) Day-Ahead');
xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X2,1:3,"061","",[0,0,20,20]);
plot2d(Y,X1,[-2,0,0]);
//xtitle('Evolution du prix spot St Day-Ahead');
xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X3,1:3,"061","",[0,0,20,20]);
plot2d(Y,X1,[-2,0,0]);
//xtitle('Evolution du log-prix spot Xt=ln(St/S0) Day-Ahead')

end




if(VG_DENSITE)
//Densite du VG
alpha = 0.33;
beta = 0;
lambda = 0.6;
mu = 0.23 ;

n = 200;
m = 10;

X1 = zeros(1,n+1);
X2 = zeros(1,n+1);
X3 = zeros(1,n+1);
Y = zeros(1,n+1);

Y = -m+((0:n)*2*m)/n;
X1 = Densite_VG(Y,lambda,alpha,beta,mu);
X2 = Densite_VG(Y,lambda,alpha-0.2,beta,mu);
X3 = Densite_VG(Y,lambda,alpha,beta,2.01);

//0: Dessin de trois graphes ft, st, xt
xbasc(1);
xset("window",1);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,exp(-Y^2/2)/(sqrt(2*%pi)),1:3,"061","",[0,0,20,20]);
plot2d(Y,X1,[-2,0,0]);
//xtitle('Evolution du forward F(0,t) Day-Ahead');
xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X2,1:3,"061","",[0,0,20,20]);
plot2d(Y,X1,[-2,0,0]);
//xtitle('Evolution du prix spot St Day-Ahead');
xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X3,1:3,"061","",[0,0,20,20]);
plot2d(Y,X1,[-2,0,0]);
//xtitle('Evolution du log-prix spot Xt=ln(St/S0) Day-Ahead')
end




if(NIG_RESID_SIMULATION)
n = 500;
alpha = 0.33;
beta = 0;
delta = 1;
mu = 0;
omega = 0.01;


X1 = zeros(1,n);
X2 = zeros(1,n);
X3 = zeros(1,n);

X1 = Var_NIG(n,alpha,beta,delta,mu);
X2 = Var_NIG(n,alpha,beta+0.3,delta,mu);
X3 = Var_NIG(n,alpha/10,beta,delta,mu);
Y = 1:n;

xbasc(2);
xset("window",2);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X1,1:3,"061","",[0,0,20,20]);

xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X2,1:3,"061","",[0,0,20,20]);

xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X3,1:3,"061","",[0,0,20,20]);

end


if(VG_RESID_SIMULATION)
n = 1000;
lambda = 0.7;
alpha = 0.33;
beta = 0;
delta = 1;
mu = 0;



X1 = zeros(1,n);
X2 = zeros(1,n);
X3 = zeros(1,n);


X1 = abs(Var_VG(n,lambda,alpha,beta,mu));
X2 = abs(Var_VG(n,lambda,alpha,(beta+0.3),mu));
X3 = abs(Var_VG(n,lambda,alpha,(beta-0.3),mu));

Y = 1:n;

xbasc(2);
xset("window",2);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X1,1:3,"061","",[0,0,20,20]);

//xtitle('Evolution du forward F(0,t) Day-Ahead');
xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X2,1:3,"061","",[0,0,20,20]);

//xtitle('Evolution du prix spot St Day-Ahead');
xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,X3,1:3,"061","",[0,0,20,20]);

end


if(SIMULATION_PRIXSPOT_NIG)

h = 1/(24*365);//horaire
alpha = 2;
beta = 0;
delta = 2;
mu = 0;
a = 100;
sigma = 1.9;
omega = 0;

//Donnees d'entree
ForwardFile = "Data\EDF_Forward_01042004_08052004.txt";//EDF_Forward_01012002_31122002.txt";

//Matrice du prix spot et forward
MF = read(ForwardFile,-1,1)';

n = length(MF);
Xt = zeros(1,n);

Xt1 = Simul_NIG_Process(n,h,a,sigma,omega,alpha,beta,delta,mu);

Y = 1:n;

xbasc(1);
xset("window",2);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,Xt1,1:3,"061","",[0,0,20,20]);
xtitle('Evolution du logarithme du prix -Modele NIG-');

xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,MF,1:3,"061","",[0,0,20,20]);
xtitle('Evolution de la courbe forward');

xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,MF.*exp(Xt1),1:3,"061","",[0,0,20,20]);
xtitle('Evolution du prix sopt');
end



if(SIMULATION_PRIXSPOT_VG)

h = 1/(24*365);//horaire
alpha = 2;
beta = 0;
lambda = 1;
mu = 0;
a = 100;
sigma = 1.9;
omega = 0;

//Donnees d'entree
ForwardFile = "Data\EDF_Forward_01042004_08052004.txt";//EDF_Forward_01012002_31122002.txt";

//Matrice du prix spot et forward
MF = read(ForwardFile,-1,1)';

n = length(MF);

Xt = zeros(1,n);
Xt2 = Simul_VG_Process(n,h,a,sigma,omega,lambda,alpha,beta,mu);
Y = 1:n;

xbasc(2);
xset("window",2);
xsetech([0,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,Xt2,1:3,"061","",[0,0,20,20]);
xtitle('Evolution du logarithme du prix -Modele VG-');

xsetech([1/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,MF,1:3,"061","",[0,0,20,20]);
xtitle('Evolution de la courbe forward');

xsetech([2/3,1/4,1/3,1/2],[-1,1,-1,1]);
plot2d(Y,MF.*exp(Xt2),1:3,"061","",[0,0,20,20]);
xtitle('Evolution du prix sopt');
end


if(NIG_OPTION_MC)
n = 10000;
m = 10000;
l = 5000;
alpha = 0.3;
beta = 0.001;
delta = 0.3;
mu = 10;
a = 0.19;
sigma = 0.13;
omega = 0;

Tf = 1+30/365;
t = 0;
T0 = 1;

C1 = zeros(1,m);
P1 = zeros(1,m);
C2 = zeros(1,m);
P2 = zeros(1,m);

x = 100;
r = 0.25;
K = 100;

XT = 0;
rhoc = 0;
rhop = 0;

XT = Simul_XT_NIG(l,n,T0,Tf,t,alpha,beta+omega,delta,mu,a,sigma);
rhoc = mean(exp(-r*(T0-t))*max(x*exp(XT)-K,0).*(1-exp(XT)));
rhop = mean(exp(-r*(T0-t))*max(K-x*exp(XT),0).*(1-exp(XT)));

varXT = exp(integrate('Phi_NIG(2*sigma*exp(-a*(Tf-y)),alpha,beta+omega,delta,mu)-2*Phi_NIG(sigma*exp(-a*(Tf-y)),alpha,beta+omega,delta,mu)','y',t,T0))-1;

rhoc = rhoc/varXT;
rhop = rhop/varXT;


disp([rhoc rhop],"rhoc rhop");

XT = Simul_XT_NIG(m,n,T0,Tf,t,alpha,beta+omega,delta,mu,a,sigma);
C1 = exp(-r*(T0-t))*(max(x*exp(XT)-K,0)+rhoc*exp(XT)-rhoc);
P1 = exp(-r*(T0-t))*(max(-x*exp(XT)+K,0)+rhop*exp(XT)-rhop);
C2 = exp(-r*(T0-t))*(max(x*exp(XT)-K,0));
P2 = exp(-r*(T0-t))*(max(-x*exp(XT)+K,0));

disp([mean(C1) mean((C1-mean(C1))^2)],"prix et variance du call Avec reduction");
disp([mean(P1) mean((P1-mean(P1))^2)],"prix et variance du put Avec reduction");
disp([mean(C2) mean((C2-mean(C2))^2)],"prix et variance du call Sans reduction");
disp([mean(P2) mean((P2-mean(P2))^2)],"prix et variance du put Sans reduction");
end



if(VG_OPTION_MC)
n = 10000;
m = 500000;
l = 5000;
lambda = 1;
alpha = 0.3;
beta = 0.001;
mu = 10;
a = 0.19;
sigma = 0.13;
omega = 0;

Tf = 1+30/365;
t = 0;
T0 = 1;

C1 = zeros(1,m);
P1 = zeros(1,m);
C2 = zeros(1,m);
P2 = zeros(1,m);

x = 100;
r = 0.25;
K = 100;

XT = 0;
rhoc = 0;
rhop = 0;

XT = Simul_XT_VG(l,n,T0,Tf,t,lambda,alpha,beta+omega,mu,a,sigma);
rhoc = mean(exp(-r*(T0-t))*max(x*exp(XT)-K,0).*(1-exp(XT)));
rhop = mean(exp(-r*(T0-t))*max(K-x*exp(XT),0).*(1-exp(XT)));

varXT = exp(integrate('Phi_VG(2*sigma*exp(-a*(Tf-y)),lambda,alpha,beta+omega,mu)-2*Phi_VG(sigma*exp(-a*(Tf-y)),lambda,alpha,beta+omega,mu)','y',t,T0))-1;

rhoc = rhoc/varXT;
rhop = rhop/varXT;


disp([rhoc rhop],"rhoc rhop");

XT = Simul_XT_VG(m,n,T0,Tf,t,lambda,alpha,beta+omega,mu,a,sigma);
C1 = exp(-r*(T0-t))*(max(x*exp(XT)-K,0)+rhoc*exp(XT)-rhoc);
P1 = exp(-r*(T0-t))*(max(-x*exp(XT)+K,0)+rhop*exp(XT)-rhop);
C2 = exp(-r*(T0-t))*(max(x*exp(XT)-K,0));
P2 = exp(-r*(T0-t))*(max(-x*exp(XT)+K,0));

disp([mean(C1) mean((C1-mean(C1))^2)],"prix et variance du call Avec reduction");
disp([mean(P1) mean((P1-mean(P1))^2)],"prix et variance du put Avec reduction");
disp([mean(C2) mean((C2-mean(C2))^2)],"prix et variance du call Sans reduction");
disp([mean(P2) mean((P2-mean(P2))^2)],"prix et variance du put Sans reduction");
end



if(NIG_OPTION_FERMEE)

alpha = 1;
beta = 0;
delta = 1;
mu = 0;
a = 0.19;
sigma = 0.13;
omega = 0;

Tf = 1+30/365;
t = 0;
T0 = 1;

x = 100;
r = 0.25;
K = 100;

eta = 0.25;
epsilon = 10^(-3);//10^(-3);
h = 0.001;//2.5*10^(-4);

C = Call_NIG_Fermee(r,K,x,eta,T0,Tf,t,alpha,beta+omega,delta,mu,a,sigma,epsilon, h);
disp([C h],"price call euro NIG");

end


if(VG_OPTION_FERMEE)

lambda = 1;
alpha = 0.3;
beta = 0.001;
mu = 10;
a = 0.19;
sigma = 0.13;
omega = 0;

Tf = 1+30/365;
t = 0;
T0 = 1;

x = 100;
r = 0.25;
K = 100;

eta = 0.25;
epsilon = 10^(-2);
h = 0.005;

C = Call_VG_Fermee(h,r,K,x,eta,T0,Tf,t,lambda,alpha,beta+omega,mu,a,sigma,epsilon);
disp([C h],"price call euro VG");

end



