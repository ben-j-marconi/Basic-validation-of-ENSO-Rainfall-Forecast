%KT Lab 09
Precip= load('CRUPrcp.tsv');
X =load('pc1_Kaplan_ssta.txt');

%Probability 1.1
long = 117.5:5:302.5;
lat = -32.5:5:32.5;
years = 1950:2009;
nx=38;
ny=14;
nt=60;

%remove NaNs
iocean=find(isfinite(X(1,:)));
Xo=X(:,iocean);

%Covariance matrix 
Cova = cov(Xo);
[EOFs,ve] = eig(Cova);
PCs=Xo*EOFs;
Var=real( diag(ve./trace(ve)) );
EOFLand=nan(nx*ny,2);

%check VAR for desc/asc order?  
%Var is in ascending order, use:
EOFLand(iocean,1)=EOFs(:,end);
PC1 = PCs(:,end);
eof1_xy = reshape(EOFLand(:,1),nx,ny);
z= (Precip-mean(Precip)/std(Precip));
AvePrecip= mean(Precip);
zlow=(424-AvePrecip)/std(Precip);
zhigh=(514-AvePrecip)/std(Precip);
PL = normcdf(zlow);
PH= 1-normcdf(zhigh);
P = 1-(PH-PL);

%Probability
years_below=length(find(Precip<424));
years_above=length(find(Precip>514));
year_average=60-years_below-years_above;
pb=years_below/60
pa=years_above/60
pave=(year_average)/60;

%correlation=corrcoef(Precip, PC1)
xbar= mean(Precip);
ybar= mean(PC1);
sd_Precip=std(Precip);
sd_PC1=std(PC1);
zPC1= (PC1-ybar)/(sd_PC1)
R=corrcoef([ z' zPC1 ]);
corr=R(1,2)
zy=corr*zPC1;

%STANDARD OF ERROR 
SEofE=sqrt(1-(corr).^2);

% EMV 1.2
X0=[0 0 0 0]';
A=[1 1 1 1];
b=[43934];
Aeq=[];
beq=[];
lb=[0;0;0;0];
ub=[43934;43934;43934;43934];
Q=fmincon(@profit,X0,A,b,Aeq,beq,lb,ub);
bnor = [126273389.2;281290981.6;274513384.5];
nor = [82788562.09;298080293;379955737.2];
anor = [18035935.59;250521893.4;450323241.4];
bnorEMV = (126273389.2*0.35)+(281290981.6*0.35)+(274513384.5*0.3);
EMV = (82788562.09*0.35)+(298080293*0.35)+(379955737.2*0.3);
anorEMV = (18035935.59*0.35)+(250521893.4*0.35)+(450323241.4*0.3);

%Perfect Forecast 
Perfect = ((450323241.4*0.3)+(298080293*0.35)+(126273389.2*0.35));
PerfectVal= Perfect -(247290000*0.35);

% 2.1 Regression
%correlation=corrcoef(Precip, PC1)
PC1= load('pc1_Kaplan_ssta.txt');
xbar= mean(Precip);
ybar= mean(PC1);
sd_Precip=std(Precip);
sd_PC1=std(PC1);
zPC1= (PC1-ybar)/(sd_PC1)
R=corrcoef([ z' zPC1 ]);
corr=R(1,2)
zy=corr*zPC1;

% Part 2 - probabalistic loop
probability=zeros(60,3);
for t=1:60
    y_Precip=0.4375*zPC1(t);
    P_lowA=normcdf(-0.431, y_Precip, 0.8992);
    P_highA=1-normcdf(0.431, y_Precip, 0.8992);
    P_normal=1-P_lowA-P_highA;
    Probability(t,3)=P_lowA;
    Probability(t,2)=P_normal;
    Probability(t,1)=P_highA;
end

value=[450323241.4 250521893.4 18035935.59;379955737.2 298080293 82788562.09;274513384.5 281290981.6 126273389.2]
% EMV
for t=1:60
 for d=1:3
     EMV(t,d)=Probability(t,:)*value(d,:)';
 end
end

% 2.2 decision
decision=zeros(60,1)
for t=1:60
    [a b]=max(EMV(t,:));
    decision(t)=b;
end

decisionBN=decision(find(Precip<424));
decisionN=decision(find(Precip>424 & Precip<514));
decisionAN=decision(find(Precip>514));