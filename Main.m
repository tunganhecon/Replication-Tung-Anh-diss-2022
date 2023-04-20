%Replication codes of Assessing the causes and effects of oil price shocks: Empirical evidence from Vietnam
%Author: Lutz Killian (1997)
%Updated by Tung Anh Nguyen

%% Preparation
clc 
clear
close all

%% Replicating figure 12
clear;

trivar

% Compute structural shocks Ehat
Ehat=inv(chol(SIGMA)')*Uhat(1:q,:);
q1=Ehat(1,:); q1=[(q1(1,1)+q1(1,2))/2 q1];
q2=Ehat(2,:); q2=[(q2(1,1)+q2(1,2))/2 q2];
q3=Ehat(3,:); q3=[(q3(1,1)+q3(1,2))/2 q3];

% Average monthly shocks by quarter
time=1976:1:2021;
for i=1:length(time)
   q1a(i)=(q1(12*(i-1)+1)+q1(12*(i-1)+2)+q1(12*(i-1)+3)+q1(12*(i-1)+4)+q1(12*(i-1)+5)+q1(12*(i-1)+6)+q1(12*(i-1)+7)+q1(12*(i-1)+8)+q1(12*(i-1)+9)+q1(12*(i-1)+10)+q1(12*(i-1)+11)+q1(12*(i-1)+12))/12;
   q2a(i)=(q2(12*(i-1)+1)+q2(12*(i-1)+2)+q2(12*(i-1)+3)+q2(12*(i-1)+4)+q2(12*(i-1)+5)+q2(12*(i-1)+6)+q2(12*(i-1)+7)+q2(12*(i-1)+8)+q2(12*(i-1)+9)+q2(12*(i-1)+10)+q2(12*(i-1)+11)+q2(12*(i-1)+12))/12;
   q3a(i)=(q3(12*(i-1)+1)+q3(12*(i-1)+2)+q3(12*(i-1)+3)+q3(12*(i-1)+4)+q3(12*(i-1)+5)+q3(12*(i-1)+6)+q3(12*(i-1)+7)+q3(12*(i-1)+8)+q3(12*(i-1)+9)+q3(12*(i-1)+10)+q3(12*(i-1)+11)+q3(12*(i-1)+12))/12;
end;

subplot(3,1,1)
plot(time,q1a,time,zeros(size(q1a))); axis tight;
title('(a) Oil Supply Shock')
axis([1976 2021 -1 +1])
%axis([1976 1983 -1 +1])

subplot(3,1,2)
plot(time,q2a,time,zeros(size(q2a))); axis tight;
title('(b) Global Demand Shock')
axis([1976 2021 -1 +1])
%axis([1976 1983 -1 +1])

subplot(3,1,3)
plot(time,q3a,time,zeros(size(q3a))); axis tight;
title('(c) Oil-Specific Demand Shock')
axis([1976 2021 -1 +1])
%axis([1976 1983 -1 +1])

%% Replicating figure 13
clear; 
trivar
% VAR Impulse response analysis 
[IRF]=irfvar(A,SIGMA(1:q,1:q),p,h);
    IRF(1,:)=cumsum(IRF(1,:));
    IRF(4,:)=cumsum(IRF(4,:));
    IRF(7,:)=cumsum(IRF(7,:));

% VAR bootstrap
randn('seed',1234);
nrep=3000;
IRFmat=zeros(nrep,(q^2)*(h+1));

[t,q]=size(y);				
y=y';
Y=y(:,p:t);	
for i=1:p-1
 	Y=[Y; y(:,p-i:t-i)];		
end;

Ur=zeros(q*p,t-p);   
Yr=zeros(q*p,t-p+1); 
U=Uhat;    
for r=1:nrep
    r    
	pos=fix(rand(1,1)*(t-p+1))+1;
	Yr(:,1)=Y(:,pos);

    % Recursive design WB
    eta=randn(1,size(Uhat,2));  eta=[eta; eta; eta];
	Ur(1:q,2:t-p+1)=U(1:q,:).*eta;	
    
	for i=2:t-p+1
		Yr(:,i)= V + A*Yr(:,i-1)+Ur(:,i); 
	end;

	yr=[Yr(1:q,:)];
	for i=2:p
		yr=[Yr((i-1)*q+1:i*q,1) yr];
    end;
    yr=yr';
    [Ar,SIGMAr]=olsvarc(yr,p);

    % Compute IRFs
    IRFr=irfvar(Ar,SIGMAr(1:q,1:q),p,h);
    IRFr(1,:)=cumsum(IRFr(1,:));
    IRFr(4,:)=cumsum(IRFr(4,:));
    IRFr(7,:)=cumsum(IRFr(7,:));
    IRFmat(r,:)=vec(IRFr)';
end;
IRFrstd=reshape((std(IRFmat)'),q^2,h+1);
CI1LO=IRF-1*IRFrstd; CI1UP=IRF+1*IRFrstd;
CI2LO=IRF-2*IRFrstd; CI2UP=IRF+2*IRFrstd;

horizon=0:h;

subplot(3,3,1)
plot(horizon,-IRF(1,:),'r-',horizon,-CI1LO(1,:),'b--',horizon,-CI1UP(1,:),'b--',horizon,-CI2LO(1,:),'b:',horizon,-CI2UP(1,:),'b:',horizon,zeros(size(horizon)));
title('(a) Oil supply shock')
ylabel('Oil production (Δprod) (%)')
axis([0 h -2 1])

subplot(3,3,2)
plot(horizon,-IRF(2,:),'r-',horizon,-CI1LO(2,:),'b--',horizon,-CI1UP(2,:),'b--',horizon,-CI2LO(2,:),'b:',horizon,-CI2UP(2,:),'b:',horizon,zeros(size(horizon)));
title('(b) Oil supply shock')
ylabel('Real activity (rea) (%)')
axis([0 h -7 22])

subplot(3,3,3)
plot(horizon,-IRF(3,:),'r-',horizon,-CI1LO(3,:),'b--',horizon,-CI1UP(3,:),'b--',horizon,-CI2LO(3,:),'b:',horizon,-CI2UP(3,:),'b:',horizon,zeros(size(horizon)));
title('(c) Oil supply shock')
ylabel('Real oil price (rpo) (%)')
axis([0 h -7 17])

subplot(3,3,4)
plot(horizon,IRF(4,:),'r-',horizon,CI1LO(4,:),'b--',horizon,CI1UP(4,:),'b--',horizon,CI2LO(4,:),'b:',horizon,CI2UP(4,:),'b:',horizon,zeros(size(horizon)));
title('(d) Global demand shock')
ylabel('Oil Production (Δprod) (%)')
axis([0 h -2 1])

subplot(3,3,5)
plot(horizon,IRF(5,:),'r-',horizon,CI1LO(5,:),'b--',horizon,CI1UP(5,:),'b--',horizon,CI2LO(5,:),'b:',horizon,CI2UP(5,:),'b:',horizon,zeros(size(horizon)));
title('(e) Global demand shock')
ylabel('Real activity (rea) (%)')
axis([0 h -7 22])

subplot(3,3,6)
plot(horizon,IRF(6,:),'r-',horizon,CI1LO(6,:),'b--',horizon,CI1UP(6,:),'b--',horizon,CI2LO(6,:),'b:',horizon,CI2UP(6,:),'b:',horizon,zeros(size(horizon)));
title('(f) Global demand shock')
ylabel('Real oil price (rpo) (%)')
axis([0 h -7 17])

subplot(3,3,7)
plot(horizon,IRF(7,:),'r-',horizon,CI1LO(7,:),'b--',horizon,CI1UP(7,:),'b--',horizon,CI2LO(7,:),'b:',horizon,CI2UP(7,:),'b:',horizon,zeros(size(horizon)));
title('(g) Oil-specific demand shock')
ylabel('Oil production (Δprod) (%)')
xlabel('Months')
axis([0 h -2 1])

subplot(3,3,8)
plot(horizon,IRF(8,:),'r-',horizon,CI1LO(8,:),'b--',horizon,CI1UP(8,:),'b--',horizon,CI2LO(8,:),'b:',horizon,CI2UP(8,:),'b:',horizon,zeros(size(horizon)));
title('(h) Oil-specific demand shock')
ylabel('Real activity (rea) (%)')
xlabel('Months')
axis([0 h -7 22])

subplot(3,3,9)
plot(horizon,IRF(9,:),'r-',horizon,CI1LO(9,:),'b--',horizon,CI1UP(9,:),'b--',horizon,CI2LO(9,:),'b:',horizon,CI2UP(9,:),'b:',horizon,zeros(size(horizon)));
title('(i) Oil-specific demand shock')
ylabel('Real oil price (rpo) (%)')
xlabel('Months')
axis([0 h -7 17])

%% Replicating figure 14
clear;
trivar

% Compute structural multipliers
IRF=irfvar(A,SIGMA(1:q,1:q),p,t-p-1);

% Compute structural shocks Ehat from reduced form shocks Uhat
Ehat=inv(chol(SIGMA)')*Uhat(1:q,:);

% Cross-multiply the weights for the effect of a given shock on the real
% oil price (given by the relevant row of IRF) with the structural shock
% in question
yhat1=zeros(t-p,1); yhat2=zeros(t-p,1); yhat3=zeros(t-p,1);
for i=1:t-p
	yhat1(i,:)=dot(IRF(3,1:i),Ehat(1,i:-1:1));
	yhat2(i,:)=dot(IRF(6,1:i),Ehat(2,i:-1:1));
	yhat3(i,:)=dot(IRF(9,1:i),Ehat(3,i:-1:1));
end;

% Contribution to real price of oil
time=1976+1/12:1/12:2021+11/12;

subplot(3,1,1)
plot(time,yhat1,'b-');
title('(a) Cumulative Effect of Oil Supply Shock on Real Oil Price')
axis([1976+1/12 2021+10/12 -100 +100])
grid on

subplot(3,1,2)
plot(time,yhat2,'b-');
title('(b) Cumulative Effect of Global Demand Shock on Real Oil Price')
axis([1976+1/12 2021+10/12 -100 +100])
grid on

subplot(3,1,3)
plot(time,yhat3,'b-');
title('(c) Cumulative Effect of Oil-Specific Demand Shock on Real Oil Price')
axis([1976+1/12 2021+10/12 -100 +100])
grid on

%% Replicating figure 15
clear;

trivar

Ehat=inv(chol(SIGMA)')*Uhat(1:q,:);
q1=Ehat(1,:); q1=[(q1(1,1)+q1(1,2))/2 q1];
q2=Ehat(2,:); q2=[(q2(1,1)+q2(1,2))/2 q2];
q3=Ehat(3,:); q3=[(q3(1,1)+q3(1,2))/2 q3];

% Average monthly shocks by quarter
time=2001+1/4:1/4:2021+4/4;
for i=1:length(time)
   q1q(i)=(q1(3*(i-1)+1)+q1(3*(i-1)+2)+q1(3*(i-1)+3))/3;
   q2q(i)=(q2(3*(i-1)+1)+q2(3*(i-1)+2)+q2(3*(i-1)+3))/3;
   q3q(i)=(q3(3*(i-1)+1)+q3(3*(i-1)+2)+q3(3*(i-1)+3))/3;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=(2001+1/4:1/4:2021+4/4)';

VNdata = xlsread('VNdata.xls'); 
% IRF for GDP
y = VNdata(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0:12;

[irf2hat,irf2a,irf2b,cumirf2hat,cumirf2a,cumirf2b]=stage2irf(y,q1q');
subplot(3,3,1)
plot(time,-cumirf2hat,'b-',time,-cumirf2a,'b--',time,-cumirf2b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -45 40])
ylabel('Real GDP')
title('(a) Oil Supply Shock')

[irf3hat,irf3a,irf3b,cumirf3hat,cumirf3a,cumirf3b]=stage2irf(y,q2q');
subplot(3,3,4)
plot(time,cumirf3hat,'b-',time,cumirf3a,'b--',time,cumirf3b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -45 40])
ylabel('Real GDP')
title('(b) Global Demand Shock')

[irf4hat,irf4a,irf4b,cumirf4hat,cumirf4a,cumirf4b]=stage2irf(y,q3q');
subplot(3,3,7)
plot(time,cumirf4hat,'b-',time,cumirf4a,'b--',time,cumirf4b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -45 40])
ylabel('Real GDP')
title('(c) Oil-Specific Demand Shock')
xlabel('Quarters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IRF for CPI
y = VNdata(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0:12;

[irf2hat,irf2a,irf2b,cumirf2hat,cumirf2a,cumirf2b]=stage2irf(y,q1q');
subplot(3,3,2)
plot(time,-cumirf2hat,'b-',time,-cumirf2a,'b--',time,-cumirf2b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('CPI')
title('(d) Oil Supply Shock')

[irf3hat,irf3a,irf3b,cumirf3hat,cumirf3a,cumirf3b]=stage2irf(y,q2q');
subplot(3,3,5)
plot(time,cumirf3hat,'b-',time,cumirf3a,'b--',time,cumirf3b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('CPI')
title('(e) Global Demand Shock')

[irf4hat,irf4a,irf4b,cumirf4hat,cumirf4a,cumirf4b]=stage2irf(y,q3q');
subplot(3,3,8)
plot(time,cumirf4hat,'b-',time,cumirf4a,'b--',time,cumirf4b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('CPI')
title('(f) Oil-Specific Demand Shock')
xlabel('Quarters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IRF for IR
y = VNdata(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0:12;

[irf2hat,irf2a,irf2b,cumirf2hat,cumirf2a,cumirf2b]=stage2irf(y,q1q');
subplot(3,3,3)
plot(time,-cumirf2hat,'b-',time,-cumirf2a,'b--',time,-cumirf2b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('Nominal Interst Rate')
title('(g) Oil Supply Shock')

[irf3hat,irf3a,irf3b,cumirf3hat,cumirf3a,cumirf3b]=stage2irf(y,q2q');
subplot(3,3,6)
plot(time,cumirf3hat,'b-',time,cumirf3a,'b--',time,cumirf3b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel(['Nominal Interst Rate'])
title('(h) Global Demand Shock')

[irf4hat,irf4a,irf4b,cumirf4hat,cumirf4a,cumirf4b]=stage2irf(y,q3q');
subplot(3,3,9)
plot(time,cumirf4hat,'b-',time,cumirf4a,'b--',time,cumirf4b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('Nominal Interst Rate')
title('(i) Oil-Specific Demand Shock')
xlabel('Quarters')

%% Replicatin figure 16
clear;

trivar

Ehat=inv(chol(SIGMA)')*Uhat(1:q,:);
q1=Ehat(1,:); q1=[(q1(1,1)+q1(1,2))/2 q1];
q2=Ehat(2,:); q2=[(q2(1,1)+q2(1,2))/2 q2];
q3=Ehat(3,:); q3=[(q3(1,1)+q3(1,2))/2 q3];

% Average monthly shocks by quarter
time=2001+1/4:1/4:2015+4/4;
for i=1:length(time)
   q1q(i)=(q1(3*(i-1)+1)+q1(3*(i-1)+2)+q1(3*(i-1)+3))/3;
   q2q(i)=(q2(3*(i-1)+1)+q2(3*(i-1)+2)+q2(3*(i-1)+3))/3;
   q3q(i)=(q3(3*(i-1)+1)+q3(3*(i-1)+2)+q3(3*(i-1)+3))/3;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=(2001+1/4:1/4:2015+4/4)';

VNdata = xlsread('VNdatasub.xls'); 
% IRF for GDP
y = VNdata(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0:12;

[irf2hat,irf2a,irf2b,cumirf2hat,cumirf2a,cumirf2b]=stage2irf(y,q1q');
subplot(3,3,1)
plot(time,-cumirf2hat,'b-',time,-cumirf2a,'b--',time,-cumirf2b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -45 40])
ylabel('Real GDP')
title('(a) Oil Supply Shock')

[irf3hat,irf3a,irf3b,cumirf3hat,cumirf3a,cumirf3b]=stage2irf(y,q2q');
subplot(3,3,4)
plot(time,cumirf3hat,'b-',time,cumirf3a,'b--',time,cumirf3b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -45 40])
ylabel('Real GDP')
title('(b) Global Demand Shock')

[irf4hat,irf4a,irf4b,cumirf4hat,cumirf4a,cumirf4b]=stage2irf(y,q3q');
subplot(3,3,7)
plot(time,cumirf4hat,'b-',time,cumirf4a,'b--',time,cumirf4b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -45 40])
ylabel('Real GDP')
title('(c) Oil-Specific Demand Shock')
xlabel('Quarters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IRF for CPI
y = VNdata(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0:12;

[irf2hat,irf2a,irf2b,cumirf2hat,cumirf2a,cumirf2b]=stage2irf(y,q1q');
subplot(3,3,2)
plot(time,-cumirf2hat,'b-',time,-cumirf2a,'b--',time,-cumirf2b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('CPI')
title('(d) Oil Supply Shock')

[irf3hat,irf3a,irf3b,cumirf3hat,cumirf3a,cumirf3b]=stage2irf(y,q2q');
subplot(3,3,5)
plot(time,cumirf3hat,'b-',time,cumirf3a,'b--',time,cumirf3b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('CPI')
title('(e) Global Demand Shock')

[irf4hat,irf4a,irf4b,cumirf4hat,cumirf4a,cumirf4b]=stage2irf(y,q3q');
subplot(3,3,8)
plot(time,cumirf4hat,'b-',time,cumirf4a,'b--',time,cumirf4b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('CPI')
title('(f) Oil-Specific Demand Shock')
xlabel('Quarters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IRF for IR
y = VNdata(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0:12;

[irf2hat,irf2a,irf2b,cumirf2hat,cumirf2a,cumirf2b]=stage2irf(y,q1q');
subplot(3,3,3)
plot(time,-cumirf2hat,'b-',time,-cumirf2a,'b--',time,-cumirf2b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('Nominal Interst Rate')
title('(g) Oil Supply Shock')

[irf3hat,irf3a,irf3b,cumirf3hat,cumirf3a,cumirf3b]=stage2irf(y,q2q');
subplot(3,3,6)
plot(time,cumirf3hat,'b-',time,cumirf3a,'b--',time,cumirf3b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel(['Nominal Interst Rate'])
title('(h) Global Demand Shock')

[irf4hat,irf4a,irf4b,cumirf4hat,cumirf4a,cumirf4b]=stage2irf(y,q3q');
subplot(3,3,9)
plot(time,cumirf4hat,'b-',time,cumirf4a,'b--',time,cumirf4b,'b:',time,zeros(size(time)),'k-'); axis([0 12 -20 30])
ylabel('Nominal Interst Rate')
title('(i) Oil-Specific Demand Shock')
xlabel('Quarters')