clear all;

cd('C:\Users\Tung Anh\Documents\Redoing codes for diss')

global h t

% Load monthly data for 1974.1 - 2021.11, in order:
% 1. Growth rate of world oil production 
% 2. Global real activity (index based on dry cargo shipping rates)
% 3. Real price of oil 

Oildata = xlsread('Oildata.xls'); 
y=Oildata; [t,q]=size(y); 

time=(1974+1/12:1/12:2021+11/12)';  % Time line
h=15;                               % Impulse response horizon
p=24;                               % VAR lag order
[A,SIGMA,Uhat,V,X]=olsvarc(y,p);	% VAR with intercept	
SIGMA=SIGMA(1:q,1:q);