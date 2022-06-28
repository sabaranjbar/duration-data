%% Reading the data
clc; clear;
data1 = readmatrix ('C:\Users\sabar\Documents\MATLAB\econ683.xls');
data2 = readmatrix ('C:\Users\sabar\Documents\MATLAB\X.xlsx');
%% Assigning the data to X
S=zeros(1661,13);
medu=zeros(1661,1);
famincome=zeros(1661,1);
neutral=zeros(1661,1);
Rural=zeros(1661,1);
nofsib=zeros(1661,1);
for i=1:1661
    S(i,:)=data1(i,29:41);
    medu(i,1)=data2(i+1,1);
    famincome(i,1)=data2(i+1,2);
    neutral(i,1)=data2(i+1,3);
    Rural(i,1)=data2(i+1,4);
    nofsib(i,1)=data2(i+1,6);
    e=ones(i,1);
end
X=[e,medu,famincome,neutral,Rural,nofsib];
%% Finding the highest grade achieved by each individual
K=zeros(1661,1);
for i=1:1661
  [~,b]=bounds(S(i,:));
  K(i,1)=b;
end
%% Using fminunc() to calculate the maximum likelihood
param_init=[3 0 0 0 0 0 -2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
options = optimoptions('fminunc','display','iter','MaxIterations',10000,'MaxFunEvals',10000);
[c,fval,E,output,grad,hessian]=fminunc(@(param) objfun(param,K,X),param_init,options);

sigma2 = inv(hessian)/1661;
std_beta=diag(sqrt(sigma2));

t=c./std_beta;