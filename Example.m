clear;clc;
n_snp=200;n=50;
causal_id=(1:10)';
beta=zeros(n_snp,1);
beta(causal_id)=1;         % define the true signal (solution) to be recovered
A=randn(n,n_snp);          % generate independent variable matrix
y=A*beta;                  % generate response variable without noise
weight=ones(n_snp,1);
weight(causal_id)=0.01;    % give true signal variants a smaller weight, so that they are easier to be selected into the model
lambda=19;         % tuning parameter of USR
lambda_E=1;        % tuning parameter of unified Elastic-net
p=0.3;             % L0.3 norm regularization problem
Kin=eye(n,n);      % assume unrelated individual data
alpha=0.5;

solution_USR=USR(A,y,lambda,Kin,p,weight);                                % solution of USR
solution_u=USR_unweighted(A,y,lambda,Kin,p);                              % solution of unweighted USR
[solution_E beta0]=Elastic_net_Fw(A,y,alpha,lambda_E,Kin);                % solution of weighted Elastic-net
[solution_Ew beta0]=Elastic_net_Fw(A,y,alpha,lambda_E,Kin,weight);        % solution of unweighted Elastic-net

sum(solution_USR~=0)
sum(solution_u~=0)
sum(solution_E~=0)
sum(solution_Ew~=0)

% This example shows how the USR method give a solution closest to the true solutuion and a proper determined
% weight can increase the detection power