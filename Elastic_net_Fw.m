function [solution beta0]=Elastic_net_Fw(A,y,alpha,lambda,Kin,weight)
% Function for unified weighted Elastic-net regularization

% Input
% A: independent matrix (genetic score matrix)
% y: response variable (phenotype data)
% alpha: tunning parameter that control the balance between L1 norm
% regularization and L2 norm regularization
% lambda: tuning parameter that control the sparsity level
% Kin: inverse of Kinship matrix
% weight: weight coefficient for weighted method

% Output
% solution: the sparse solution of weighted Elastic-net regularization
[n n_snp]=size(A);
if nargin < 2
    error('more input arguments needed.');
end
if nargin < 3
    alpha=0.5;
end
if nargin < 4
    lambda=0.5*max(A'*y)/(n);
end
if nargin < 5
    Kin=eye(n,n);
end
if nargin < 6
    weight=ones(n_snp,1);
end

esp=1;step=1;
beta=zeros(n_snp,1);beta0=sum(y)/n;
AK=A'*Kin;
AKA=zeros(n_snp,1);
for j=1:n_snp
    AKA(j)=AK(j,:)*A(:,j)/n;
end
while(esp>0.001 && step<50)
    step=step+1;
    beta_old=beta;
    for j=1:n_snp
        r=y-A*beta-beta0;
        z=AK(j,:)*r/n+AKA(j)*beta(j);
        beta(j)=soft_threshold(z,lambda*alpha*weight(j))/(1+lambda*(1-alpha));
    end
    beta0=(sum(y)-sum(A*beta))/n;
    esp=norm(beta-beta_old);
%     obj_e=(y-A*beta-beta0)'*Kin*(y-A*beta-beta0)/2+lambda*alpha*norm(beta,1)+lambda*(1-alpha)/2*norm(beta)^2
end
solution=beta;