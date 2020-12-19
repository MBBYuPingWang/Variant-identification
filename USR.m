function solution=USR(A,y,lambda,Kin,p,weight)
% Function for USR method

% Input
% A: independent matrix (genetic score matrix)
% y: response variable (phenotype data)
% lambda: tuning parameter that control the sparsity level
% Kin: inverse of Kinship matrix
% p: Lp norm based regularization
% weight: weight coefficient for weighted method

% Output
% solution: the sparse solution of USR

[n n_snp]=size(A);
if nargin < 4
    Kin=eye(n,n);
end
if nargin < 5
    p=0.5;
end
if nargin < 6
    weight=ones(n_snp,1);
end

esp0=1e-3;
r=1;
ro=0.8;
delta=0.3;
mu=0.1;

x0=L_half_Fw(A,y,lambda,Kin,weight);
support=find(x0~=0);
B=A(:,support);
xb=x0(support);
weight_supp=weight(support);
mu_loop=0;
x_mu=xb;
d_mu=1;

% Lower bound
L1=zeros(length(support),1);
L2=(lambda*p*sqrt(length(support))/(2*norm(B)*norm(Kin)*sqrt(y'*Kin*y)))^(1/(1-p));
for i=1:length(support)
    L1(i)=(lambda*weight_supp(i)*p*(1-p)/(2*B(:,i)'*Kin*B(:,i)))^(1/(2-p));
end

while norm(d_mu)>1e-5*n_snp
    x_mu_old=x_mu;
    g0=smooth_fun_dev(x_mu,B,Kin,y,lambda,p,mu,weight_supp);
    d0=-g0;
    k=0;
    x=x_mu;
    d=d0;
    g=g0;
    sk=1;
    while norm(sk)>1e-5*n_snp
        d_old=d;
        g_old=g;
        x_old=x;
        i=0;
        while smooth_fun(x+(ro^i)*d,B,Kin,y,lambda,p,mu,weight_supp)>smooth_fun(x,B,Kin,y,lambda,p,mu,weight_supp)+delta*(ro^i)*g'*d
            i=i+1;
        end
        x=x+(ro^i)*d;
        
        g=smooth_fun_dev(x,B,Kin,y,lambda,p,mu,weight_supp);
        yk=g-g_old;
        sk=x-x_old;
        tk=esp0*norm(g)^r+max(0,-sk'*yk/(sk'*sk));
        zk=yk-tk*sk;
        beta=g'*zk/(d_old'*zk)-2*(zk'*zk)*(g'*d_old)/(d_old'*zk)^2;
        theta=g'*d_old/(d_old'*zk);
        d=-g+beta*d_old+theta*zk;
        k=k+1;
    end
    x_mu=x;
    for i=1:length(support)
        if abs(x_mu(i))<max(L1(i),L2)
            x_mu(i)=0;
        end
    end
    d_mu=x_mu-x_mu_old;
    mu=mu/5;
    mu_loop=mu_loop+1;
end
solution=zeros(n_snp,1);
solution(support)=x_mu;