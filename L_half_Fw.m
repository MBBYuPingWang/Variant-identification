function solution=L_half_Fw(A,y,lambda,Kin,weight,option,k)
% Function for unified weighted L1/2 thresholding algorithm

% Input
% A: independent matrix (genetic score matrix)
% y: response variable (phenotype data)
% lambda: tuning parameter that control the sparsity level
% Kin: inverse of Kinship matrix
% weight: weight coefficient for weighted method
% option: if option='typical', solve the unified L1/2 thresholding
% algorithm; if option~='typical', solve the L1/2 thresholding under
% predetermined sparsity level
% k: predetermined sparsity level (only works when option~='typical')

% Output
% solution: the sparse solution of weighted L1/2 norm based regularization

[in n_snp]=size(A);
if nargin < 2
    error('more input arguments needed.');
end
mu=1/(norm(A))^2;
NA=(norm(A))^2;
if nargin < 3
    lambda=0.5*max(A'*y)/(in);
end
if nargin < 4
    Kin=eye(in,in);
end
if nargin < 5
    weight=ones(n_snp,1);
end
if nargin < 6
    option='typical';
end
if nargin < 7
    k=floor(in/2);
end


x=zeros(n_snp,1);
esp=1;step=0;
if strcmp(option,'typical')
%     disp('typical')
    while esp>0.001 && step<100
        x_old=x;
        B=x_old+mu*A'*Kin*(y-A*x_old);
        x=half_threshold_w(B,lambda*weight,mu);
        step=step+1;
        esp=norm(x-x_old);
        %     obj=norm(y-A*x)+lambda*sqrt(norm(x,0.5))
    end
else
    %     disp(['sparsity level: 'num2str(k)])
    while esp>0.001 && step<1000
        x_old=x;
        B=x_old+mu*A'*Kin*(y-A*x_old);
        Babs=abs(B)./weight;
        SB=sort(Babs,'descend');
        lambda=sqrt(96)/9*NA*SB(k)^1.5;
        x=half_threshold_w(B,lambda*weight,mu);
        step=step+1;
        esp=norm(x-x_old);
        %     obj=norm(y-A*x)+lambda*norm(x,0.5)^2
    end
end
solution=x;