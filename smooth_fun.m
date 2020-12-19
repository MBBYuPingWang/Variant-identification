function f=smooth_fun(x,A,Kin,y,lambda,p,mu,weight)
f=(y-A*x)'*Kin*(y-A*x)+lambda*sum(Si_v(x.*weight,p,mu));