function f=smooth_fun_o(x,A,Kin,y,lambda,p,mu)
f=(y-A*x)'*Kin*(y-A*x)+lambda*sum(Si_v(x,p,mu));