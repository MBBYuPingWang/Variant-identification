function df=smooth_fun_dev_o(x,A,Kin,y,lambda,p,mu)
df=2*A'*Kin*(A*x-y)+lambda*sum(Si_dev_v(x,p,mu));