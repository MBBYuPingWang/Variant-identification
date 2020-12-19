function df=smooth_fun_dev(x,A,Kin,y,lambda,p,mu,weight)
df=2*A'*Kin*(A*x-y)+lambda*sum(Si_dev_v(x.*weight,p,mu));