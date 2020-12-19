function y=Si_v(x,p,mu)
y=x;
id=(x>mu);
y(id)=abs(x(id)).^p;
y(~id)=(x(~id).^2/(2*mu)+mu/2).^p;