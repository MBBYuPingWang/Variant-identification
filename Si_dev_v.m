function y=Si_dev_v(x,p,mu)
y=x;
id=(x>mu);
y(id)=p*sign(x(id)).*abs(x(id)).^(p-1);
y(~id)=p*(x(~id).^2/(2*mu)+mu/2).^(p-1).*x(~id)/mu;