function y=half_threshold_w(x,lambda,mu)
n=length(x);
y=zeros(size(x));
t=54^(1/3)/4*(mu*lambda).^(2/3);
for i=1:n
    if abs(x(i))>t(i)
        y(i)=f_l_mu(x(i),lambda(i),mu);
    end
end