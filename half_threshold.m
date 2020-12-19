function y=half_threshold(x,lambda,mu)
n=length(x);
y=zeros(size(x));
t=54^(1/3)/4*(lambda*mu)^(2/3);
for i=1:n
    if abs(x(i))>t
        y(i)=f_l_mu(x(i),lambda,mu);
    end
end