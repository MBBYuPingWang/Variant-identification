function z=soft_threshold(x,y)
if x>0 && y<x
    z=x-y;
elseif x<0 && y<abs(x)
    z=x+y;
else
    z=0;    
end
