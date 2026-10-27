function f = SE( x )

ntotal=length(x);
for i=1:ntotal
    f(i)=1+0.5*x(i)-log(x(i))+log(1-exp(-x(i)))-0.5*x(i)*(1+exp(-x(i)))/(1-exp(-x(i)));
end

end

