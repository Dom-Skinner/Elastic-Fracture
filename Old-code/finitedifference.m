function [derivatives] = finitedifference(x,n,l)

k = length(x);

inversevand = vandermonde(x);

deriv = zeros(1,k);
for i = 1:k-n
    deriv(i) = factorial(k-i)/factorial(k-i-n)*x(l)^(k-i-n);
end

derivatives =  deriv*inversevand;

end