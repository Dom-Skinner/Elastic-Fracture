function [v,w] = spacing(x,z,L,r)

%L*(sin(pi*0.5*(0:r-1)/r).^3 -1);
v = [L*sin(pi*0.5*(0:r-1)/r).^2 , x+L] ;
w = [L*sin(pi*0.5*(0.5:r-1.5)/r).^2 - L , 0.5*(v(r)-L), z]+L;

end