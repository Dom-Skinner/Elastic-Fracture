function [K11aint, K11bint, K12cint, K12dint, K21aint, K21bint, K22cint, K22dint] = linkernelmatricesadjusted(x,z)

%finds the kernel K11

K11aint = real(4.*(4+3.*x.^2+(-4).*x.*z+z.^2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2));

  
  
K11bint = real(8.*(4+(x+(-1).*z).^2).^(-2).*(x+(-1).*z));

%finds the kernel K12

K12cint = real(2.*(4+(x+(-1).*z).^2).^(-2).*((-8).*x+(-1).*(4+(x+(-1).*z).^2).*z) ...
  +(1/2).*z.*log(4+(x+(-1).*z).^2)+(-1).*z.*log(x+(-1).*z));
  
  
K12dint = real((-2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2).*(12+x.^2+(-2).*x.*z+z.^2)+( ...
  1/2).*log(4+(x+(-1).*z).^2)+(-1).*log(x+(-1).*z));



%finds the kernel K21

K21aint = -real(16.*x.*(4+(x+(-1).*z).^2).^(-2)+(-8).*x.*(4+(x+(-1).*z).^2).^(-1)+ ...
  2.*(4+(x+(-1).*z).^2).^(-1).*z+4.*atan((1/2).*(x+(-1).*z))+(-1/2) ...
  .*z.*log(4+(x+(-1).*z).^2)+z.*log(x+(-1).*z));
  
  
  
K21bint = -real(2.*(8+(-3).*(4+(x+(-1).*z).^2)).*(4+(x+(-1).*z).^2).^(-2)+(-1/2).* ...
  log(4+(x+(-1).*z).^2)+log(x+(-1).*z));




%finds the kernel K22

K22cint = -real(4.*(4+3.*x.^2+(-4).*x.*z+z.^2).*(4+x.^2+(-2).*x.*z+z.^2).^(-2));

  
  
K22dint = -real(8.*(4+(x+(-1).*z).^2).^(-2).*(x+(-1).*z));
  
  
  
return
end