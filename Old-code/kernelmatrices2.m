function [K11aint, K11bint, K12cint, K12dint, K21aint, K21bint, K22cint, K22dint] = kernelmatrices2(x,z)

%finds the kernel K11

K11aint = real((-2).*(x.^(-1)).^(-1/2).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2) ...
...
  .^(-2).*(x.^3+(-4).*x.^2.*z+x.*((-12)+z.^2)+(-2).*(x.^(-1)).^(1/2) ...
...
  .*(4+z.^2).*(4+3.*x.^2+(-1).*(x.^(-1)).^(-1/2).*z+(-4).*x.*z+z.^2) ...
...
  )+(sqrt(-1)*(-1/2)).*((sqrt(-1)*2)+(-1).*z).^(-3/2).*atan((x.^(-1) ...
...
  ).^(-1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+(sqrt(-1)*(1/2)).*(( ...
...
  sqrt(-1)*2)+z).^(-3/2).*atanh((x.^(-1)).^(-1/2).*((sqrt(-1)*2)+z) ...
...
  .^(-1/2)));
  
  
  
K11bint = real(8.*(x+(-1).*z).*(4+x.^2+(-2).*x.*z+z.^2).^(-2)+(x.^(-1)).^(1/2).*( ...
...
  12.*z.*(4+z.^2).^(-2)+8.*(x+(-1).*z).*(4+x.^2+(-2).*x.*z+z.^2).^( ...
...
  -2)+2.*(x+(-2).*z).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2).^(-1) ...
...
  )+(sqrt(-1)*(-3/2)).*((sqrt(-1)*2)+(-1).*z).^(-5/2).*atan((x.^(-1) ...
...
  ).^(-1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+(sqrt(-1)*(-3/2)).*(( ...
...
  sqrt(-1)*2)+z).^(-5/2).*atanh((x.^(-1)).^(-1/2).*((sqrt(-1)*2)+z) ...
...
  .^(-1/2)));

%finds the kernel K12

K12cint = -real((1/2).*((-32).*x.*(4+(x+(-1).*z).^2).^(-2)+(-4).*(4+(x+(-1).*z) ...
  ...
.^2).^(-1).*z+2.*(x.^(-1)).^(1/2).*((-16).*x.*(4+(x+(-1).*z).^2) ...
 ...
 .^(-2)+(-1).*(4+(x+(-1).*z).^2).^(-1).*(x+z)+z.*(4+z.^2).^(-1))+( ...
...
  -1).*((sqrt(-1)*2)+(-1).*z).^(-3/2).*((-3)+(sqrt(-1)*(-6)).*z+2.* ...
...
  z.^2).*atan((x.^(-1)).^(-1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+4.* ...
...
  z.^(1/2).*atanh((x.^(-1)).^(-1/2).*z.^(-1/2))+(-1).*((sqrt(-1)*2)+ ...
...
  z).^(-3/2).*((-3)+(sqrt(-1)*6).*z+2.*z.^2).*atanh((x.^(-1)).^( ...
 ...
 -1/2).*((sqrt(-1)*2)+z).^(-1/2))+z.*log(4+(x+(-1).*z).^2)+(-2).* ...
...
  z.*log(x+(-1).*z)));
  
  
K12dint = -real((1/2).*((-32).*(4+(x+(-1).*z).^2).^(-2)+(-4).*(4+(x+(-1).*z).^2) ...
 ...
 .^(-1)+2.*(x.^(-1)).^(1/2).*((-16).*(4+(x+(-1).*z).^2).^(-2)+(4+ ...
...
  z.^2).^(-2).*(28+z.^2)+(-1).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+ ...
...
  z.^2).^(-1).*(12+x.*z+z.^2))+((sqrt(-1)*2)+(-1).*z).^(-5/2).*(( ...
...
  -15)+(sqrt(-1)*(-10)).*z+2.*z.^2).*atan((x.^(-1)).^(-1/2).*((sqrt( ...
...
  -1)*2)+(-1).*z).^(-1/2))+4.*z.^(-1/2).*atanh((x.^(-1)).^(-1/2).* ...
...
  z.^(-1/2))+(-1).*((sqrt(-1)*2)+z).^(-5/2).*((-15)+(sqrt(-1)*10).* ...
...
  z+2.*z.^2).*atanh((x.^(-1)).^(-1/2).*((sqrt(-1)*2)+z).^(-1/2))+ ...
...
  log(4+(x+(-1).*z).^2)+(-2).*log(x+(-1).*z)));



%finds the kernel K21

K21aint = -real((1/4).*(64.*x.*(4+(x+(-1).*z).^2).^(-2)+8.*((-4).*x+z).*(4+x.^2+( ...
...
  -2).*x.*z+z.^2).^(-1)+64.*(x.^(-1)).^(1/2).*(x.*(4+(x+(-1).*z).^2) ...
...
  .^(-2)+(-1/16).*z.*(4+z.^2).^(-1)+(1/16).*((-7).*x+z).*(4+x.^2+( ...
...
  -2).*x.*z+z.^2).^(-1))+2.*((sqrt(-1)*2)+(-1).*z).^(-3/2).*((-11)+( ...
...
  sqrt(-1)*(-10)).*z+2.*z.^2).*atan((x.^(-1)).^(-1/2).*((sqrt(-1)*2) ...
...
  +(-1).*z).^(-1/2))+2.*(4+(sqrt(-1)*(-1)).*z).*atan((1/2).*(x+(-1) ...
...
  .*z))+2.*(4+sqrt(-1).*z).*atan((1/2).*(x+(-1).*z))+(-8).*z.^(1/2) ...
...
  .*atanh((x.^(-1)).^(-1/2).*z.^(-1/2))+2.*((sqrt(-1)*2)+z).^(-3/2) ...
...
  .*((-11)+(sqrt(-1)*10).*z+2.*z.^2).*atanh((x.^(-1)).^(-1/2).*(( ...
...
  sqrt(-1)*2)+z).^(-1/2))+(-1).*((sqrt(-1)*(-4))+z).*log(4+(x+(-1).* ...
...
  z).^2)+(-1).*((sqrt(-1)*4)+z).*log(4+(x+(-1).*z).^2)+4.*z.*log(x+( ...
...
  -1).*z)));
  
  
  
K21bint = -real((1/2).*(32.*(4+(x+(-1).*z).^2).^(-2)+(-12).*(4+(x+(-1).*z).^2).^( ...
...
  -1)+32.*(x.^(-1)).^(1/2).*((4+(x+(-1).*z).^2).^(-2)+(1/16).*((-20) ...
...
  +x.*z+(-7).*z.^2).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2).^(-1)+ ...
...
  (1/16).*(4+z.^2).^(-2).*(4+7.*z.^2))+(-1).*((sqrt(-1)*2)+(-1).*z) ...
...
  .^(-5/2).*((-7)+(sqrt(-1)*(-6)).*z+2.*z.^2).*atan((x.^(-1)).^( ...
 ...
 -1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+(-4).*z.^(-1/2).*atanh(( ...
 ...
 x.^(-1)).^(-1/2).*z.^(-1/2))+((sqrt(-1)*2)+z).^(-5/2).*((-7)+( ...
  ...
sqrt(-1)*6).*z+2.*z.^2).*atanh((x.^(-1)).^(-1/2).*((sqrt(-1)*2)+z) ...
...
  .^(-1/2))+(-1).*log(4+(x+(-1).*z).^2)+2.*log(x+(-1).*z)));




%finds the kernel K22

K22cint = real((-2).*(x.^(-1)).^(-1/2).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2) ...
...
  .^(-2).*(x.^3+(-4).*x.^2.*z+x.*((-12)+z.^2)+(-2).*(x.^(-1)).^(1/2) ...
...
  .*(4+z.^2).*(4+3.*x.^2+(-1).*(x.^(-1)).^(-1/2).*z+(-4).*x.*z+z.^2) ...
...
  )+(sqrt(-1)*(-1/2)).*((sqrt(-1)*2)+(-1).*z).^(-3/2).*atan((x.^(-1) ...
...
  ).^(-1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+(sqrt(-1)*(1/2)).*(( ...
...
  sqrt(-1)*2)+z).^(-3/2).*atanh((x.^(-1)).^(-1/2).*((sqrt(-1)*2)+z) ...
...
  .^(-1/2)));
  
  
  
K22dint = real(8.*(x+(-1).*z).*(4+x.^2+(-2).*x.*z+z.^2).^(-2)+(x.^(-1)).^(1/2).*( ...
...
  12.*z.*(4+z.^2).^(-2)+8.*(x+(-1).*z).*(4+x.^2+(-2).*x.*z+z.^2).^( ...
...
  -2)+2.*(x+(-2).*z).*(4+z.^2).^(-1).*(4+x.^2+(-2).*x.*z+z.^2).^(-1) ...
...
  )+(sqrt(-1)*(-3/2)).*((sqrt(-1)*2)+(-1).*z).^(-5/2).*atan((x.^(-1) ...
...
  ).^(-1/2).*((sqrt(-1)*2)+(-1).*z).^(-1/2))+(sqrt(-1)*(-3/2)).*(( ...
...
  sqrt(-1)*2)+z).^(-5/2).*atanh((x.^(-1)).^(-1/2).*((sqrt(-1)*2)+z) ...
...
  .^(-1/2)));
  
  
  
return
end