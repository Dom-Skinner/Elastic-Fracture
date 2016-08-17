function [w_coeff, e_coeff, r_coeff] = get_integral_coeff(x1,x2,H1,H2,alpha,beta)

w_co1 = (x2^(1+3*beta)/H2^3 - x1^(1+3*beta)/H1^3) ; 
e_co1 = (x2^(3*beta)/H2^3 - x1^(3*beta)/H1^3);
r_co1 = (x2^(3*beta-alpha)/H2^3 - x1^(3*beta-alpha)/H1^3);

w_co2 = (x1^(3*beta)/H1^3 - x2^(3*beta)/H2^3);
e_co2 = (x1^(3*beta-1)/H1^3 - x2^(3*beta-1)/H2^3);
r_co2 = (x1^(3*beta-alpha-1)/H1^3 - x2^(3*beta-alpha-1)/H2^3);


c1 = (x2^(alpha-3*beta+2)-x1^(alpha-3*beta+2))/(alpha-3*beta+2)/(x2-x1);
c2 = x2*x1*(x2^(alpha-3*beta+1)-x1^(alpha-3*beta+1))/(alpha-3*beta+1)/(x2-x1);

w_coeff = 2*(c1*w_co1 + c2*w_co2);
e_coeff = 2*(c1*e_co1 + c2*e_co2);
r_coeff = 2*(c1*r_co1 + c2*r_co2);
    
end