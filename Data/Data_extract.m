
clear
%x_val = [20,24,32,38,45,48];
x_val = [20,30,40,50,60,70];

for k = 1:numel(x_val)
    clearvars -except x_val k

    n = 500;
    t = 250;  
    xmax = x_val(k);
    x = tan((0:n-1)*atan(sqrt(xmax))/(n-1)).^2;
    z = tan((0.5:1:n-1.5)*atan(sqrt(xmax))/(n-1)).^2;
    %some values of lambda to try
    lambda = 0.056:0.0004:0.0588;
    file = strcat('n',num2str(n),'x',num2str(xmax),'-modified');
    scaled_K_of_c_march

    s = 0.138673;
    %[kernel_matrix_s, ~] = pressure_shear_matrix_s(x,z,s);

    
    %save(file,'kernel_matrix_s','s','z','-append')
    save(file)
end


%{

clear
xend_val = [20,30,40,50,60, 70];
l=1;

file = strcat('n200','x',num2str(xend_val(l)));
load(file)

u = 4 - 6*s;
K = 3*sqrt(2*pi)*KI;
p1 = polyfit(K(end-1:end).^u,lambda(end-1:end),1);
l0 = p1(2);

[h0_prime ,h0_prime_LEFM] = interpolate_hprime(x,n,hprime_data,K,0.5,l0);

h0_prime_LEFM_23 = convert(0.5,2/3,n,t,x,h0_prime_LEFM);
h_coefficient_matrix_23 = hprime_to_h_s(x,2/3);
h0_LEFM_23 = h_integrate(h0_prime_LEFM_23,x,n,t, ... 
    h_coefficient_matrix_23,2/3);
h0_LEFM_23_z = h_integrate(h0_prime_LEFM_23,z,n-1,t, ...
    h_coefficient_matrix_23,2/3);

[~,H_LEFM_23] = linear_perturbation_solve(n,t,xmax, h0_LEFM_23,...
    h0_LEFM_23_z,l0,kernel_matrix_s,s);
save(file ,'H_LEFM_23', '-append')
%}