%n_val = [639, 720, 830, 962, 1124, 1286, 1440 ];
%x_val = [893, 906, 946, 921, 934, 944, 893 ];
%{  
n_val = [350,407,465,524,815];
x_val = [873,822,819,846,846];
s = 0.138673;
l0 = 0.0591;

for k = 1:numel(n_val)
    file = strcat('n', num2str(n_val(k)), 'x', num2str(x_val(k)));
    load(file)
    K = 3*sqrt(2*pi)*KI;
            
    [~ ,h0_prime_LEFM] = interpolate_hprime(x,n,hprime_data,K,0.5,l0);
    
    h0_prime_LEFM_23 = convert(0.5,2/3,n,t,x,h0_prime_LEFM);
    h_coefficient_matrix_23 = hprime_to_h_s(x,2/3,t);
    %
    h0_LEFM_23 = h_integrate(h0_prime_LEFM_23,x,n,t, ... 
    h_coefficient_matrix_23,2/3);
    h0_LEFM_23_z = h_integrate(h0_prime_LEFM_23,z,n-1,t, ...
    h_coefficient_matrix_23,2/3);

    %plot(x,h0_LEFM_23,'o-',z,h0_LEFM_23_z);
   [kernel_matrix, ~] = pressure_shear_matrix_s(x,z,s,t);
    save(file, '-append', 'h0_LEFM_23','h0_LEFM_23_z','l0','kernel_matrix');
end
%}
%{
clear 
load n605x833
n605_h0 = h0_LEFM_23;
x_605 = x;
load n639x893
n639_h0 = h0_LEFM_23;
x_639 = x;
load n698x832
plot(x,h0_LEFM_23,x_639,n639_h0,'o-',x_605,n605_h0,'o-')
%}
clear 
load n605x833
xmax = x(n);
%some values of lambda to try
lambda =0:0.002:0.0588;
scaled_K_of_c_march
s = 0.138673;
l0 = 0.0591;
K = 3*sqrt(2*pi)*KI;
u = 4-6*s;
plot(K.^u,hprime_data(1,:),'o-')
