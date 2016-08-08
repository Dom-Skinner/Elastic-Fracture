n_val = [639, 720, 830, 962, 1124, 1286, 1440 ];
x_val = [893, 906, 946, 921, 934, 944, 893 ];
s = 0.138673;
l0 = 0.591;

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

    plot(x,h0_LEFM_23,'o-',z,h0_LEFM_23_z);
    %[kernel_matrix, ~] = pressure_shear_matrix_s(x,z,s,t);
    save(file, '-append', 'h0_LEFM_23','h0_LEFM_23_z','l0');
end