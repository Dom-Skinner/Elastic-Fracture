function interp = find_h1_prime(n,hprime_data,K,h0_prime)
% This function assumes we accurately know h0_prime to a reasonable 
% degree of accuracy. It then interpolates via the equation
% \[ K^{-u}(h'-h_0') = C h_1' + \epsilon(K) C h_2' + \dots \]
% By varying K, one hopes to obtain a reasonable estimate for Ch_1'
% Won't be as good of an approxiamation as our estimate of h_0'
s = 0.138673;
u = 4 - 6*s;

h1_prime_data = zeros(n,numel(K));
for k = 1:n
    h1_prime_data(k,:) = (hprime_data(n+k,:)-h0_prime(k)).*K.^(-u);
end

interp = zeros(1,n);
er = zeros(1,n);
pc_er = zeros(1,n);

for l = 1:n
    p1 = polyfit(K(end-3:end-2).^u , h1_prime_data(l,end-3:end-2),1);
    p2 = polyfit(K(end-4:end-2).^u , h1_prime_data(l,end-4:end-2),2);
    interp(l) = p1(2);
    er(l)    = abs(p1(2)-p2(3));
    pc_er(l) = abs(p1(2)-p2(3))/p1(2);
end

fprintf('Approximate interpolation error = %.2e%%\n',100*max(pc_er(40:end)))

end
