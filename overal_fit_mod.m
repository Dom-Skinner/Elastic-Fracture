load overall_fit
clearvars -except x z t n xmax
lambda = 0:0.0007:0.0588;
scaled_K_of_c_march
K = 3*sqrt(2*pi)*KI;
K2 = 3*sqrt(2*pi)*hprime_data(1,:);
clearvars hprime KI hprime_data hprime_start 
s = 0.138673;
u = 4 - 6*s;
p = polyfit(K.^u,lambda,2);

fid = fopen('overall-fit-mod.csv','w');
fprintf(fid,'KI,     KII,    V\n');
for j = 1:numel(K)
    fprintf(fid, '%.5e,    %.5e,   %.5e\n', ...
        K(j), K2(j), (36/pi)*lambda(j));
end
fclose(fid);