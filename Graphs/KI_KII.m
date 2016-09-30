clear
l_val = [2,4,6,8,95];

fid1 = fopen('KI-KII-1.csv','w');
fprintf(fid1,'lambda,  L,     KI,   KII,\n');

for k = 1:numel(l_val)
    file = strcat('n995x846l00', num2str(l_val(k)));
    
    load(file)
    subplot(2,2,1)
    plot(L, 3*sqrt(2*pi)*KI,'+-');
    xlabel('L');
    ylabel('K_I');
    hold on

    subplot(2,2,2)
    plot(L, KII,'o-');
    xlabel('L');
    ylabel('K_{II}');
    hold on
    
    for j = 1:numel(L)
        fprintf(fid1, '%.5e,    %.5e,    %.5e,     %.5e \n',...
            lambda, L(j), 3*sqrt(2*pi)*KI(j), 3*sqrt(2*pi)*KII(j));
    end

end
fclose(fid1);

clear
fid1 = fopen('KI-KII-2.csv','w');
fprintf(fid1,'lambda,  L,     KI,   KII,\n');
L_val=[1,4,8,16,32];
for k = 1:numel(L_val)
    if k < 4
        file = strcat('n995x846L0', num2str(L_val(k)));
    else
        file = strcat('n995x846L', num2str(L_val(k)));
    end
    
    load(file)
    subplot(2,2,3)
    plot(lambda, KI,'+-');
    xlabel('\lambda');
    ylabel('K_I');
    hold on

    subplot(2,2,4)
    plot(lambda, KII,'o-');
    xlabel('\lambda');
    ylabel('K_{II}');
    hold on
    
    for j = 1:numel(lambda)
        fprintf(fid1, '%.5e,    %.5e,    %.5e,     %.5e \n',...
            lambda(j), L, 3*sqrt(2*pi)*KI(j), 3*sqrt(2*pi)*KII(j));
    end

end
fclose(fid1);
return