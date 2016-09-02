clear
load K_lambda
K2 = 3*sqrt(2*pi)*hprime_data(1,:);
plot(K,K2,'o-')

hold on
plot([0,2.5],1.385*[1,1],'k')


plot([1.9322,2.5],1.5054*[1,1],'k')
plot(1.9322*[1,1],[1.5054,2.5],'k')


axis([0,2.5,1,2])


fid = fopen('catagory.csv','w');
fprintf(fid,'KI,     KII\n');
for j = 1:numel(K)
    fprintf(fid, '%.5e,    %.5e\n', K(j), K2(j));
end

fclose(fid);