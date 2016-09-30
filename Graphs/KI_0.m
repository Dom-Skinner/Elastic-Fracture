% Data in this file is made in quite an inefficient way. Since for a given
% L, we wish to push lambda as high as possible for a good extrapolation.
% But, this then depends on the value of L (in a way that we wish to
% determine from this data!). Therefore, its done with guesswork and by
% three different data sets.

clear

load n995x846-early
L1 = L; lam1 = l0; K2 = 3*sqrt(2*pi)*KII0;
load n995x846-mid
L1 = [L1 L]; lam1 = [lam1 l0]; K2 = [K2 3*sqrt(2*pi)*KII0];
load n995x846-late
L1 = [L1 L]; lam1 = [lam1 l0]; K2 = [K2 3*sqrt(2*pi)*KII0];



subplot(1,3,1)
hold on

plot(lam1,K2.^2,'o-')
xlabel('\lambda_0')
ylabel('K_{II}')


subplot(1,3,2)
hold on

plot(L1,lam1,'o-')
xlabel('L')
ylabel('\lambda_0')


subplot(1,3,3)
hold on

plot(L1,K2,'o-')
xlabel('L')
ylabel('K_{II}')

fid = fopen('KI-0.csv','w');
fprintf(fid,'KII0,     l0,    KII2,      L\n');
for j = 1:numel(L1)
    fprintf(fid, '%.5e,    %.5e,   %.5e,     %.5e \n', ...
        K2(j), lam1(j), K2(j)^2, L1(j));
end

fclose(fid);
