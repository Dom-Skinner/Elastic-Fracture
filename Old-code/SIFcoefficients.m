%we find the stress intensity factors by plotting n against stuff.
%we use a fixed x_n = 30 for all our calculations.
%we use n = 100,140,200,280,400.

%choice of n
n = [100,140,200,280,400];
endpoint = 10;

KIbend = zeros(1,length(n));
KIIbend = zeros(1,length(n));
KIcomp = zeros(1,length(n));
KIIcomp = zeros(1,length(n));

%finds the compression (P) and bending(M) solutions
for i = 1:length(n)
    [KIbend(i), KIIbend(i), ~] = nofluidsolution3(n(i), endpoint,0,1,0);
    [KIcomp(i), KIIcomp(i), ~] = nofluidsolution3(n(i), endpoint,1,0,0);
end

figure(1);
plot(1./n.^2, KIbend, '.-')
xlabel('number of data points')
ylabel('KI for bending')
print -depsc '../latex/plots/KIbendconvergence3.eps'

figure(2);
plot(1./n.^2, KIIbend, '.-')
xlabel('number of data points')
ylabel('KII for bending')
print -depsc '../latex/plots/KIIbendconvergence3.eps'

figure(3);
plot(1./n.^2, KIcomp, '.-')
xlabel('number of data points')
ylabel('KI for compression')
print -depsc '../latex/plots/KIcompconvergence3.eps'

figure(4);
plot(1./n.^2, KIIcomp, '.-')
xlabel('number of data points')
ylabel('KII for compression')
print -depsc '../latex/plots/KIIcompconvergence3.eps'

for i=1:length(n)
    fprintf('%13.16f,%13.16f,%13.16f,%13.16f,%13.16f\n', ...
        n(i), KIbend(i), KIIbend(i), KIcomp(i), KIIcomp(i))
end