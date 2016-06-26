%we find the stress intensity factors by plotting n against stuff.
%we use a fixed x_n = 30 for all our calculations.
%we use n = 100,140,200,280,400.

%choice of n
n = [100,140,200,280,400];
endpoint = [10,20,30,40];

KIbend = zeros(length(endpoint),length(n));
KIIbend = zeros(length(endpoint),length(n));
KIcomp = zeros(length(endpoint),length(n));
KIIcomp = zeros(length(endpoint),length(n));

%finds the compression (P) and bending(M) solutions
for i = 1:length(n)
    for j = 1:length(endpoint)
        [KIbend(j,i), KIIbend(j,i), ~] = linear_spline_simpleinfty(n(i), endpoint(j),0,1);
        [KIcomp(j,i), KIIcomp(j,i), ~] = linear_spline_simpleinfty(n(i), endpoint(j),1,0);
    end
end

figure(1);
plot(n, KIbend, '.-')
xlabel('number of data points')
ylabel('KI for bending')
legend('x_n = 10', 'x_n = 20', 'x_n = 30', 'x_n = 40')
print -depsc '../latex/plots/KIbendconvergence_simpleinfty.eps'

figure(2);
plot(n, KIIbend, '.-')
xlabel('number of data points')
ylabel('KII for bending')
legend('x_n = 10', 'x_n = 20', 'x_n = 30', 'x_n = 40')
print -depsc '../latex/plots/KIIbendconvergence_simpleinfty.eps'

figure(3);
plot(n, KIcomp, '.-')
xlabel('number of data points')
ylabel('KI for compression')
legend('x_n = 10', 'x_n = 20', 'x_n = 30', 'x_n = 40')
print -depsc '../latex/plots/KIcompconvergence_simpleinfty.eps'

figure(4);
plot(n, KIIcomp, '.-')
xlabel('number of data points')
ylabel('KII for compression')
legend('x_n = 10', 'x_n = 20', 'x_n = 30', 'x_n = 40')
print -depsc '../latex/plots/KIIcompconvergence_simpleinfty.eps'

for j=1:length(endpoint)
for i=1:length(n)
    fprintf('%7f & %7f & %13.16f & %13.16f & %13.16f & %13.16f \\\\ \n', ...
        n(i), endpoint(j), KIbend(j,i), KIIbend(j,i), KIcomp(j,i), KIIcomp(j,i))
end
fprintf('\n')
end