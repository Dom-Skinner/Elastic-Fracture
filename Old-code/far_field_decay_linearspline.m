%we look at how fast the solution sets into its far field behaviour
%for different values of n
%looks at just the M solution

n = [100,140,200,280,400];
endpoint = 50;

xpoints = zeros(length(n),max(n));
M_g_soln = zeros(length(n),max(n));
M_h_soln = zeros(length(n),max(n));
M_g_decay = zeros(length(n),max(n));
M_h_decay = zeros(length(n),max(n));

P_g_soln = zeros(length(n),max(n));
P_h_soln = zeros(length(n),max(n));
P_g_decay = zeros(length(n),max(n));
P_h_decay = zeros(length(n),max(n));

matrix_cond = zeros(1,length(n));

for i=1:length(n)
    [KI, KII, soln, gprime, hprime, inpoint, ...
        interpol, gcoeffmatrix, hcoeffmatrix, conditioning] ...
        = nofluidsolution3(n(i), endpoint, 0,1,0);
    g_ffield_coeffs = gcoeffmatrix*soln(n(i)-3:n(i));
    h_ffield_coeffs = hcoeffmatrix*soln(2*n(i)-4:2*n(i));
    
    xpoints(i,1:n(i)) = inpoint;
    M_g_soln(i,1:n(i)) = gprime;
    M_h_soln(i,1:n(i)) = hprime;
    
    M_g_decay(i,1:n(i)) = transpose(gprime) - g_ffield_coeffs(1)*inpoint ...
        - g_ffield_coeffs(2)*ones(1,n(i)) - g_ffield_coeffs(3)*inpoint.^(-1) ...
        - g_ffield_coeffs(4)*inpoint.^(-2);
    M_h_decay(i,1:n(i)) = transpose(hprime) - h_ffield_coeffs(1)*inpoint.^2 ...
        - h_ffield_coeffs(2).*inpoint - h_ffield_coeffs(3)*ones(1,n(i)) ...
        - h_ffield_coeffs(4)*inpoint.^(-1) - h_ffield_coeffs(5)*inpoint.^(-2);
    
    [KI, KII, soln, gprime, hprime, inpoint, ...
        interpol, gcoeffmatrix, hcoeffmatrix, conditioning] ...
        = nofluidsolution3(n(i), endpoint, 1,0,0);
    g_ffield_coeffs = gcoeffmatrix*soln(n(i)-3:n(i));
    h_ffield_coeffs = hcoeffmatrix*soln(2*n(i)-4:2*n(i));
    
    xpoints(i,1:n(i)) = inpoint;
    P_g_soln(i,1:n(i)) = gprime;
    P_h_soln(i,1:n(i)) = hprime;
    
    P_g_decay(i,1:n(i)) = transpose(gprime) - g_ffield_coeffs(1)*inpoint ...
        - g_ffield_coeffs(2)*ones(1,n(i)) - g_ffield_coeffs(3)*inpoint.^(-1) ...
        - g_ffield_coeffs(4)*inpoint.^(-2);
    P_h_decay(i,1:n(i)) = transpose(hprime) - h_ffield_coeffs(1)*inpoint.^2 ...
        - h_ffield_coeffs(2).*inpoint - h_ffield_coeffs(3)*ones(1,n(i)) ...
        - h_ffield_coeffs(4)*inpoint.^(-1) - h_ffield_coeffs(5)*inpoint.^(-2);
    
    matrix_cond(i) = conditioning;
end

colours = ['y';'g';'k';'b';'r'];

%the decay for horizontal displacement for M soln
figure(1);
hold on;
for i = 1:length(n)
    plot(xpoints(i,1:n(i)-5), log(abs(M_g_decay(i,1:n(i)-5))), strcat(colours(i),'.-'));
end
legend('n=100','n=140','n=200','n=280','n=400');
xlabel('x')
ylabel('log_{10}|far field representation - numerical soln|')
hold off;
print -depsc '../latex/plots/ffield_decay_linearspline/M_g_decay.eps'

%decay for vertical displacement for M soln
figure(2);
hold on;
for i = 1:length(n)
    plot(xpoints(i,1:n(i)-5), log(abs(M_h_decay(i,1:n(i)-5))), strcat(colours(i),'.-'));
end
legend('n=100','n=140','n=200','n=280','n=400');
xlabel('x')
ylabel('log_{10}|far field representation - numerical soln|')
hold off;
print -depsc '../latex/plots/ffield_decay_linearspline/M_h_decay.eps'

%the decay for horizontal displacement for P soln
figure(3);
hold on;
for i = 1:length(n)
    plot(xpoints(i,1:n(i)-5), log(abs(P_g_decay(i,1:n(i)-5))), strcat(colours(i),'.-'));
end
legend('n=100','n=140','n=200','n=280','n=400');
xlabel('x')
ylabel('log_{10}|far field representation - numerical soln|')
hold off;
print -depsc '../latex/plots/ffield_decay_linearspline/P_g_decay.eps'

%decay for vertical displacement for P soln
figure(4);
hold on;
for i = 1:length(n)
    plot(xpoints(i,1:n(i)-5), log(abs(P_h_decay(i,1:n(i)-5))), strcat(colours(i),'.-'));
end
legend('n=100','n=140','n=200','n=280','n=400');
xlabel('x')
ylabel('log_{10}|far field representation - numerical soln|')
hold off;
print -depsc '../latex/plots/ffield_decay_linearspline/P_h_decay.eps'

%the conditioning!
figure(5);
plot(n, log10(matrix_cond), '.-');
xlabel('number of data points');
ylabel('log_{10} (condition number of matrix)');
print -depsc '../latex/plots/ffield_decay_linearspline/conditioning.eps'
    