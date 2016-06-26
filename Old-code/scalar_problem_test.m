%we analyse the convergence in the scalar problem

%first we analyze convergence with n of inverse problem, using endpoint =
%10
n = [100,140,200,280,400];
endpoint = 10;

colors = 'gkyrb';
figure(1);

gcondition = zeros(1,length(n));
hcondition = zeros(1,length(n));

K_g = zeros(1,length(n));
K_h = zeros(1,length(n));

for i=1:length(n)
    [x, z, soln, Pmatrix, equation_matrix, gcondition(i)] = diagonal_matrix_test_g(n(i),endpoint);
    K_g(i) = soln(1);
    plot(x, soln, strcat(colors(i),'.-'))
    hold on;
end
legend('n=100','n=140','n=200','n=280','n=400')
xlabel('x');
ylabel('dg/dx');
print -depsc '../latex/plots/scalar_problem_test/g_with_n_plot.eps'
hold off;

figure(2);
plot(n, log10(gcondition),'.-');
xlabel('n')
legend('conditioning of g problem');
print -depsc '../latex/plots/scalar_problem_test/g_with_n_condition.eps'

figure(3);
for i=1:length(n)
    [x, z, soln, Pmatrix, equation_matrix, hcondition(i)] = diagonal_matrix_test_h(n(i),endpoint);
    K_h(i) = soln(1);
    plot(x, soln, strcat(colors(i),'.-'))
    hold on;
end
legend('n=100','n=140','n=200','n=280','n=400')
xlabel('x');
ylabel('dh/dx');
hold off;
print -depsc '../latex/plots/scalar_problem_test/h_with_n_plot.eps'

figure(4);
plot(n, log10(hcondition), '.-')
xlabel('n')
legend('conditioning of h problem')
print -depsc '../latex/plots/scalar_problem_test/h_with_n_condition.eps'

figure(5);
plot(1./sqrt(n), K_g, '.-')
xlabel('n^{-1/2}')
legend('SIF in g problem')
print -depsc '../latex/plots/scalar_problem_test/Kg_with_n.eps'

figure(6);
plot(1./sqrt(n), K_h, '.-')
xlabel('n^{-1/2}')
legend('SIF in h problem')
print -depsc '../latex/plots/scalar_problem_test/Kh_with_n.eps'

clear

%now we analyse the dependence on endpoint. we use n = 280 for all
%calculations.
n = 280;
endpoint = [10,14,20,28,40];

colors = 'gkyrb';

figure(7);
for i=1:length(endpoint)
    [x, z, soln, Pmatrix, equation_matrix, gcondition(i)] = diagonal_matrix_test_g(n, endpoint(i));
    K_g(i) = soln(1);
    plot(x, soln, strcat(colors(i),'.-'))
    hold on;
end
legend('x_n = 10','x_n = 14','x_n = 20','x_n = 28','x_n = 40')
xlabel('x');
ylabel('dg/dx');
hold off;
print -depsc '../latex/plots/scalar_problem_test/g_with_end_plot.eps'

figure(8);
plot(endpoint, log10(gcondition),'.-');
xlabel('x_n')
legend('conditioning of g problem');
print -depsc '../latex/plots/scalar_problem_test/g_with_end_condition.eps'

figure(9);
for i=1:length(endpoint)
    [x, z, soln, Pmatrix, equation_matrix, hcondition(i)] = diagonal_matrix_test_h(n,endpoint(i));
    K_h(i) = soln(1);
    plot(x, soln, strcat(colors(i),'.-'))
    hold on;
end
legend('x_n = 10','x_n = 14','x_n = 20','x_n = 28','x_n = 40')
xlabel('x');
ylabel('dh/dx');
print -depsc '../latex/plots/scalar_problem_test/h_with_end_plot.eps'
hold off;

figure(10);
plot(endpoint, log10(hcondition), '.-')
xlabel('x_n')
legend('conditioning of h problem')
print -depsc '../latex/plots/scalar_problem_test/h_with_end_condition.eps'
