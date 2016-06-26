%we analyse the forward transform properties of the scalar problem
%we pick a known function that transforms well.
%we only do this for h.

%we pick h'(x) = x + 1/sqrt(x)

h = @(x) x + 1./sqrt(x);

%the functions for the kernel

K12x1int = @(x,z) real(2.*(4+(x+(-1).*z).^2).^(-2).*((-8).*x+(-1).*(4+(x+(-1).*z).^2).*z) ...
  +(1/2).*z.*log(4+(x+(-1).*z).^2)+(-1).*z.*log(x+(-1).*z));

K12sqrtx0int = @(x,z) real((1/2).*(2.*x.^(1/2).*(4+(x+(-1).*z).^2).^(-2).*(4+z.^2).^(-2).*( ...
 ...
 16.*(x+(-2).*z).*(4+z.^2)+(4+x.^2+(-2).*x.*z+z.^2).*((-3).*z.*(20+ ...
...
  z.^2)+x.*(28+z.^2)))+((sqrt(-1)*2)+(-1).*z).^(-5/2).*((-15)+(sqrt( ...
...
  -1)*(-10)).*z+2.*z.^2).*atan(x.^(1/2).*((sqrt(-1)*2)+(-1).*z).^( ...
...
  -1/2))+4.*z.^(-1/2).*atanh(x.^(1/2).*z.^(-1/2))+(-1).*((sqrt(-1)* ...
...
  2)+z).^(-5/2).*((-15)+(sqrt(-1)*10).*z+2.*z.^2).*atanh(x.^(1/2).*( ...
...
  (sqrt(-1)*2)+z).^(-1/2))));

%the analytic pressure
p_an = @(z) K12x1int(10^(10),z) - K12x1int(0,z) ...
    + K12sqrtx0int(10^(10),z) - K12sqrtx0int(0,z);

colors = 'gkyrb';

%number of data points we test on
n = [100,140,200,280,400];
t = round(n/2);
endpoint = 20;

error = zeros(1,length(n));

%measures the error near the tip and the middle
tip_error = zeros(1,length(n));
mid_error = zeros(1,length(n));

tip_error_point = 0.05;
mid_error_point = 3;

figure(1);
for i = 1:length(n)
    [x,z,soln,Pmatrix, equation_matrix, conditioning] = diagonal_matrix_test_h(n(i), endpoint);

    %vector for input h
    input = zeros(n(i),1);
    input(1:t(i)) = x(1:t(i))'.^(3/2) + ones(t(i),1);
    input(t(i)+1:n(i)) = h(x(t(i)+1:n(i))');

    %the numerical pressure
    p_num = Pmatrix*input;

    plot(z,p_an(z)'-p_num,strcat(colors(i),'.-'))
    hold on;
    
    %we evaluate an error by emulating an L^2 integration.
    for j=1:n(i)-1
        error(i) = error(i) + (x(j+1)-x(j))*(p_an(z(j))-p_num(j))^2;
    end
    error(i) = sqrt(error(i));
    
    %we now find the smallest x(j) >= tip_error_point
    j = 1;
    while z(j) < tip_error_point
        j = j+1;
    end
    tip_error(i) = p_an(z(j))-p_num(j);
    
    %we now find the smallest x(j) >= mid_error_point
    j = 1;
    while z(j) < mid_error_point
        j = j+1;
    end
    mid_error(i) = p_an(z(j))-p_num(j);    
end
xlabel('x')
ylabel('p_{analytic} - p_{numeric}')
legend('n=100','n=140','n=200','n=280','n=400')
axis([0,2,-3.5*10^(-4),0.5*10^(-4)])
print -depsc '../latex/plots/scalar_problem_test/h_forward_pressure_plot.eps'
hold off;

figure(2);
plot(1./n, error, '.-')
xlabel('n^{-1}')
ylabel('L2 difference between p_{analytic} and p_{numeric}')
print -depsc '../latex/plots/scalar_problem_test/h_forward_pressure_error.eps'

for i=1:length(n)
    fprintf('%4d & %13.16e & %13.16e \\\\ \n', n(i), tip_error(i), mid_error(i))
end






%we now test for the reverse transform! we use our knowledge of the
%analytic pressure to try reverse to obtain a numerical h.

%number of data points we test on
n = [100,140,200,280,400];
t = round(n/2);
endpoint = 20;

conditioning = zeros(1,length(n));
K = zeros(1,length(n));

h_error = zeros(1,length(n));

figure(3);
for i = 1:length(n)
    
    [x,z,soln,Pmatrix, equation_matrix, conditioning(i)] = diagonal_matrix_test_h(n(i), endpoint);    
    
    % extract from the actual h what the boundary condition is
    hcoeffmatrix = [x(n(i)-3), 1, 1/x(n(i)-3), 1/x(n(i)-3)^2; x(n(i)-2), 1, 1/x(n(i)-2), 1/x(n(i)-2)^2; ...
        x(n(i)-1), 1, 1/x(n(i)-1), 1/x(n(i)-1)^2; x(n(i)), 1, 1/x(n(i)), 1/x(n(i))^2]^(-1);
    
    inhomog = zeros(n(i),1);
    %sets the pressure
    inhomog(1:n(i)-1) = p_an(z)';
    %sets the BC
    inhomog(n(i)) = hcoeffmatrix(1,:)*h(x(n(i)-3:n(i))');
    
    %finds analytic h, with funny factor adjustment
    h_an = zeros(n(i),1);
    h_an(1:t(i)) = x(1:t(i))'.^(3/2) + ones(t(i),1);
    h_an(t(i)+1:n(i)) = h(x(t(i)+1:n(i))');
    
    %finds the numerical h!
    h_num = equation_matrix\inhomog;
    
    plot(x, h_an - h_num, strcat(colors(i),'.-'));
    hold on;
    
    %the stress intensity factor
    K(i) = h_num(1);
    
    %finds the L2 error
    for j = 1:n(i)-1
        h_error(i) = h_error(i)+(x(j+1)-x(j))*(h_an(j)-h_num(j))^2;
    end
    h_error(i) = sqrt(h_error(i));
end
xlabel('x');
ylabel('h_{analytic} - h_{numeric}')
legend('n=100','n=140','n=200','n=280','n=400')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_plot.eps'
hold off;

figure(4);
plot(n, h_error,'.-')
xlabel('n')
ylabel('L2 difference between h_{analytic} and h_{numeric}')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_error.eps'

figure(5);
plot(n, log10(conditioning), '.-')
xlabel('n')
ylabel('the condition number, to log base 10')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_conditioning.eps'

figure(6);
plot(n, K, '.-')
xlabel('n')
ylabel('The numerical stress intensity')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_sif.eps'






%we redo it with more n values
%number of data points we test on
n = [280,400,560,800,1100];
t = round(n/2);
endpoint = 20;

conditioning = zeros(1,length(n));
K = zeros(1,length(n));

h_error = zeros(1,length(n));

figure(7);
for i = 1:length(n)
    
    [x,z,soln,Pmatrix, equation_matrix, conditioning(i)] = diagonal_matrix_test_h(n(i), endpoint);    
    
    % extract from the actual h what the boundary condition is
    hcoeffmatrix = vandermonde(x(n(i)-3:n(i)));
    hcoeffmatrix = hcoeffmatrix*[x(n(i)-3)^2,0,0,0;...
    0,x(n(i)-2)^2,0,0;0,0,x(n(i)-1)^2,0;...
    0,0,0,x(n(i))^2];
    
    inhomog = zeros(n(i),1);
    %sets the pressure
    inhomog(1:n(i)-1) = p_an(z)';
    %sets the BC
    inhomog(n(i)) = hcoeffmatrix(1,:)*h(x(n(i)-3:n(i))');
    
    %finds analytic h, with funny factor adjustment
    h_an = zeros(n(i),1);
    h_an(1:t(i)) = x(1:t(i))'.^(3/2) + ones(t(i),1);
    h_an(t(i)+1:n(i)) = h(x(t(i)+1:n(i))');
    
    %finds the numerical h!
    h_num = equation_matrix\inhomog;
    
    plot(x, h_an - h_num, strcat(colors(i),'.-'));
    hold on;
    
    %the stress intensity factor
    K(i) = h_num(1);
    
    %finds the L2 error
    for j = 1:n(i)-1
        h_error(i) = h_error(i)+(x(j+1)-x(j))*(h_an(j)-h_num(j))^2;
    end
    h_error(i) = sqrt(h_error(i));
end
xlabel('x');
ylabel('h_{analytic} - h_{numeric}')
legend('n=280','n=400','n=560','n=800','n=1100')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_plot_2.eps'
hold off;

figure(8);
plot(n, h_error,'.-')
xlabel('n')
ylabel('L2 difference between h_{analytic} and h_{numeric}')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_error_2.eps'

figure(9);
plot(n, log10(conditioning), '.-')
xlabel('n')
ylabel('the condition number, to log base 10')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_conditioning_2.eps'

figure(10);
plot(n, K, '.-')
xlabel('n')
ylabel('The numerical stress intensity')
print -depsc '../latex/plots/scalar_problem_test/h_forward_backward_sif_2.eps'

