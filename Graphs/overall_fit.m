%%{
clear 
load n815x846
xmax = x(n);
%some values of lambda to try
lambda =0:0.0012:0.0588;
scaled_K_of_c_march
%}

s = 0.138673;
l0 = 0.0591;
D = -0.00826;
K = 3*sqrt(2*pi)*KI;
u = 4-6*s;
p3 = polyfit(K.^(1.5*u),lambda-l0-D*K.^u ,1);
disp(p3(1))
er = max(abs(lambda-l0-D*K.^u-p3(1)*K.^(1.5*u)));
disp(er)
plot(K.^u,lambda,K.^u,l0+D*K.^u,'*-',K.^u,l0+D*K.^u+p3(1)*K.^(1.5*u),'o-')
legend({'numerical results','linear', 'linear + K^{3u/2} term'})
xlabel('K^{u}')
ylabel('\lambda')
