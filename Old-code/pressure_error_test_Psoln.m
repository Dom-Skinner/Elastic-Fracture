%this tests for the error in the pressure and shear, both at the tip and at
%the far field. We just analyse the M solution here.

n = [100, 140, 200, 280, 400];
endpoint = 100;

measurepoint = zeros(length(n),10*max(n));
pressure = zeros(length(n),10*max(n));
shear = zeros(length(n),10*max(n));

for i = 1:length(n)
    measurepoint(i,1:10*n(i)) = tan((0.05:0.1:n(i)-0.05)*(atan(endpoint))/n(i));
    
    [KI, KII, soln, gprime, hprime, inpoint, interpol] = nofluidsolution3(n(i), endpoint, 1,0,0);
    
    for j=1:10*n(i)
        x = measurepoint(i,j);
        [pressurerow, shearrow] = measurepressure(x, inpoint);
        pressure(i,j) = pressurerow*interpol*soln;
        shear(i,j) = shearrow*interpol*soln;
    end
end

%plot of comparative errors

colours = ['b'; 'r'; 'k'; 'y'; 'g'];

%picture of the pressure error near tip
figure(1);
for i = 1:length(n)
    plot(measurepoint(i,1:10*n(i)),pressure(i,1:10*n(i)), colours(i))
    hold on;
end
legend('n=100','n=140','n=200','n=280','n=400');
axis([0, 0.25, -0.2*10^(-3), 0.2*10^(-3)]);
hold off;

print -depsc '../latex/plots/pressure_error_test/Psoln_tip_picture_pressure.eps'

%plot of the shear error near tip
figure(2);
for i = 1:length(n)
    plot(measurepoint(i,1:10*n(i)),shear(i,1:10*n(i)), colours(i))
    hold on;
end
legend('n=100','n=140','n=200','n=280','n=400');
axis([0, 0.25, -0.2*10^(-3), 0.2*10^(-3)]);
hold off;

print -depsc '../latex/plots/pressure_error_test/Psoln_tip_picture_shear.eps'

%picture of the pressure error naround far field
figure(3);
for i = 1:length(n)
    plot(measurepoint(i,5*n(i)+10:10*n(i)),pressure(i,5*n(i)+10:10*n(i)), colours(i))
    hold on;
end
legend('n=100','n=140','n=200','n=280','n=400');
hold off;

print -depsc '../latex/plots/pressure_error_test/Psoln_ffield_picture_pressure.eps'

%picture of the shear error naround far field
figure(4);
for i = 1:length(n)
    plot(measurepoint(i,5*n(i)+10:10*n(i)),shear(i,5*n(i)+10:10*n(i)), colours(i))
    hold on;
end
legend('n=100','n=140','n=200','n=280','n=400');
hold off;

print -depsc '../latex/plots/pressure_error_test/Psoln_ffield_picture_shear.eps'

%plot of the pressure error in the tip solution
figure(5);
tip_pressure_error = zeros(1,length(n));
for i = 1:length(n)
    tip_pressure_error(i) = max(pressure(i,1:5*n(i)));
end
plot(n, tip_pressure_error, '.-');
legend('error in the pressure at the tip')

print -depsc '../latex/plots/pressure_error_test/Psoln_tip_max_pressure.eps'

%plot of the shear error at the tip
figure(6);
tip_shear_error = zeros(1,length(n));
for i = 1:length(n)
    tip_shear_error(i) = max(shear(i,1:5*n(i)));
end
plot(n, tip_shear_error, '.-');
legend('error in the shear at the tip')

print -depsc '../latex/plots/pressure_error_test/Psoln_tip_max_shear.eps'

%plot of the pressure error in far field
figure(7);
ffield_pressure_error = zeros(1,length(n));
for i = 1:length(n)
    ffield_pressure_error(i) = max(pressure(i,5*n(i)+10:10*n(i)));
end
plot(n, ffield_pressure_error, '.-');
legend('error in the pressure in far field')

print -depsc '../latex/plots/pressure_error_test/Psoln_ffield_max_pressure.eps'

%plot of the shear error in far field
figure(8);
ffield_shear_error = zeros(1,length(n));
for i = 1:length(n)
    ffield_shear_error(i) = max(shear(i,5*n(i)+10:10*n(i)));
end
plot(n, ffield_shear_error, '.-');
legend('error in the shear in far field')

print -depsc '../latex/plots/pressure_error_test/Psoln_ffield_max_shear.eps'

