function [x, z, soln, Pmatrix, equation_matrix, condition] = diagonal_matrix_test_g(n, endpoint)

%we only use one of the kernels, to solve for just g.
%uses just the pressure kernel

%n = 400;
%endpoint = 10;
%Mload = 1;

%number of panels using 1/sqrt(x) factor
t = round(n/2);

%the value we integrate up to
infinity = 10^(20);

%the data points
x = tan((0:n-1)*atan(endpoint)/(n-1));
z = tan((0.5:1:n-1.5)*atan(endpoint)/(n-1));

%also we have linear panels up to panel nfin.
nfin = n-1;

%also must add special coefficients at infinity. put these at the end.
interpol = zeros(2*nfin+3,n);

for i = 1:t-1
    %to find the gradients
    interpol(i,i) = -1/(x(i+1)-x(i));
    interpol(i,i+1) = 1/(x(i+1)-x(i));
    %constant terms
    interpol(nfin+i,i) = x(i+1)/(x(i+1)-x(i));
    interpol(nfin+i,i+1) = -x(i)/(x(i+1)-x(i));
end

%matching the two parts
%to find the gradients
interpol(t,t) = -1/((x(t+1)-x(t))*sqrt(x(t)));
interpol(t,t+1) = 1/(x(t+1)-x(t));
%to find the constant terms
interpol(nfin+t,t) = x(t+1)/((x(t+1)-x(t))*sqrt(x(t)));
interpol(nfin+t,t+1) = -x(t)/(x(t+1)-x(t));

%for the purely linear part
for i = t+1:n-1
    %to find the gradients
    interpol(i,i) = -1/(x(i+1)-x(i));
    interpol(i,i+1) = 1/(x(i+1)-x(i));
    %to find the constant terms
    interpol(nfin+i,i) = x(i+1)/(x(i+1)-x(i));
    interpol(nfin+i,i+1) = -x(i)/(x(i+1)-x(i));
end

%use bx + c + d/x + e/x^2 for h'
gcoeffmatrix = vandermonde(x(n-2:n));
gcoeffmatrix = gcoeffmatrix*[x(n-2)^2,0,0;...
    0,x(n-1)^2,0;0,0,x(n)^2];
interpol(2*nfin+1:2*nfin+3,n-2:n) = gcoeffmatrix;

%creates a matrix to find pressure
Kmatrix = zeros(n-1, 2*nfin+3);
for i = 1:n-1
    %does the 1/sqrt(x) panels
    for j = 1:t-1
        lx = x(j);
        ux = x(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatricesadjusted(lx,z(i));
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatricesadjusted(ux,z(i));
            
        Kmatrix(i,j) = uK21aint - lK21aint;
        Kmatrix(i,nfin+j) = uK21bint - lK21bint;
    end
    
    %does the linear panels
    for j = t:nfin
        lx = x(j);
        ux = x(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,z(i));
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,z(i));
            
        Kmatrix(i,j) = uK21aint - lK21aint;
        Kmatrix(i,nfin+j) = uK21bint - lK21bint;
    end
    
        %does the panel at infinity
    lx = x(nfin+1);
    ux = infinity;
    [lK11aint, lK11bint, lK11eint, lK11gint, lK11iint, ...
    lK12cint, lK12dint, lK12fint, lK12hint, lK12jint, ...
    lK21aint, lK21bint, lK21eint, lK21gint, lK21iint, ...
    lK22cint, lK22dint, lK22fint, lK22hint, lK22jint] = quadinvandinvsquarekernelmatrices(lx,z(i));

    [uK11aint, uK11bint, uK11eint, uK11gint, uK11iint, ...
    uK12cint, uK12dint, uK12fint, uK12hint, uK12jint, ...
    uK21aint, uK21bint, uK21eint, uK21gint, uK21iint, ...
    uK22cint, uK22dint, uK22fint, uK22hint, uK22jint] = quadinvandinvsquarekernelmatrices(ux,z(i));

    Kmatrix(i, 2*nfin+1) = uK21bint - lK21bint;
    Kmatrix(i, 2*nfin+2) = uK21gint - lK21gint;
    Kmatrix(i, 2*nfin+3) = uK21iint - lK21iint;
end

Pmatrix = Kmatrix*interpol;

%the linear coeff of h'
ginfty = zeros(1,n);
ginfty(n-2:n) = gcoeffmatrix(1,:);

equation_matrix = zeros(n,n);
equation_matrix(1:n-1,:) = Kmatrix*interpol;
equation_matrix(n,:) = ginfty;

inhomog = zeros(n,1);
inhomog(n) = 1;

condition = rcond(equation_matrix);
soln = equation_matrix\inhomog;





