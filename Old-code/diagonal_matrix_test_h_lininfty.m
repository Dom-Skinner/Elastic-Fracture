function [x, z, soln, Pmatrix, equation_matrix, condition] = diagonal_matrix_test_h_lininfty(n, xmax)

%we only use one of the kernels, to solve for just h.
%uses just the pressure kernel

%n = 400;
%endpoint = 10;

%number of panels using 1/sqrt(x) factor
t = round(n/2);

%the value we integrate up to
infinity = 10^(20);

%the data points
%x = tan((0:n-1)*atan(xmax)/(n-1));
x(1:t) = tan((0:t-1)*(pi/4)/(t-1));
x(t+1:n) = (1:n-t)*(xmax-1)/(n-t) + ones(1,n-t);

%z = tan((0.5:1:n-1.5)*atan(xmax)/(n-1));
z(1:t-1) = tan((0.5:t-1.5)*(pi/4)/(t-1));
z(t:n-1) = (0.5:n-t-0.5)*(xmax-1)/(n-t) + ones(1,n-t);

%also we have linear panels up to panel nfin.
nfin = n-1;

%also must add special coefficients at infinity. put these at the end.
interpol = zeros(2*nfin,n);

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

%creates a matrix to find pressure
Kmatrix = zeros(n-1, 2*nfin);
for i = 1:n-1
    %does the 1/sqrt(x) panels
    for j = 1:t-1
        lx = x(j);
        ux = x(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatricesadjusted(lx,z(i));
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatricesadjusted(ux,z(i));
            
        Kmatrix(i,j) = uK12cint - lK12cint;
        Kmatrix(i,nfin+j) = uK12dint - lK12dint;
    end
    
    %does the linear panels
    for j = t:nfin
        if j == n-1
            lx = x(n-1);
            ux = infinity;
        else
            lx = x(j);
            ux = x(j+1);
        end

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,z(i));
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,z(i));
            
        Kmatrix(i,j) = uK12cint - lK12cint;
        Kmatrix(i,nfin+j) = uK12dint - lK12dint;
    end
    end

Pmatrix = Kmatrix*interpol;

%the linear coeff of h'
hinfty = zeros(1,n);
hinfty(n-1:n) = interpol(n-1,n-1:n);

equation_matrix = zeros(n,n);
equation_matrix(1:n-1,:) = Kmatrix*interpol;
equation_matrix(n,:) = hinfty;

inhomog = zeros(n,1);
inhomog(n) = 1;

condition = rcond(equation_matrix);

soln = equation_matrix\inhomog;





