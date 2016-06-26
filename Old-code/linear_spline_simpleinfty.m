function [KI, KII, x, z, soln, conditioning] ...
    = linear_spline_simpleinfty(n, xmax, Pload,Mload)

%n = 250;
%endpoint = 30;
%Pload = 1;
%Mload = 0;

%uses linear representation (with 1/sqrt(x) near tip) between x_1 and x_n.
%x_n-1 is fixed to be endpoint!!!!!
%uses for h': linear representation at in-1ity extrapolated from x_n-1,
%x_n; and for g': uses constant representation from x_n onwards.

t = round(n/2);
infinity = 10^10;

%the data points
%x = tan((0:n-1)*atan(xmax)/(n-1));
x(1:t) = tan((0:t-1)*(pi/4)/(t-1));
x(t+1:n) = (1:n-t)*(xmax-1)/(n-t) + ones(1,n-t);

%z = tan((0.5:1:n-1.5)*atan(xmax)/(n-1));
z(1:t-1) = tan((0.5:t-1.5)*(pi/4)/(t-1));
z(t:n-1) = (0.5:n-t-0.5)*(xmax-1)/(n-t) + ones(1,n-t);

%creates a matrix for 4(n-1) coefficients, plus 4 infinity coefficients
interpolate_matrix = zeros(4*(n-1)+4, 2*n);
for i = 1:t-1
    %to find the gradients
    interpolate_matrix(i,i) = -1/(x(i+1)-x(i));
    interpolate_matrix(i,i+1) = 1/(x(i+1)-x(i));
    interpolate_matrix(2*(n-1)+i,n+i) = -1/(x(i+1)-x(i));
    interpolate_matrix(2*(n-1)+i,n+i+1) = 1/(x(i+1)-x(i));
    %to find the constant terms
    interpolate_matrix(n-1+i,i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(n-1+i,i+1) = -x(i)/(x(i+1)-x(i));
    interpolate_matrix(3*(n-1)+i,n+i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(3*(n-1)+i,n+i+1) = -x(i)/(x(i+1)-x(i));
    
end

%matching the two parts
%to find the gradients
interpolate_matrix(t,t) = -1/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(t,t+1) = 1/(x(t+1)-x(t));
interpolate_matrix(2*(n-1)+t,n+t) = -1/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(2*(n-1)+t,n+t+1) = 1/(x(t+1)-x(t));
%to find the constant terms
interpolate_matrix(n-1+t,t) = x(t+1)/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(n-1+t,t+1) = -x(t)/(x(t+1)-x(t));
interpolate_matrix(3*(n-1)+t,n+t) = x(t+1)/((x(t+1)-x(t))*sqrt(x(t)));
interpolate_matrix(3*(n-1)+t,n+t+1) = -x(t)/(x(t+1)-x(t));

%for the purely linear part
for i = t+1:n-1
    %to find the gradients
    interpolate_matrix(i,i) = -1/(x(i+1)-x(i));
    interpolate_matrix(i,i+1) = 1/(x(i+1)-x(i));
    interpolate_matrix(2*(n-1)+i,n+i) = -1/(x(i+1)-x(i));
    interpolate_matrix(2*(n-1)+i,n+i+1) = 1/(x(i+1)-x(i));
    %to find the constant terms
    interpolate_matrix((n-1)+i,i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix((n-1)+i,i+1) = -x(i)/(x(i+1)-x(i));
    interpolate_matrix(3*(n-1)+i,n+i) = x(i+1)/(x(i+1)-x(i));
    interpolate_matrix(3*(n-1)+i,n+i+1) = -x(i)/(x(i+1)-x(i));
end

%linear term of g'
interpolate_matrix(4*n-3,:) = 0;
%constant term of g'
interpolate_matrix(4*n-2,n) = 1;
%linear term of h'
interpolate_matrix(4*n-1,:) = interpolate_matrix(3*(n-1),:);
%constant term of h'
interpolate_matrix(4*n, :) = interpolate_matrix(4*(n-1),:);



%creates a matrix to find pressure
Kmatrix = zeros(2*n-2, 4*(n-1)+4);
for i = 1:n-1
    %does the 1/sqrt(x) panels
    for j = 1:t-1
        lx = x(j);
        ux = x(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatricesadjusted(lx,z(i));
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatricesadjusted(ux,z(i));
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,n-1+j) = uK11bint - lK11bint;
        Kmatrix(i,2*(n-1)+j) = uK12cint - lK12cint;
        Kmatrix(i,3*(n-1)+j) = uK12dint - lK12dint;
        Kmatrix(n-1+i,j) = uK21aint - lK21aint;
        Kmatrix(n-1+i,n-1+j) = uK21bint - lK21bint;
        Kmatrix(n-1+i,2*(n-1)+j) = uK22cint - lK22cint;
        Kmatrix(n-1+i,3*(n-1)+j) = uK22dint - lK22dint;
    end
    
    %does the linear panels
    for j = t:n-1
        lx = x(j);
        ux = x(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatricesadjusted(lx,z(i));
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatricesadjusted(ux,z(i));
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,n-1+j) = uK11bint - lK11bint;
        Kmatrix(i,2*(n-1)+j) = uK12cint - lK12cint;
        Kmatrix(i,3*(n-1)+j) = uK12dint - lK12dint;
        Kmatrix(n-1+i,j) = uK21aint - lK21aint;
        Kmatrix(n-1+i,n-1+j) = uK21bint - lK21bint;
        Kmatrix(n-1+i,2*(n-1)+j) = uK22cint - lK22cint;
        Kmatrix(n-1+i,3*(n-1)+j) = uK22dint - lK22dint;
    end
    
    %does the panel at in-1ity
    lx = x(n);
    ux = infinity;

    [lK11aint, lK11bint, lK12cint, lK12dint, ...
        lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatricesadjusted(lx,z(i));
    [uK11aint, uK11bint, uK12cint, uK12dint, ...
        uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatricesadjusted(ux,z(i));

    Kmatrix(i, 4*(n-1)+1) = uK11aint - lK11aint;
    Kmatrix(i, 4*(n-1)+2) = uK11bint - lK11bint;
    Kmatrix(i, 4*(n-1)+3) = uK12cint - lK12cint;
    Kmatrix(i, 4*(n-1)+4) = uK12dint - lK12dint;
    
    Kmatrix(n-1+i, 4*(n-1)+1) = uK21aint - lK21aint;
    Kmatrix(n-1+i, 4*(n-1)+2) = uK21bint - lK21bint;
    Kmatrix(n-1+i, 4*(n-1)+3) = uK22cint - lK22cint;
    Kmatrix(n-1+i, 4*(n-1)+4) = uK22dint - lK22dint;
end

%matrix for the pressure
Pmatrix = Kmatrix*interpolate_matrix;

%creates the infinity conditions: g' and h''
g_infty_condition = interpolate_matrix(4*(n-1)+2, :);
h_infty_condition = interpolate_matrix(4*(n-1)+3,:);

%creates the equation
equation_matrix = zeros(2*n,2*n);
equation_matrix(1:2*n-2,:) = Pmatrix;
equation_matrix(2*n-1,:) = g_infty_condition;
equation_matrix(2*n,:) = h_infty_condition;

%inhomogenous conditions
inhomog = zeros(2*n,1);
inhomog(2*n-1) = (-Pload+6*Mload); %constant coeff of g'
inhomog(2*n) = (12*Mload); %x coeff of h'

%adds normalization factor
inhomog = (2*pi)^(3/2)/(8*pi)*inhomog;

%soln
conditioning = rcond(equation_matrix);
soln = equation_matrix\inhomog;

KI = soln(n+1);
KII = soln(1);

return
end