function [KI, KII, conditioning, soln, gprime, hprime, inpoint, outpoint, Pmatrix, interpol, gcoeffmatrix, hcoeffmatrix] ...
    = linear_spline_linear_spacing(n, endpoint, Pload,Mload)

%n = 200;
%endpoint = 30;
%Pload = 0;
%Mload = 1;
%Nload = 0;

%only measure the pressure up to a fixed point endpoint, to increase the resolution.

%at infinity uses an ax^2 + bx + c + d/x + e/x^2 expansion for h'
%at infinity uses an ax + b + c/x + d/x^2 expansion for g'
%IMPORTANT!! contrast with nofluidsolution.m

%does NOT impose conditions on the crack
%measures at n-2 points
%uses a mixed funny factor: 1/sqrt(x) for the first t panels.
%around the crack tip and at infinity, uses LINEAR not constant plot.
%fixes g'' and h''' at infinity as well.

%number of panels using 1/sqrt(x) factor
t = round(n/5);
%t = round(10*pi*(n+1)/(8*(endpoint-1)+pi));
%small windows around crack tip
epsilon = 0;
%the value we integrate up to
infinity = 10^(20);


%vector of h' points
inpoint = zeros(1,n);
inpoint(1:t-1) = tan((0:t-2)*(pi/4)/(t-1));
inpoint(t:n) = ((endpoint-1)/(n-t))*(t:n)+(t*endpoint-n)/(t-n)*ones(1,n-t+1);

%inpoint = tan((0:n-1)*atan(endpoint)/(n-1));

%vector of h points
outpoint(1:t-1) = tan((0.5:t-1.5)*(pi/4)/(t-1));
outpoint(t:n-1) = ((endpoint-1)/(n-t))*(t+0.5:n-0.5)+(t*endpoint-n)/(t-n)*ones(1,n-t);

%outpoint = tan((0.5:1:n-1.5)*atan(endpoint)/(n-1));

%we give a vector of 2n points for h', and interpolate h' linearly.
%get eqns for pressure at 1, 2, ..., n-1: so 2n-2 eqns
%add in 2 eqns for the conditions at infinity.

%matrix to interpolate h' for datapoints
%the matrix must be different at the t^th panel: we must ensure continuity
%between the two representations for the solution

%also we have linear panels up to panel nfin.
nfin = n-1;

%also must add special coefficients at infinity. put these at the end.
interpol = zeros(4*nfin+10,2*n);
%for the 1/sqrt(x) part
for i = 1:t-1
    %to find the gradients
    interpol(i,i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(i,i+1) = 1/(inpoint(i+1)-inpoint(i));
    interpol(2*nfin+i,n+i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(2*nfin+i,n+i+1) = 1/(inpoint(i+1)-inpoint(i));
    %to find the constant terms
    interpol(nfin+i,i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(nfin+i,i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    interpol(3*nfin+i,n+i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(3*nfin+i,n+i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    
end

%matching the two parts
%to find the gradients
interpol(t,t) = -1/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(t,t+1) = 1/(inpoint(t+1)-inpoint(t));
interpol(2*nfin+t,n+t) = -1/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(2*nfin+t,n+t+1) = 1/(inpoint(t+1)-inpoint(t));
%to find the constant terms
interpol(nfin+t,t) = inpoint(t+1)/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(nfin+t,t+1) = -inpoint(t)/(inpoint(t+1)-inpoint(t));
interpol(3*nfin+t,n+t) = inpoint(t+1)/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(3*nfin+t,n+t+1) = -inpoint(t)/(inpoint(t+1)-inpoint(t));

%for the purely linear part
for i = t+1:n-1
    %to find the gradients
    interpol(i,i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(i,i+1) = 1/(inpoint(i+1)-inpoint(i));
    interpol(2*nfin+i,n+i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(2*nfin+i,n+i+1) = 1/(inpoint(i+1)-inpoint(i));
    %to find the constant terms
    interpol(nfin+i,i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(nfin+i,i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    interpol(3*nfin+i,n+i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(3*nfin+i,n+i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    
end

%to get the terms at infinity

%use ax + b + c/x + d/x^2 for g'
gcoeffmatrix = vandermonde(inpoint(n-2:n));
gcoeffmatrix = gcoeffmatrix*[inpoint(n-2)^2,0,0;0,inpoint(n-1)^2,0;0,0,inpoint(n)^2];
interpol(4*nfin+3:4*nfin+5,n-2:n) = gcoeffmatrix;

%use ax^2 + bx + c + d/x + e/x^2 for h'
hcoeffmatrix = vandermonde(inpoint(n-3:n));
hcoeffmatrix = hcoeffmatrix*[inpoint(n-3)^2,0,0,0;...
    0,inpoint(n-2)^2,0,0;0,0,inpoint(n-1)^2,0;...
    0,0,0,inpoint(n)^2];
interpol(4*nfin+7:4*nfin+10,2*n-3:2*n) = hcoeffmatrix;
    


%creates a matrix to find pressure
Kmatrix = zeros(2*n-2, 4*nfin+8);
for i = 1:n-1
    x = outpoint(i);
    %does the 1/sqrt(x) panels
    for j = 1:t-1
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatricesadjusted(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatricesadjusted(ux,x);
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,nfin+j) = uK11bint - lK11bint;
        Kmatrix(i,2*nfin+j) = uK12cint - lK12cint;
        Kmatrix(i,3*nfin+j) = uK12dint - lK12dint;
        Kmatrix(n-1+i,j) = uK21aint - lK21aint;
        Kmatrix(n-1+i,nfin+j) = uK21bint - lK21bint;
        Kmatrix(n-1+i,2*nfin+j) = uK22cint - lK22cint;
        Kmatrix(n-1+i,3*nfin+j) = uK22dint - lK22dint;
    end
    
    %does the linear panels
    for j = t:nfin
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,nfin+j) = uK11bint - lK11bint;
        Kmatrix(i,2*nfin+j) = uK12cint - lK12cint;
        Kmatrix(i,3*nfin+j) = uK12dint - lK12dint;
        Kmatrix(n-1+i,j) = uK21aint - lK21aint;
        Kmatrix(n-1+i,nfin+j) = uK21bint - lK21bint;
        Kmatrix(n-1+i,2*nfin+j) = uK22cint - lK22cint;
        Kmatrix(n-1+i,3*nfin+j) = uK22dint - lK22dint;
    end
    
    %does the panel at infinity
    lx = inpoint(nfin+1);
    ux = infinity;
    [lK11aint, lK11bint, lK11eint, lK11gint, lK11iint, ...
    lK12cint, lK12dint, lK12fint, lK12hint, lK12jint, ...
    lK21aint, lK21bint, lK21eint, lK21gint, lK21iint, ...
    lK22cint, lK22dint, lK22fint, lK22hint, lK22jint] = quadinvandinvsquarekernelmatrices(lx,x);

    [uK11aint, uK11bint, uK11eint, uK11gint, uK11iint, ...
    uK12cint, uK12dint, uK12fint, uK12hint, uK12jint, ...
    uK21aint, uK21bint, uK21eint, uK21gint, uK21iint, ...
    uK22cint, uK22dint, uK22fint, uK22hint, uK22jint] = quadinvandinvsquarekernelmatrices(ux,x);

    %Kmatrix(i, 4*nfin+1) = uK11eint - lK11eint;
    Kmatrix(i, 4*nfin+2) = uK11aint - lK11aint;
    Kmatrix(i, 4*nfin+3) = uK11bint - lK11bint;
    Kmatrix(i, 4*nfin+4) = uK11gint - lK11gint;
    Kmatrix(i, 4*nfin+5) = uK11iint - lK11iint;
    
    Kmatrix(i, 4*nfin+6) = uK12fint - lK12fint;
    Kmatrix(i, 4*nfin+7) = uK12cint - lK12cint;
    Kmatrix(i, 4*nfin+8) = uK12dint - lK12dint;
    Kmatrix(i, 4*nfin+9) = uK12hint - lK12hint;
    Kmatrix(i, 4*nfin+10) = uK12jint - lK12jint;

    %Kmatrix(n-2+i, 4*nfin+1) = uK21eint - lK21eint;
    Kmatrix(n-1+i, 4*nfin+2) = uK21aint - lK21aint;
    Kmatrix(n-1+i, 4*nfin+3) = uK21bint - lK21bint;
    Kmatrix(n-1+i, 4*nfin+4) = uK21gint - lK21gint;
    Kmatrix(n-1+i, 4*nfin+5) = uK21iint - lK21iint;
    
    Kmatrix(n-1+i, 4*nfin+6) = uK22fint - lK22fint;
    Kmatrix(n-1+i, 4*nfin+7) = uK22cint - lK22cint;
    Kmatrix(n-1+i, 4*nfin+8) = uK22dint - lK22dint;
    Kmatrix(n-1+i, 4*nfin+9) = uK22hint - lK22hint;
    Kmatrix(n-1+i, 4*nfin+10) = uK22jint - lK22jint;

end

%matrix for the pressure
Pmatrix = Kmatrix*interpol;

%creates the infinity conditions

%conditions on g' = bx+c+d/x, i.e. just on b, c
ginfty = zeros(1, 2*n);

%the constant coeff
ginfty(1,n-2:n) = gcoeffmatrix(1,:);

%conditions on h' = ax^2 + bx + c + d/x, i.e. on an a,b
hinfty = zeros(1,2*n);

%the linear coeff
hinfty(1,2*n-3:2*n) = hcoeffmatrix(1,:);

%creates matrix for the conditions
eqnmatrix = zeros(2*n,2*n);
eqnmatrix(1:2*n-2,:) = Pmatrix;
eqnmatrix(2*n-1,:) = ginfty;
eqnmatrix(2*n,:) = hinfty;

%eqnmatrix2 = eqnmatrix(1:2*n-1,1:2*n-1);

%the forcing
inhomog = zeros(2*n,1);

%
inhomog(2*n-1) = (-Pload+6*Mload); %constant coeff of g'
inhomog(2*n) = (12*Mload); %x coeff of h'

%adds normalization factor
inhomog = (2*pi)^(3/2)/(8*pi)*inhomog;

%soln
conditioning = rcond(eqnmatrix);
soln = eqnmatrix\inhomog;

%modification factor
modif = ones(2*n,1);
for i = 1:t
    modif(i) = sqrt(1/inpoint(i));
    modif(n+i) = sqrt(1/inpoint(i));
end

gprime = soln(1:n).*modif(1:n);
hprime = soln(n+1:2*n).*modif(n+1:2*n);

%x1 = inpoint(1);
%x2 = inpoint(2);

%KI = (x2*soln(n+1)-soln(n+2)*x1)/(x2-x1);
%KII = (x2*soln(1)-soln(2)*x1)/(x2-x1);

KI = soln(n+1);
KII = soln(1);
