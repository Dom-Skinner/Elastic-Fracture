%we only use one of the kernels, to solve for just h.
%uses just the pressure kernel

n = 400;
endpoint = 10;
Mload = 1;

%number of panels using 1/sqrt(x) factor
t = round(n/2);

%the value we integrate up to
infinity = 10^(20);

%the data points
inpoint = tan((0:n-1)*atan(endpoint)/(n-1));
outpoint = tan((0.5:1:n-1.5)*atan(endpoint)/(n-1));

%also we have linear panels up to panel nfin.
nfin = n-1;

%also must add special coefficients at infinity. put these at the end.
interpol = zeros(2*nfin+4,n);

for i = 1:t-1
    %to find the gradients
    interpol(i,i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(i,i+1) = 1/(inpoint(i+1)-inpoint(i));
    %constant terms
    interpol(nfin+i,i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(nfin+i,i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
end

%matching the two parts
%to find the gradients
interpol(t,t) = -1/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(t,t+1) = 1/(inpoint(t+1)-inpoint(t));
%to find the constant terms
interpol(nfin+t,t) = inpoint(t+1)/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(nfin+t,t+1) = -inpoint(t)/(inpoint(t+1)-inpoint(t));

%for the purely linear part
for i = t+1:n-1
    %to find the gradients
    interpol(i,i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(i,i+1) = 1/(inpoint(i+1)-inpoint(i));
    %to find the constant terms
    interpol(nfin+i,i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(nfin+i,i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
end

%use bx + c + d/x + e/x^2 for h'
hcoeffmatrix = vandermonde(inpoint(n-3:n));
hcoeffmatrix = hcoeffmatrix*[inpoint(n-3)^2,0,0,0;...
    0,inpoint(n-2)^2,0,0;0,0,inpoint(n-1)^2,0;...
    0,0,0,inpoint(n)^2];
interpol(2*nfin+1:2*nfin+4,n-3:n) = hcoeffmatrix;

%creates a matrix to find pressure
Kmatrix = zeros(n-1, 2*nfin+4);
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
            
        Kmatrix(i,j) = uK12cint - lK12cint;
        Kmatrix(i,nfin+j) = uK12dint - lK12dint;
    end
    
    %does the linear panels
    for j = t:nfin
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
            
        Kmatrix(i,j) = uK12cint - lK12cint;
        Kmatrix(i,nfin+j) = uK12dint - lK12dint;
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

    Kmatrix(i, 2*nfin+1) = uK12cint - lK12cint;
    Kmatrix(i, 2*nfin+2) = uK12dint - lK12dint;
    Kmatrix(i, 2*nfin+3) = uK12hint - lK12hint;
    Kmatrix(i, 2*nfin+4) = uK12jint - lK12jint;
end

%the linear coeff of h'
hinfty = zeros(1,n);
hinfty(n-3:n) = hcoeffmatrix(1,:);

equation_matrix = zeros(n,n);
equation_matrix(1:n-1,:) = Kmatrix*interpol;
equation_matrix(n,:) = hinfty;

inhomog = zeros(n,1);
inhomog(n) = 1;

soln = equation_matrix\inhomog;





