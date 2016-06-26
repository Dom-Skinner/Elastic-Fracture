%we integrate from 0 to 2^(4k), dividing each region into 2^k segments.
%thus we have 2^(5k) data points.

%number of data points
n = 100;
%small windows around crack tip and for Cauchy integral
epsilon = 1/(n^4);

%the compression at infinity
compress = 1;

%vector of h' points
inpoint = tan((0.5:1:n-0.5)*(pi/(2*(n+1)))).^2;
%vector of h points
outpoint = tan((1:n)*(pi/(2*(n+1)))).^2;

%we give a vector of 2n points for h', and interpolate h' linearly.
%get eqns for pressure at 1, 2, ..., n-1: so 2n-2 eqns
%add in 2 eqns for the conditions at infinity.

%matrix to interpolate h' for datapoints
interpol = zeros(4*n,2*n);
for i = 1:n-1
    %to find the gradients
    interpol(i,i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(i,i+1) = 1/(inpoint(i+1)-inpoint(i));
    interpol(2*(n-1)+i,n+i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(2*(n-1)+i,n+i+1) = 1/(inpoint(i+1)-inpoint(i));
    %to find the constant terms
    interpol(n-1+i,i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(n-1+i,i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    interpol(3*(n-1)+i,n+i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(3*(n-1)+i,n+i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    
    %infinite terms
    interpol(4*n-3,n) = 1;
    interpol(4*n-2,2*n) = 1;
    
    %terms at zero
    interpol(4*n-1, 1) = 1;
    interpol(4*n, n+1) = 1;
end
    

%creates a matrix to find pressure
Kmatrix = zeros(2*(n-1), 4*n);
for i = 1:n-1
    x = outpoint(i+1);
    for j = 1:n-1
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,(n-1)+j) = uK11bint - lK11bint;
        Kmatrix(i,2*(n-1)+j) = uK12cint - lK12cint;
        Kmatrix(i,3*(n-1)+j) = uK12dint - lK12dint;
        Kmatrix(n-1+i,j) = uK21aint - lK21aint;
        Kmatrix(n-1+i,(n-1)+j) = uK21bint - lK21bint;
        Kmatrix(n-1+i,2*(n-1)+j) = uK22cint - lK22cint;
        Kmatrix(n-1+i,3*(n-1)+j) = uK22dint - lK22dint;
    end
     
    %does the panel at infinity
    lx = inpoint(n);
    ux = 10^(20);
    [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
    
    Kmatrix(i,4*n-3) = uK11bint - lK11bint;
    Kmatrix(i,4*n-2) = uK12dint - lK12dint;
    Kmatrix(n-1+i,4*n-3) = uK21bint - lK21bint;
    Kmatrix(n-1+i,4*n-2) = uK22dint - lK22dint;
    
    %does panel at the crack itself
    lx = epsilon;
    ux = inpoint(1);
    [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
    Kmatrix(i,4*n-1) = uK11bint - lK11bint;
    Kmatrix(i,4*n) = uK12dint - lK12dint;
    Kmatrix(n-1+i,4*n-1) = uK21bint - lK21bint;
    Kmatrix(n-1+i,4*n) = uK22dint - lK22dint;
end

%creates the infinity conditions
ginfty = zeros(1, 2*n);
ginfty(1,n) = 1+1/sqrt(inpoint(n));
hinfty = zeros(1, 2*n);
hinfty(1,2*n-4) = -1/4;
hinfty(1,2*n-3) = 4/3;
hinfty(1,2*n-2) = -3;
hinfty(1,2*n-1) = 4;
hinfty(1,2*n) = -25/12;

%creates matrix for the conditions
eqnmatrix = zeros(2*n,2*n);
eqnmatrix(1:2*(n-1),1:2*n) = Kmatrix*interpol;
eqnmatrix(2*n-1,:) = ginfty;
eqnmatrix(2*n,:) = hinfty;

%eqnmatrix2 = eqnmatrix(1:2*n-1,1:2*n-1);

%the forcing
inhomog = zeros(2*n-1,1);
inhomog(2*n-1) = 1;
inhomog(2*n) = 0;

%soln
soln = eqnmatrix\inhomog;
