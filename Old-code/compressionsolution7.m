%we check if specifying a square root singularity / finiteness of soln at
%origin is enough to uniquely give a boundary condition.

%imposes conditions on the crack
%measures at n-2 points
%uses 1+1/sqrt(x) as funny factor

%number of data points
n = 100;
%small windows around crack tip and for Cauchy integral
epsilon = 1/(n^4);

%the compression at infinity
compress = 1;

%vector of h' points
inpoint = (0.005:0.01:0.995)./(1-(0.005:0.01:0.995));
%vector of h points
outpoint = (0.01:0.01:0.99)./(1-(0.01:0.01:0.99));

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
Kmatrix = zeros(2*(n-2), 4*n);
for i = 1:n-2
    x = outpoint(i+1);
    for j = 1:n-1
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices2(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices2(ux,x);
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,(n-1)+j) = uK11bint - lK11bint;
        Kmatrix(i,2*(n-1)+j) = uK12cint - lK12cint;
        Kmatrix(i,3*(n-1)+j) = uK12dint - lK12dint;
        Kmatrix(n-2+i,j) = uK21aint - lK21aint;
        Kmatrix(n-2+i,(n-1)+j) = uK21bint - lK21bint;
        Kmatrix(n-2+i,2*(n-1)+j) = uK22cint - lK22cint;
        Kmatrix(n-2+i,3*(n-1)+j) = uK22dint - lK22dint;
    end
    
    %does the panel at infinity
    lx = inpoint(n);
    ux = 10^(20);
    [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices2(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices2(ux,x);
    
    Kmatrix(i,4*n-3) = uK11bint - lK11bint;
    Kmatrix(i,4*n-2) = uK12dint - lK12dint;
    Kmatrix(n-2+i,4*n-3) = uK21bint - lK21bint;
    Kmatrix(n-2+i,4*n-2) = uK22dint - lK22dint;
    
    %does panel at the crack itself
    lx = epsilon;
    ux = inpoint(1);
    [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices2(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices2(ux,x);
    Kmatrix(i,4*n-1) = uK11bint - lK11bint;
    Kmatrix(i,4*n) = uK12dint - lK12dint;
    Kmatrix(n-2+i,4*n-1) = uK21bint - lK21bint;
    Kmatrix(n-2+i,4*n) = uK22dint - lK22dint;
end

%matrix for the pressure
Pmatrix = Kmatrix*interpol;

%creates the infinity conditions
ginfty = zeros(1, 2*n);
ginfty(1,n) = 1+1/sqrt(inpoint(n));
%gfirstderiv = finitedifference(inpoint(n-2:n),1,3);
%ginfty(2,n-2:n) = gfirstderiv.*(1+1./sqrt(inpoint(n-2:n)));

hinfty = zeros(1,2*n);
hfirstderiv = finitedifference(inpoint(n-1:n),1,2);
%hsecondderiv = finitedifference(inpoint(n-3:n),2,4);
hinfty(1,2*n-1:2*n) = hfirstderiv;%.*(1+1./sqrt(inpoint(n-2:n)));
%hinfty(2,2*n-3:2*n) = hsecondderiv.*(1+1./sqrt(inpoint(n-3:n)));

%creates the crack condition
k=2;
crackcondition = zeros(2,2*n);
crackcondition(1,1:k) = finitedifference(inpoint(1:k),1,1);
crackcondition(2,n+1:n+k) = finitedifference(inpoint(1:k),1,1);

%ginfty(2,n-1) = -1-1/sqrt(inpoint(n-1));
%ginfty(2,n) = 1+1/sqrt(inpoint(n));
%hinfty = zeros(2, 2*n);

%hinfty(1,2*n-1) = -1-1/sqrt(inpoint(n-1));
%hinfty(1,2*n) = 1+1/sqrt(inpoint(n));
%hinfty(2,2*n-k+1:2*n) = derivcond.*(1+1./sqrt(inpoint(n-k+1:n)));

%creates matrix for the conditions
eqnmatrix = zeros(2*n,2*n);
eqnmatrix(1:2*(n-2),:) = Pmatrix;
eqnmatrix(2*n-3,:) = ginfty;
eqnmatrix(2*n-2,:) = hinfty;
eqnmatrix(2*n-1:2*n,:) = crackcondition;

%eqnmatrix2 = eqnmatrix(1:2*n-1,1:2*n-1);

%the forcing
inhomog = zeros(2*n,1);
inhomog(2*n-3) = 1;

%soln
soln = eqnmatrix\inhomog;