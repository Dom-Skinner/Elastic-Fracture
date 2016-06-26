%we check if specifying a square root singularity / finiteness of soln at
%origin is enough to uniquely give a boundary condition.

%imposes conditions on the crack
%measures at n-2 points
%uses a mixed funny factor: 1/sqrt(x) for the first t panels.
%around the crack tip and at infinity, uses LINEAR not constant plot.

%number of data points
n = 100;
%number of panels using 1/sqrt(x) factor
t = round(n/2);
%small windows around crack tip
epsilon = 10^(-10);

%the compression at infinity
compress = (2*pi)^(3/2)/(8*pi);

%vector of h' points
inpoint = (0.5:1:n-0.5)./(n+1-(0.5:1:n-0.5));
%vector of h points
outpoint = (1:n)./(n+1-(1:n));

%we give a vector of 2n points for h', and interpolate h' linearly.
%get eqns for pressure at 1, 2, ..., n-1: so 2n-2 eqns
%add in 2 eqns for the conditions at infinity.

%matrix to interpolate h' for datapoints
%the matrix must be different at the t^th panel: we must ensure continuity
%between the two representations for the solution
interpol = zeros(4*(n-1),2*n);
%for the 1/sqrt(x) part
for i = 1:t-1
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
    
end

%matching the two parts
%to find the gradients
interpol(t,t) = -1/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(t,t+1) = 1/(inpoint(t+1)-inpoint(t));
interpol(2*(n-1)+t,n+t) = -1/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(2*(n-1)+t,n+t+1) = 1/(inpoint(t+1)-inpoint(t));
%to find the constant terms
interpol(n-1+t,t) = inpoint(t+1)/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(n-1+t,t+1) = -inpoint(t)/(inpoint(t+1)-inpoint(t));
interpol(3*(n-1)+t,n+t) = inpoint(t+1)/((inpoint(t+1)-inpoint(t))*sqrt(inpoint(t)));
interpol(3*(n-1)+t,n+t+1) = -inpoint(t)/(inpoint(t+1)-inpoint(t));

%for the purely linear part
for i = t+1:n-1
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
    
end

%infinite terms
%interpol(4*n-3,n) = 1;
%interpol(4*n-2,2*n) = 1;
    
%terms at zero
%interpol(4*n-1, 1) = 1;
%interpol(4*n, n+1) = 1;
    


%creates a matrix to find pressure
Kmatrix = zeros(2*(n-2), 4*(n-1));
for i = 1:n-2
    x = outpoint(i+1);
    %does the 1/sqrt(x) panels
    for j = 2:t-1
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatrices(ux,x);
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,(n-1)+j) = uK11bint - lK11bint;
        Kmatrix(i,2*(n-1)+j) = uK12cint - lK12cint;
        Kmatrix(i,3*(n-1)+j) = uK12dint - lK12dint;
        Kmatrix(n-2+i,j) = uK21aint - lK21aint;
        Kmatrix(n-2+i,(n-1)+j) = uK21bint - lK21bint;
        Kmatrix(n-2+i,2*(n-1)+j) = uK22cint - lK22cint;
        Kmatrix(n-2+i,3*(n-1)+j) = uK22dint - lK22dint;
    end
    
    %does the linear panels
    for j = t:n-2
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
            
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
    lx = inpoint(n-1);
    ux = 10^(20);
    [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
    
    %we use the linear interpolation from the final panel
    j = n-1;
    Kmatrix(i,j) = uK11aint - lK11aint;
    Kmatrix(i,(n-1)+j) = uK11bint - lK11bint;
    Kmatrix(i,2*(n-1)+j) = uK12cint - lK12cint;
    Kmatrix(i,3*(n-1)+j) = uK12dint - lK12dint;
    Kmatrix(n-2+i,j) = uK21aint - lK21aint;
    Kmatrix(n-2+i,(n-1)+j) = uK21bint - lK21bint;
    Kmatrix(n-2+i,2*(n-1)+j) = uK22cint - lK22cint;
    Kmatrix(n-2+i,3*(n-1)+j) = uK22dint - lK22dint;
    
    %does panel at the crack itself
    lx = epsilon;
    ux = inpoint(2);
    [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatrices(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatrices(ux,x);
    
    %we use the same linear interpolation as for the first panel
    j = 1;
    Kmatrix(i,j) = uK11aint - lK11aint;
    Kmatrix(i,(n-1)+j) = uK11bint - lK11bint;
    Kmatrix(i,2*(n-1)+j) = uK12cint - lK12cint;
    Kmatrix(i,3*(n-1)+j) = uK12dint - lK12dint;
    Kmatrix(n-2+i,j) = uK21aint - lK21aint;
    Kmatrix(n-2+i,(n-1)+j) = uK21bint - lK21bint;
    Kmatrix(n-2+i,2*(n-1)+j) = uK22cint - lK22cint;
    Kmatrix(n-2+i,3*(n-1)+j) = uK22dint - lK22dint;
end

%matrix for the pressure
Pmatrix = Kmatrix*interpol;

%creates the infinity conditions
ginfty = zeros(1, 2*n);
ginfty(1,n) = 1;
%gfirstderiv = finitedifference(inpoint(n-2:n),1,3);
%ginfty(2,n-2:n) = gfirstderiv;

hinfty = zeros(1,2*n);
hfirstderiv = finitedifference(inpoint(n-2:n),1,3);
%hsecondderiv = finitedifference(inpoint(n-4:n),2,5);
hinfty(1,2*n-2:2*n) = hfirstderiv;
%hinfty(1,2*n-4:2*n) = hsecondderiv;

%creates the crack condition
k=3;
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

Pmatrix2 = zeros(2*n-2,2*n);
Pmatrix2(1:2*(n-2),:) = Pmatrix;
Pmatrix2(2*n-3:2*n-2,:) = crackcondition;

%the forcing
inhomog = zeros(2*n,1);
inhomog(2*n-3) = compress;

%soln
soln = eqnmatrix\inhomog;

%modification factor
modif = ones(2*n,1);
for i = 1:t
    modif(i) = sqrt(1/inpoint(i));
    modif(n+i) = sqrt(1/inpoint(i));
end

gprime = soln(1:n).*modif(1:n);
hprime = soln(n+1:2*n).*modif(n+1:2*n);