%we check if specifying a square root singularity / finiteness of soln at
%origin is enough to uniquely give a boundary condition.

%for the second half of the panels uses QUADRATIC SPLINE!!!!
%has t-1 initial panels

%does NOT impose conditions on the crack
%measures at n-2 points
%uses a mixed funny factor: 1/sqrt(x) for the first t panels.
%around the crack tip and at infinity, uses LINEAR not constant plot.
%fixes g'' and h''' at infinity as well.

%number of panels using 1/sqrt(x) factor (actually t-1 panels)
t = 200;
%number of paenls using quadratic
s = 100;

%number of data points
n = t+2*s;

%small windows around crack tip
epsilon = 10^(-10);

%the compression at infinity
Pload = 1; %(2*pi)^(3/2)/(8*pi);
Mload = 0; %(2*pi)^(3/2)/(8*pi);
Nload = 0;

%vector of h' points
inpoint = tan((0.5:1:n-0.5)*(pi/(2*(n+1))));
%vector of h points
outpoint = tan((1:n-1)*(pi/(2*(n+1))));

%we give a vector of 2n points for h', and interpolate h' linearly.
%get eqns for pressure at 1, 2, ..., n-1: so 2n-2 eqns
%add in 2 eqns for the conditions at infinity.

%we use a (ax+b)/sqrt(x) for panels 1,2,...,t-1
%use LINEAR for panel t
%use QUADRATIC for panel t+1,...,

%matrix to interpolate h' for datapoints
%the matrix must be different at the t^th panel: we must ensure continuity
%between the two representations for the solution
m = 2*(t-1) + 3*s;
interpol = zeros(2*m,2*n);
%for the 1/sqrt(x) part
for i = 1:t-1
    %to find the gradients
    interpol(i,i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(i,i+1) = 1/(inpoint(i+1)-inpoint(i));
    interpol(m+i,n+i) = -1/(inpoint(i+1)-inpoint(i));
    interpol(m+i,n+i+1) = 1/(inpoint(i+1)-inpoint(i));
    %to find the constant terms
    interpol(t-1+i,i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(t-1+i,i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    interpol(m+t-1+i,n+i) = inpoint(i+1)/(inpoint(i+1)-inpoint(i));
    interpol(m+t-1+i,n+i+1) = -inpoint(i)/(inpoint(i+1)-inpoint(i));
    
end

%for the quadratic part

%does the join
joinmatrix = eye(3);
joinmatrix(1,1) = (inpoint(t))^(-1/2);
powermatrix = vandermonde(inpoint(t:t+2));
matrix = [0,1,0;0,0,1;1,0,0]*powermatrix*joinmatrix;
%g' part
interpol(2*(t-1)+1,t:t+2) = matrix(1,:);
interpol(2*(t-1)+s+1,t:t+2) = matrix(2,:);
interpol(2*(t-1)+2*s+1,t:t+2) = matrix(3,:);
%h' part
interpol(m+2*(t-1)+1,n+t:n+t+2) = matrix(1,:);
interpol(m+2*(t-1)+s+1,n+t:n+t+2) = matrix(2,:);
interpol(m+2*(t-1)+2*s+1,n+t:n+t+2) = matrix(3,:);

%for the remaining quadratic panels
for i = 1:s-1
    powermatrix = vandermonde(inpoint(t+2*i:t+2*i+2));
    matrix = [0,1,0;0,0,1;1,0,0]*powermatrix;
    
    %g' part
    interpol(2*(t-1)+i+1,t+2*i:t+2*i+2) = matrix(1,:);
    interpol(2*(t-1)+s+i+1,t+2*i:t+2*i+2) = matrix(2,:);
    interpol(2*(t-1)+2*s+i+1,t+2*i:t+2*i+2) = matrix(3,:);
    %h' part
    interpol(m+2*(t-1)+i+1,n+t+2*i:n+t+2*i+2) = matrix(1,:);
    interpol(m+2*(t-1)+s+i+1,n+t+2*i:n+t+2*i+2) = matrix(2,:);
    interpol(m+2*(t-1)+2*s+i+1,n+t+2*i:n+t+2*i+2) = matrix(3,:);
end

%creates a matrix to find pressure
Kmatrix = zeros(2*(n-2), 2*m);
for i = 1:n-2
    x = outpoint(i);
    %does the 1/sqrt(x) panels
    for j = 2:t-1
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatrices(ux,x);
            
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,(t-1)+j) = uK11bint - lK11bint;
        Kmatrix(i,m+j) = uK12cint - lK12cint;
        Kmatrix(i,m+(t-1)+j) = uK12dint - lK12dint;
        Kmatrix(n-2+i,j) = uK21aint - lK21aint;
        Kmatrix(n-2+i,(t-1)+j) = uK21bint - lK21bint;
        Kmatrix(n-2+i,m+j) = uK22cint - lK22cint;
        Kmatrix(n-2+i,m+(t-1)+j) = uK22dint - lK22dint;
    end
    
    %does the quadratic panels
    for j = 1:s-1
        lx = inpoint(t+2*(j-1));
        ux = inpoint(t+2*j);

        [lK11aint, lK11bint, lK11eint, lK12cint, lK12dint, lK12fint, lK21aint, lK21bint, lK21eint, ...
            lK22cint, lK22dint, lK22fint] = quadkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK11eint, uK12cint, uK12dint, uK12fint, uK21aint, uK21bint, uK21eint, ...
            uK22cint, uK22dint, uK22fint] = quadkernelmatrices(ux,x);
            
        Kmatrix(i,2*(t-1)+j) = uK11aint - lK11aint;
        Kmatrix(i,2*(t-1)+s+j) = uK11bint - lK11bint;
        Kmatrix(i,2*(t-1)+2*s+j) = uK11eint - lK11eint;
        
        Kmatrix(i,m+2*(t-1)+j) = uK12cint - lK12cint;
        Kmatrix(i,m+2*(t-1)+s+j) = uK12dint - lK12dint;
        Kmatrix(i,m+2*(t-1)+2*s+j) = uK12fint - lK12fint;
        
        Kmatrix(n-2+i,2*(t-1)+j) = uK21aint - lK21aint;
        Kmatrix(n-2+i,2*(t-1)+s+j) = uK21bint - lK21bint;
        Kmatrix(n-2+i,2*(t-1)+2*s+j) = uK21eint - lK21eint;
        
        Kmatrix(n-2+i,m+2*(t-1)+j) = uK22cint - lK22cint;
        Kmatrix(n-2+i,m+2*(t-1)+s+j) = uK22dint - lK22dint;
        Kmatrix(n-2+i,m+2*(t-1)+2*s+j) = uK22fint - lK22fint;
    end
    
    %does the panel at infinity
    lx = inpoint(n-2);
    ux = 10^(20);
    [lK11aint, lK11bint, lK11eint, lK12cint, lK12dint, lK12fint, lK21aint, lK21bint, lK21eint, ...
        lK22cint, lK22dint, lK22fint] = quadkernelmatrices(lx,x);
    [uK11aint, uK11bint, uK11eint, uK12cint, uK12dint, uK12fint, uK21aint, uK21bint, uK21eint, ...
        uK22cint, uK22dint, uK22fint] = quadkernelmatrices(ux,x);    
    %we use the linear interpolation from the final panel
    j = s;
        Kmatrix(i,2*(t-1)+j) = uK11aint - lK11aint;
        Kmatrix(i,2*(t-1)+s+j) = uK11bint - lK11bint;
        Kmatrix(i,2*(t-1)+2*s+j) = uK11eint - lK11eint;
        
        Kmatrix(i,m+2*(t-1)+j) = uK12cint - lK12cint;
        Kmatrix(i,m+2*(t-1)+s+j) = uK12dint - lK12dint;
        Kmatrix(i,m+2*(t-1)+2*s+j) = uK12fint - lK12fint;
        
        Kmatrix(n-2+i,2*(t-1)+j) = uK21aint - lK21aint;
        Kmatrix(n-2+i,2*(t-1)+s+j) = uK21bint - lK21bint;
        Kmatrix(n-2+i,2*(t-1)+2*s+j) = uK21eint - lK21eint;
        
        Kmatrix(n-2+i,m+2*(t-1)+j) = uK22cint - lK22cint;
        Kmatrix(n-2+i,m+2*(t-1)+s+j) = uK22dint - lK22dint;
        Kmatrix(n-2+i,m+2*(t-1)+2*s+j) = uK22fint - lK22fint;
    
    %does panel at the crack itself
    lx = epsilon;
    ux = inpoint(2);
    [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatrices(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatrices(ux,x);
    
    %we use the same linear interpolation as for the first panel
    j = 1;
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,(t-1)+j) = uK11bint - lK11bint;
        Kmatrix(i,m+j) = uK12cint - lK12cint;
        Kmatrix(i,m+(t-1)+j) = uK12dint - lK12dint;
        Kmatrix(n-2+i,j) = uK21aint - lK21aint;
        Kmatrix(n-2+i,(t-1)+j) = uK21bint - lK21bint;
        Kmatrix(n-2+i,m+j) = uK22cint - lK22cint;
        Kmatrix(n-2+i,m+(t-1)+j) = uK22dint - lK22dint;
end

%matrix for the pressure
Pmatrix = Kmatrix*interpol;

%creates the infinity conditions
ginfty = zeros(2, 2*n);
ginfty(1,n) = 1;
gfirstderiv = finitedifference(inpoint(n-2:n),1,3);
ginfty(2,n-2:n) = gfirstderiv;

hinfty = zeros(2,2*n);
hfirstderiv = finitedifference(inpoint(n-2:n),1,3);
hsecondderiv = finitedifference(inpoint(n-2:n),2,3);
hinfty(1,2*n-2:2*n) = hfirstderiv;
hinfty(2,2*n-2:2*n) = hsecondderiv;


%ginfty(2,n-1) = -1-1/sqrt(inpoint(n-1));
%ginfty(2,n) = 1+1/sqrt(inpoint(n));
%hinfty = zeros(2, 2*n);
%hinfty(1,2*n-1) = -1-1/sqrt(inpoint(n-1));
%hinfty(1,2*n) = 1+1/sqrt(inpoint(n));
%hinfty(2,2*n-k+1:2*n) = derivcond.*(1+1./sqrt(inpoint(n-k+1:n)));

%creates matrix for the conditions
eqnmatrix = zeros(2*n,2*n);
eqnmatrix(1:2*(n-2),:) = Pmatrix;
eqnmatrix(2*n-3:2*n-2,:) = ginfty;
eqnmatrix(2*n-1:2*n,:) = hinfty;

%eqnmatrix2 = eqnmatrix(1:2*n-1,1:2*n-1);

%the forcing
inhomog = zeros(2*n,1);
inhomog(2*n-3) = (2*pi)^(3/2)/(8*pi)*(Pload+6*Mload-6*Nload*inpoint(n)); %g'
inhomog(2*n-2) = (2*pi)^(3/2)/(8*pi)*(-6*Nload); %g''
inhomog(2*n-1) = (2*pi)^(3/2)/(8*pi)*(12*Mload+12*Nload*inpoint(n)); %h''
inhomog(2*n) = (2*pi)^(3/2)/(8*pi)*(12*Nload); %h'''

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