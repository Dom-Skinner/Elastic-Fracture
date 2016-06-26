%we integrate from 0 to 2^(4k), dividing each region into 2^k segments.
%thus we have 2^(5k) data points.

k = 8;

%number of data points
n = 2^(k);
%small windows around crack tip and for Cauchy integral
epsilon = 2^(-2*k);

%the compression at infinity
compress = 1;

%constituent parts for the full matrix
Kmatrix = zeros(2*(n-1),4*n);
Ctymatrix = zeros(2*(n-1),4*n);
ginfty = zeros(2,4*n);
hinfty = zeros(2,4*n);

%vector of h' points
inpoint = (0.5:1:n-0.5);
%vector of h points
outpoint = (1:n);

%creates matrix of continuity

for i = 1:n-1
    x = outpoint(i);
    Ctymatrix(i,i) = x;
    Ctymatrix(i,i+1) = -x;
    Ctymatrix(i,n+i) = 1;
    Ctymatrix(i,n+i+1) = -1;
    Ctymatrix(n-1+i,2*n+i) = x;
    Ctymatrix(n-1+i,2*n+i+1) = -x;
    Ctymatrix(n-1+i,3*n+i) = 1;
    Ctymatrix(n-1+i,3*n+i+1) = -1;
end

%creates matrix for the pressure condition

for i = 1:n-1
    x = outpoint(i);
    for j = 1:n
        if j == 1
            if i == 1
                [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices(x,outpoint(1)-epsilon);
                [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices(x,epsilon);
            else
                [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices(x,outpoint(1));
                [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices(x,epsilon);
            end
        elseif j == i
            [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices(x,outpoint(j)-epsilon);
            [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices(x,outpoint(j-1));
        elseif j == i+1
            [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices(x,outpoint(j));
            [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices(x,outpoint(j-1)-epsilon);
        else
            [uK11aint, uK11bint, uK12cint, uK12dint, uK21aint, uK21bint, uK22cint, uK22dint] = kernelmatrices(x,outpoint(j));
            [lK11aint, lK11bint, lK12cint, lK12dint, lK21aint, lK21bint, lK22cint, lK22dint] = kernelmatrices(x,outpoint(j-1));
        end
        
        Kmatrix(i,j) = uK11aint - lK11aint;
        Kmatrix(i,n+j) = uK11bint - lK11bint;
        Kmatrix(i,2*n+j) = uK12cint - lK12cint;
        Kmatrix(i,3*n+j) = uK12dint - lK12dint;
        Kmatrix(n-1+i,j) = uK21aint - lK21aint;
        Kmatrix(n-1+i,n+j) = uK21bint - lK21bint;
        Kmatrix(n-1+i,2*n+j) = uK22cint - lK22cint;
        Kmatrix(n-1+i,3*n+j) = uK22dint - lK22dint;
    end
end

%deals with the conditions at infinity
ginfty(1,n) = 1;
ginfty(2,2*n) = 1;

hinfty(1,3*n) = 1;
hinfty(2,4*n-1) = -1;
hinfty(2,4*n) = 1;

%creates the full matrix
fmatrix = zeros(4*n,4*n);
fmatrix(1:2*(n-1),:) = Ctymatrix;
fmatrix(2*(n-1)+1:4*(n-1),:) = Kmatrix;
fmatrix(4*n-3:4*n-2,:) = ginfty;
fmatrix(4*n-1:4*n,:) = hinfty;

%creates the RHS of the matrix equation
inhomog = zeros(4*n,1);
inhomog(4*n-2,1) = compress;

%solves the matrix equation
soln = fmatrix\inhomog;        
        
        
        
        
        
        
        
        
        
        
        
    
