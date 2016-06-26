%creates a row to measure the pressure at point z
%takes in as input n points

function pressure_row = diagonal_pressure_row_h(x,z)

n = length(x);
t = round(n/2);

infinity = 10^(10);


%matrix for the interpolation
%also we have linear panels up to panel nfin.
nfin = n-1;

%also must add special coefficients at infinity. put these at the end.
interpol = zeros(2*nfin+4,n);

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
hcoeffmatrix = [x(n-3), 1, 1/x(n-3), 1/x(n-3)^2; x(n-2), 1, 1/x(n-2), 1/x(n-2)^2; ...
        x(n-1), 1, 1/x(n-1), 1/x(n-1)^2; x(n), 1, 1/x(n), 1/x(n)^2]^(-1);
interpol(2*nfin+1:2*nfin+4,n-3:n) = hcoeffmatrix;


%now for the actual entries
kernel_row = zeros(1,2*nfin+4);

    %does the 1/sqrt(x) panels
    for j = 1:t-1
        lx = x(j);
        ux = x(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatricesadjusted(lx,z);
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatricesadjusted(ux,z);
            
        kernel_row(j) = uK12cint - lK12cint;
        kernel_row(nfin+j) = uK12dint - lK12dint;
    end
    
    %does the linear panels
    for j = t:nfin
        lx = x(j);
        ux = x(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,z);
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,z);
            
        kernel_row(j) = uK12cint - lK12cint;
        kernel_row(nfin+j) = uK12dint - lK12dint;
    end
    
        %does the panel at infinity
    lx = x(nfin+1);
    ux = infinity;
    [lK11aint, lK11bint, lK11eint, lK11gint, lK11iint, ...
    lK12cint, lK12dint, lK12fint, lK12hint, lK12jint, ...
    lK21aint, lK21bint, lK21eint, lK21gint, lK21iint, ...
    lK22cint, lK22dint, lK22fint, lK22hint, lK22jint] = quadinvandinvsquarekernelmatrices(lx,z);

    [uK11aint, uK11bint, uK11eint, uK11gint, uK11iint, ...
    uK12cint, uK12dint, uK12fint, uK12hint, uK12jint, ...
    uK21aint, uK21bint, uK21eint, uK21gint, uK21iint, ...
    uK22cint, uK22dint, uK22fint, uK22hint, uK22jint] = quadinvandinvsquarekernelmatrices(ux,z);

    kernel_row(2*nfin+1) = uK12cint - lK12cint;
    kernel_row(2*nfin+2) = uK12dint - lK12dint;
    kernel_row(2*nfin+3) = uK12hint - lK12hint;
    kernel_row(2*nfin+4) = uK12jint - lK12jint;

pressure_row = kernel_row*interpol;

return
end