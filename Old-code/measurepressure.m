function [pressurerow, shearrow] = measurepressure(x, inpoint)

n = length(inpoint);
t = round(n/2);
nfin = n-1;
%small windows around crack tip
epsilon = 10^(-10);
%the value we integrate up to
infinity = 10^(20);

pressurerow = zeros(1,4*nfin+10);
shearrow = zeros(1,4*nfin+10);



%does the 1/sqrt(x) panels
    for j = 2:t-1
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatrices(ux,x);
            
        pressurerow(j) = uK11aint - lK11aint;
        pressurerow(nfin+j) = uK11bint - lK11bint;
        pressurerow(2*nfin+j) = uK12cint - lK12cint;
        pressurerow(3*nfin+j) = uK12dint - lK12dint;
        shearrow(j) = uK21aint - lK21aint;
        shearrow(nfin+j) = uK21bint - lK21bint;
        shearrow(2*nfin+j) = uK22cint - lK22cint;
        shearrow(3*nfin+j) = uK22dint - lK22dint;
    end
    
    %does the linear panels
    for j = t:nfin
        lx = inpoint(j);
        ux = inpoint(j+1);

        [lK11aint, lK11bint, lK12cint, lK12dint, ...
            lK21aint, lK21bint, lK22cint, lK22dint] = linkernelmatrices(lx,x);
        [uK11aint, uK11bint, uK12cint, uK12dint, ...
            uK21aint, uK21bint, uK22cint, uK22dint] = linkernelmatrices(ux,x);
            
        pressurerow(j) = uK11aint - lK11aint;
        pressurerow(nfin+j) = uK11bint - lK11bint;
        pressurerow(2*nfin+j) = uK12cint - lK12cint;
        pressurerow(3*nfin+j) = uK12dint - lK12dint;
        shearrow(j) = uK21aint - lK21aint;
        shearrow(nfin+j) = uK21bint - lK21bint;
        shearrow(2*nfin+j) = uK22cint - lK22cint;
        shearrow(3*nfin+j) = uK22dint - lK22dint;
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

    %pressurerow( 4*nfin+1) = uK11eint - lK11eint;
    pressurerow( 4*nfin+2) = uK11aint - lK11aint;
    pressurerow( 4*nfin+3) = uK11bint - lK11bint;
    pressurerow( 4*nfin+4) = uK11gint - lK11gint;
    pressurerow( 4*nfin+5) = uK11iint - lK11iint;
    
    pressurerow( 4*nfin+6) = uK12fint - lK12fint;
    pressurerow( 4*nfin+7) = uK12cint - lK12cint;
    pressurerow( 4*nfin+8) = uK12dint - lK12dint;
    pressurerow( 4*nfin+9) = uK12hint - lK12hint;
    pressurerow( 4*nfin+10) = uK12jint - lK12jint;

    %shearrow( 4*nfin+1) = uK21eint - lK21eint;
    shearrow( 4*nfin+2) = uK21aint - lK21aint;
    shearrow( 4*nfin+3) = uK21bint - lK21bint;
    shearrow( 4*nfin+4) = uK21gint - lK21gint;
    shearrow( 4*nfin+5) = uK21iint - lK21iint;
    
    shearrow( 4*nfin+6) = uK22fint - lK22fint;
    shearrow( 4*nfin+7) = uK22cint - lK22cint;
    shearrow( 4*nfin+8) = uK22dint - lK22dint;
    shearrow( 4*nfin+9) = uK22hint - lK22hint;
    shearrow( 4*nfin+10) = uK22jint - lK22jint;

    
    %does panel at the crack itself
    lx = epsilon;
    ux = inpoint(2);
    [lK11aint, lK11bint, lK12cint, lK12dint, ...
        lK21aint, lK21bint, lK22cint, lK22dint] = sqrtkernelmatrices(lx,x);
    [uK11aint, uK11bint, uK12cint, uK12dint, ...
        uK21aint, uK21bint, uK22cint, uK22dint] = sqrtkernelmatrices(ux,x);
    
    %we use the same linear interpolation as for the first panel
    j = 1;
    pressurerow(j) = uK11aint - lK11aint;
    pressurerow(nfin+j) = uK11bint - lK11bint;
    pressurerow(2*nfin+j) = uK12cint - lK12cint;
    pressurerow(3*nfin+j) = uK12dint - lK12dint;
    shearrow(j) = uK21aint - lK21aint;
    shearrow(nfin+j) = uK21bint - lK21bint;
    shearrow(2*nfin+j) = uK22cint - lK22cint;
    shearrow(3*nfin+j) = uK22dint - lK22dint;
    
return
end