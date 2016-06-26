function [pressurerow, shearrow] = pressure_cubic_recipe_alt_infty(z,x)

%creates two rows, for measuring the pressure and shear stress at a point
%z. x is understood to be inpoint
%uses the mixed cubic representation earlier.

%uses same representation at crack tip as in first panel
%uses same representation at infinity as in last panel

%takes as input just a vector of coefficients

n = length(x);
t = round(n/2);
%small windows around crack tip
epsilon = 10^(-10);
%the value we integrate up to
infinity = 10^(10);

pressurerow = zeros(1,8*(n-1)+9);
shearrow = zeros(1,8*(n-1)+9);

%does the x^(-1/2)*(ax^3+bx^2+cx+d) panels
for i = 1:t-1
    %handles the crack tip
    if i == 1
        lx = epsilon;
        ux = x(2);
    else
        lx = x(i);
        ux = x(i+1);
    end
    
    [lK11sqrtx3, lK11sqrtx2, lK11sqrtx1, lK11sqrtx0, ...
    lK12sqrtx3, lK12sqrtx2, lK12sqrtx1, lK12sqrtx0, ...
    lK21sqrtx3, lK21sqrtx2, lK21sqrtx1, lK21sqrtx0, ...
    lK22sqrtx3, lK22sqrtx2, lK22sqrtx1, lK22sqrtx0] = sqrt_x3_to_x0_kernel(lx,z);

    [uK11sqrtx3, uK11sqrtx2, uK11sqrtx1, uK11sqrtx0, ...
    uK12sqrtx3, uK12sqrtx2, uK12sqrtx1, uK12sqrtx0, ...
    uK21sqrtx3, uK21sqrtx2, uK21sqrtx1, uK21sqrtx0, ...
    uK22sqrtx3, uK22sqrtx2, uK22sqrtx1, uK22sqrtx0] = sqrt_x3_to_x0_kernel(ux,z);

    pressurerow(4*(i-1)+1) = uK11sqrtx3 - lK11sqrtx3;
    pressurerow(4*(i-1)+2) = uK11sqrtx2 - lK11sqrtx2;
    pressurerow(4*(i-1)+3) = uK11sqrtx1 - lK11sqrtx1;
    pressurerow(4*(i-1)+4) = uK11sqrtx0 - lK11sqrtx0;
    
    pressurerow(4*(n-1)+4*(i-1)+1) = uK12sqrtx3 - lK12sqrtx3;
    pressurerow(4*(n-1)+4*(i-1)+2) = uK12sqrtx2 - lK12sqrtx2;
    pressurerow(4*(n-1)+4*(i-1)+3) = uK12sqrtx1 - lK12sqrtx1;
    pressurerow(4*(n-1)+4*(i-1)+4) = uK12sqrtx0 - lK12sqrtx0;
    
    shearrow(4*(i-1)+1) = uK21sqrtx3 - lK21sqrtx3;
    shearrow(4*(i-1)+2) = uK21sqrtx2 - lK21sqrtx2;
    shearrow(4*(i-1)+3) = uK21sqrtx1 - lK21sqrtx1;
    shearrow(4*(i-1)+4) = uK21sqrtx0 - lK21sqrtx0;
    
    shearrow(4*(n-1)+4*(i-1)+1) = uK22sqrtx3 - lK22sqrtx3;
    shearrow(4*(n-1)+4*(i-1)+2) = uK22sqrtx2 - lK22sqrtx2;
    shearrow(4*(n-1)+4*(i-1)+3) = uK22sqrtx1 - lK22sqrtx1;
    shearrow(4*(n-1)+4*(i-1)+4) = uK22sqrtx0 - lK22sqrtx0;
end

%does the (cx + d + e/x + f/x^2) panels

for i = t:n-1
    lx = x(i);
    ux = x(i+1);
    
    [~, lK11x1, lK11x0, lK11xmin1, lK11xmin2, ...
    ~, lK12x1, lK12x0, lK12xmin1, lK12xmin2, ...
    ~, lK21x1, lK21x0, lK21xmin1, lK21xmin2, ...
    ~, lK22x1, lK22x0, lK22xmin1, lK22xmin2] = x2_to_xminus2_kernel(lx,z);

    [~, uK11x1, uK11x0, uK11xmin1, uK11xmin2, ...
    ~, uK12x1, uK12x0, uK12xmin1, uK12xmin2, ...
    ~, uK21x1, uK21x0, uK21xmin1, uK21xmin2, ...
    ~, uK22x1, uK22x0, uK22xmin1, uK22xmin2] = x2_to_xminus2_kernel(ux,z);

    pressurerow(4*(i-1)+1) = uK11x1 - lK11x1;
    pressurerow(4*(i-1)+2) = uK11x0 - lK11x0;
    pressurerow(4*(i-1)+3) = uK11xmin1 - lK11xmin1;
    pressurerow(4*(i-1)+4) = uK11xmin2 - lK11xmin2;
    
    pressurerow(4*(n-1)+4*(i-1)+1) = uK12x1 - lK12x1;
    pressurerow(4*(n-1)+4*(i-1)+2) = uK12x0 - lK12x0;
    pressurerow(4*(n-1)+4*(i-1)+3) = uK12xmin1 - lK12xmin1;
    pressurerow(4*(n-1)+4*(i-1)+4) = uK12xmin2 - lK12xmin2;
    
    shearrow(4*(i-1)+1) = uK21x1 - lK21x1;
    shearrow(4*(i-1)+2) = uK21x0 - lK21x0;
    shearrow(4*(i-1)+3) = uK21xmin1 - lK21xmin1;
    shearrow(4*(i-1)+4) = uK21xmin2 - lK21xmin2;
    
    shearrow(4*(n-1)+4*(i-1)+1) = uK22x1 - lK22x1;
    shearrow(4*(n-1)+4*(i-1)+2) = uK22x0 - lK22x0;
    shearrow(4*(n-1)+4*(i-1)+3) = uK22xmin1 - lK22xmin1;
    shearrow(4*(n-1)+4*(i-1)+4) = uK22xmin2 - lK22xmin2;
end

%does the infinity panel
lx = x(n);
ux = infinity;

[~, lK11x1, lK11x0, lK11xmin1, lK11xmin2, ...
lK12x2, lK12x1, lK12x0, lK12xmin1, lK12xmin2, ...
~, lK21x1, lK21x0, lK21xmin1, lK21xmin2, ...
lK22x2, lK22x1, lK22x0, lK22xmin1, lK22xmin2] = x2_to_xminus2_kernel(lx,z);

[~, uK11x1, uK11x0, uK11xmin1, uK11xmin2, ...
uK12x2, uK12x1, uK12x0, uK12xmin1, uK12xmin2, ...
~, uK21x1, uK21x0, uK21xmin1, uK21xmin2, ...
uK22x2, uK22x1, uK22x0, uK22xmin1, uK22xmin2] = x2_to_xminus2_kernel(ux,z);

pressurerow(8*(n-1)+1) = uK11x1 - lK11x1;
pressurerow(8*(n-1)+2) = uK11x0 - lK11x0;
pressurerow(8*(n-1)+3) = uK11xmin1 - lK11xmin1;
pressurerow(8*(n-1)+4) = uK11xmin2 - lK11xmin2;
    
pressurerow(8*(n-1)+5) = uK12x2 - lK12x2;
pressurerow(8*(n-1)+6) = uK12x1 - lK12x1;
pressurerow(8*(n-1)+7) = uK12x0 - lK12x0;
pressurerow(8*(n-1)+8) = uK12xmin1 - lK12xmin1;
pressurerow(8*(n-1)+9) = uK12xmin2 - lK12xmin2;

shearrow(8*(n-1)+1) = uK21x1 - lK21x1;
shearrow(8*(n-1)+2) = uK21x0 - lK21x0;
shearrow(8*(n-1)+3) = uK21xmin1 - lK21xmin1;
shearrow(8*(n-1)+4) = uK21xmin2 - lK21xmin2;
    
shearrow(8*(n-1)+5) = uK22x2 - lK22x2;
shearrow(8*(n-1)+6) = uK22x1 - lK22x1;
shearrow(8*(n-1)+7) = uK22x0 - lK22x0;
shearrow(8*(n-1)+8) = uK22xmin1 - lK22xmin1;
shearrow(8*(n-1)+9) = uK22xmin2 - lK22xmin2;

return
end


