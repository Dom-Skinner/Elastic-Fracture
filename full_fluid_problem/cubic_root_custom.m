function w = cubic_root_custom(a,b,c)
% Finds the roots of the cubic polynomial ax^3 + bx + c
% Faster than matlab, since we know its a (special) cubic

d0 = -3*a.*b;
if ~all(d0) 
    disp('Error in custom root function')
end
d1 = 27*a.^2.*c;
unit_root = [1 ; (-0.5 + 0.5*sqrt(-3)); (-0.5 - 0.5*sqrt(-3))];
C = unit_root*((d1 + sqrt(d1.^2 - 4*d0.^3))/2).^(1/3);
w(1,:) = -(1/3./a).*(C(1,:) + d0./C(1,:));
w(2,:) = -(1/3./a).*(C(2,:) + d0./C(2,:));
w(3,:) = -(1/3./a).*(C(3,:) + d0./C(3,:));
end