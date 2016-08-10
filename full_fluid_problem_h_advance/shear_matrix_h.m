function kernel_matrix = shear_matrix_h(x,z,t,lambda)

nx = length(x);
nz = length(z);
infinity = 10^10;

interpolate_matrix = linear_simpleinfty_interpolate_h(x,t,lambda);

function_matrix = zeros(nz,2*nx);


%FINDS ALL THE KERNEL FUNCTIONS

%finds the kernel K22

K22x1 = @(x,z) -real(4.*(4+3.*x.^2+(-4).*x.*z+z.^2).* ...
    ...
    (4+x.^2+(-2).*x.*z+z.^2).^(-2));

  
  
K22x0 = @(x,z) -real(8.*(4+(x+(-1).*z).^2).^(-2).*(x+(-1).*z));


%finds the kernel K22

  
K22sqrtx0 = @(x,z,at,at2) -real(2.*x.^(1/2).*(4+z.^2).^(-2).* ...
    (4+x.^2+(-2).*x.*z+z.^2).^(-2).*(6.* ...
    ...
    x.^3.*z+x.^2.*(4+(-23).*z.^2)+32.*x.*(z+z.^3)+(-5).*((-16)+8.* ...
    ...
    z.^2+3.*z.^4))+(sqrt(-1)*(3/2)).*((sqrt(-1)*(-2))+(-1).*z).^(-5/2) ...
    ...
    .*at2+(sqrt(-1)*(-3/2)).*((sqrt(-1)*2)+(-1).*z).^(-5/2).*at);


for i=1:nz
    %goes over the 1/sqrt(x) panels
    % All code vectorised for speed. Common expresions computed first.
    x_t = x(1:t);
    %
    atan_store_1 = atan(x_t.^(1/2).*((sqrt(-1)*2)+(-1)*z(i)).^(-1/2));
    atan_store_2 = atan(x_t.^(1/2).*((sqrt(-1)*(-2))+(-1)*z(i)).^(-1/2));
    %
    K22x0_store = K22x0(x_t,z(i));
    K22sqrtx0_store = K22sqrtx0(x_t,z(i),atan_store_1,atan_store_2);
    
    j=1:t-1;
    %for K22
    function_matrix(i,j) = K22x0_store(j+1)-K22x0_store(j);
    function_matrix(i,nx+j) = K22sqrtx0_store(j+1)-K22sqrtx0_store(j);
   
    %goes over the linear panels including infinity
    x_tmp = [x(t:nx) infinity];
    
    K22x1_store = K22x1(x_tmp,z(i));
    K22x0_store = K22x0(x_tmp,z(i));
    %
    j=1:nx-t+1;

    %for K22
    function_matrix(i,j+t-1)   = K22x1_store(j+1)-K22x1_store(j);
    function_matrix(i,nx+j+t-1) = K22x0_store(j+1)-K22x0_store(j);
  
end        


% The following lines are equivalent to:
%h_co_mat = h_coeff_mat*interpolate_matrix;
% But we know that the interpolate matrix is essentially diagonal, so we
% can save ourselves time by doing the multiplication ourselves

%%{
kernel_matrix = zeros(nz,nx);

kernel_matrix(:,1) = ...
    interpolate_matrix(1,1)*function_matrix(:,1)+...
    + interpolate_matrix(nx+1,1)*function_matrix(:,nx+1);
for k = 2:nx-2
    kernel_matrix(:,k) = ...
        interpolate_matrix(k,k)  *function_matrix(:,k)   + ...
        interpolate_matrix(k-1,k)*function_matrix(:,k-1) + ...
        interpolate_matrix(k+nx,k)  *function_matrix(:,nx+k) + ...
        interpolate_matrix(k-1+nx,k)*function_matrix(:,nx+k-1);
end
kernel_matrix(:,nx-1) = ...
    interpolate_matrix(nx-2,nx-1)*function_matrix(:,nx-2)   + ...
    interpolate_matrix(nx-1,nx-1)*function_matrix(:,nx-1)   + ...
    interpolate_matrix(nx,nx-1)  *function_matrix(:,nx)     + ...
    interpolate_matrix(2*nx-2,nx-1)*function_matrix(:,2*nx-2) + ...
    interpolate_matrix(2*nx-1,nx-1)*function_matrix(:,2*nx-1) + ...
    interpolate_matrix(2*nx,nx-1)  *function_matrix(:,2*nx);

kernel_matrix(:,nx) = ...
    interpolate_matrix(nx-1,nx)*function_matrix(:,nx-1)     + ...
    interpolate_matrix(nx,nx)  *function_matrix(:,nx)       + ...
    interpolate_matrix(2*nx-1,nx)*function_matrix(:,2*nx-1)   + ...
    interpolate_matrix(2*nx,nx)  *function_matrix(:,2*nx);
%}
return
end
    