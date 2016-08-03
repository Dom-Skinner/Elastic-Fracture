%this function, for a given h', analytically integrates to give h
%use 1/sqrt(x) and linear representation for h'

%produces a matrix, which when multiplied by the n values of h',
%gives 3n coefficients for h.

function h_co_mat = hprime_to_h(x,t,lambda)

n = length(x);

%gets the coefficients of g', h'
interpolate_matrix = linear_simpleinfty_interpolate(x,t,lambda);
interpolate_matrix = interpolate_matrix(2*n+1:4*n,n+1:2*n);

h_coefficient_matrix = zeros(3*n,2*n);


%gets the a and b coefficients - easy
for j=1:n
    if j <= t-1
        h_coefficient_matrix(3*(j-1)+1,j) = 2/3;
        h_coefficient_matrix(3*(j-1)+2,n+j) = 2;
    else
        h_coefficient_matrix(3*(j-1)+1,j) = 1/2;
        h_coefficient_matrix(3*(j-1)+2,n+j) = 1;
    end
end
%gets the c coefficients.
% now vectorised for speed at some readibility cost...
x_diff_save_sqrt_1 = (2/3)*(x(2:t).^(3/2) - x(1:t-1).^(3/2));
x_diff_save_sqrt_2 = 2*(x(2:t).^(1/2) - x(1:t-1).^(1/2));

for j = 2:t-1
     h_coefficient_matrix(3*j,j) = -(2/3)*x(j)^(3/2);
     h_coefficient_matrix(3*j,n+j) = -2*x(j)^(1/2);
     %
     i = 1:j-1;
     h_coefficient_matrix(3*j,i) = x_diff_save_sqrt_1(i);
     h_coefficient_matrix(3*j,n+i) = x_diff_save_sqrt_2(i);
    
end
for j = t:n
      h_coefficient_matrix(3*j,j) = -(1/2)*x(j)^2;
      h_coefficient_matrix(3*j,n+j) = -x(j);
   
      h_coefficient_matrix(3*j,1:t-1) = x_diff_save_sqrt_1;
      h_coefficient_matrix(3*j,n+1:t-1+n) = x_diff_save_sqrt_2;
      if j > t
          i = t:j-1;
          h_coefficient_matrix(3*j,i) = (1/2)*(x(i+1).^2 - x(i).^2);
          h_coefficient_matrix(3*j,n+i) = x(i+1) - x(i);
      end
end

%h_coefficient_matrix = h_coefficient_matrix*interpolate_matrix;
% To save time, we do the sparse calculation manually
h_co_mat(:,1) = interpolate_matrix(1,1)*h_coefficient_matrix(:,1)+...
    + interpolate_matrix(n+1,1)*h_coefficient_matrix(:,n+1);
for k = 2:n-2
    h_co_mat(:,k) = interpolate_matrix(k,k)*h_coefficient_matrix(:,k)+...
                    interpolate_matrix(k-1,k)*h_coefficient_matrix(:,k-1)+...
                    interpolate_matrix(n+k,k)*h_coefficient_matrix(:,n+k)+...
                    interpolate_matrix(n+k-1,k)*h_coefficient_matrix(:,n+k-1);
end
h_co_mat(:,n-1) = interpolate_matrix(n-2,n-1)*h_coefficient_matrix(:,n-2)+...
                  interpolate_matrix(n-1,n-1)*h_coefficient_matrix(:,n-1)+...
                  interpolate_matrix(n,n-1)*h_coefficient_matrix(:,n)+...
                  interpolate_matrix(2*n-2,n-1)*h_coefficient_matrix(:,2*n-2)+...
                  interpolate_matrix(2*n-1,n-1)*h_coefficient_matrix(:,2*n-1)+...
                  interpolate_matrix(2*n,n-1)*h_coefficient_matrix(:,2*n);

h_co_mat(:,n) = interpolate_matrix(n-1,n)*h_coefficient_matrix(:,n-1)+...
                interpolate_matrix(n,n)*h_coefficient_matrix(:,n)+...
                interpolate_matrix(2*n-1,n)*h_coefficient_matrix(:,2*n-1)+...
                interpolate_matrix(2*n,n)*h_coefficient_matrix(:,2*n);
                




