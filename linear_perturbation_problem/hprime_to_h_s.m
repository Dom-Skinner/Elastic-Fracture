%this function, for a given h', analytically integrates to give h

%produces a matrix, which when multiplied by the n values of h', in the 
%s representation, gives 3n coefficients for h.

function h_coeff_mat = hprime_to_h_s(x,s,t)

n = length(x);

%gets the coefficients of g', h'
interpolate_matrix = linear_simpleinfty_interpolate_s(x,s);
interpolate_matrix = interpolate_matrix(2*n+1:4*n,n+1:2*n);

h_coeff_mat = zeros(3*n,2*n);


%gets the a and b coefficients - easy
for j=1:n
    if j <= t-1
        h_coeff_mat(3*(j-1)+1,j) = 1/(s+1);
        h_coeff_mat(3*(j-1)+2,n+j) = 1/s;
    else
        h_coeff_mat(3*(j-1)+1,j) = 1/2;
        h_coeff_mat(3*(j-1)+2,n+j) = 1;
    end
end
%gets the c coefficients.
% now vectorised for speed at some readibility cost...
x_diff_save_s_1 = (1/(s+1))*(x(2:t).^(s+1) - x(1:t-1).^(s+1));
x_diff_save_s_2 = (1/s)*(x(2:t).^(s) - x(1:t-1).^(s));

for j = 2:t-1
     h_coeff_mat(3*j,j) = -(1/(s+1))*x(j)^(s+1);
     h_coeff_mat(3*j,n+j) = -(1/s)*x(j)^(s);
     %
     i = 1:j-1;
     h_coeff_mat(3*j,i) = x_diff_save_s_1(i);
     h_coeff_mat(3*j,n+i) = x_diff_save_s_2(i);
    
end
for j = t:n
      h_coeff_mat(3*j,j) = -(1/2)*x(j)^2;
      h_coeff_mat(3*j,n+j) = -x(j);
   
      h_coeff_mat(3*j,1:t-1) = x_diff_save_s_1;
      h_coeff_mat(3*j,n+1:t-1+n) = x_diff_save_s_2;
      if j > t
          i = t:j-1;
          h_coeff_mat(3*j,i) = (1/2)*(x(i+1).^2 - x(i).^2);
          h_coeff_mat(3*j,n+i) = x(i+1) - x(i);
      end
end
%{
%gets the c coefficients.
for j=2:n
    for i = 1:j-1
        if i <= t-1
            h_coeff_mat(3*j,i) = (1/(s+1))*(x(i+1)^(s+1) - x(i)^(s+1));
            h_coeff_mat(3*j,n+i) = (1/s)*(x(i+1)^(s) - x(i)^(s));
        else
            h_coeff_mat(3*j,i) = (1/2)*(x(i+1)^2 - x(i)^2);
            h_coeff_mat(3*j,n+i) = x(i+1) - x(i);
        end
    end

    if j <= t-1
        h_coeff_mat(3*j,j) = -(1/(s+1))*x(j)^(s+1);
        h_coeff_mat(3*j,n+j) = -(1/s)*x(j)^(s);
    else
        h_coeff_mat(3*j,j) = -(1/2)*x(j)^2;
        h_coeff_mat(3*j,n+j) = -x(j);
    end
end
%}
%
% The following lines are equivalent to:
%h_co_mat = h_coeff_mat*interpolate_matrix;
% But we know that the interpolate matrix is essentially diagonal, so we
% can save ourselves time by doing the multiplication ourselves
h_co_mat(:,1) = interpolate_matrix(1,1)*h_coeff_mat(:,1)+...
    + interpolate_matrix(n+1,1)*h_coeff_mat(:,n+1);
for k = 2:n-2
    h_co_mat(:,k) = interpolate_matrix(k,k)*h_coeff_mat(:,k)+...
                    interpolate_matrix(k-1,k)*h_coeff_mat(:,k-1)+...
                    interpolate_matrix(n+k,k)*h_coeff_mat(:,n+k)+...
                    interpolate_matrix(n+k-1,k)*h_coeff_mat(:,n+k-1);
end
h_co_mat(:,n-1) = interpolate_matrix(n-2,n-1)*h_coeff_mat(:,n-2)+...
                  interpolate_matrix(n-1,n-1)*h_coeff_mat(:,n-1)+...
                  interpolate_matrix(n,n-1)*h_coeff_mat(:,n)+...
                  interpolate_matrix(2*n-2,n-1)*h_coeff_mat(:,2*n-2)+...
                  interpolate_matrix(2*n-1,n-1)*h_coeff_mat(:,2*n-1)+...
                  interpolate_matrix(2*n,n-1)*h_coeff_mat(:,2*n);

h_co_mat(:,n) = interpolate_matrix(n-1,n)*h_coeff_mat(:,n-1)+...
                interpolate_matrix(n,n)*h_coeff_mat(:,n)+...
                interpolate_matrix(2*n-1,n)*h_coeff_mat(:,2*n-1)+...
                interpolate_matrix(2*n,n)*h_coeff_mat(:,2*n);
                

h_coeff_mat = h_co_mat;