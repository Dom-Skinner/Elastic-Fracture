function scaling_2 = convert(s1,s2,n,t,x,scaling_1)
% Converts a vector 'scaling_1' i.e. think hprime or gprime into
% another vector with different scalings. 
% I.e unscaled(1:t-1) = x(1:t-1).^(s1-1).*scaling_1 =
% x(1:t-1).^(s2-1).*scaling_2 
%
% N.B. a purely linear scaling gets a scale of s1=1

%
scaling_2 = zeros(n,1);
scaling_2(1:t-1) = x(1:t-1).^(s1-s2).*scaling_1(1:t-1);
scaling_2(t:n) = scaling_1(t:n);
if s1 < s2 
    disp('Warning, could be singularities in convert function')
    scaling_2(1) = scaling_2(2); % Better than infinity right?...
end
end