
x = ones(3,8000);
for n = 1:4000
    x(:,n)=roots([1,0,3,2]);
    x(:,2*n)=cubic_root_custom([1,3,2]);
end