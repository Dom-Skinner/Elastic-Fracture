n_val = [350,407,465,524,605,698,815];
x_val = [873,822,819,846,833,832,846];

for l = 1:numel(n_val)
    %{
    clearvars '-except' l x_val n_val
    nstr = num2str(n_val(l));
    xstr = num2str(x_val(l));
    file = strcat('n',nstr,'x',xstr);
    load(file)
    clearvars '-except'  x z n t nstr xstr l x_val n_val
    xmax = x(end);
    lambda = 0.057:0.0001:0.0589;
    scaled_K_of_c_march
    file = strcat('n',nstr,'x',xstr,'-mod');
    save(file)
    %}
    clearvars '-except' l x_val n_val
    nstr = num2str(n_val(l));
    xstr = num2str(x_val(l));
    file = strcat('n',nstr,'x',xstr,'-mod');
    load(file)
    file = strcat('n',nstr,'x',xstr,'-mod');
    save(file, '-append', 'kernel_matrix')
end
