n_val = [639, 720, 830, 962, 1124, 1286, 1440, 1661];
x_val = [893, 906, 946, 921, 934, 944, 893, 940];
s = 0.138673;

for k = 1:numel(n_val)
    file = strcat('n', num2str(n_val(k)), 'x', num2str(x_val(k)));
    load(file)
    [kernel_matrix, ~] = pressure_shear_matrix_s(x,z,s,t);
    save(file, '-append', 'kernel_matrix','s');
end