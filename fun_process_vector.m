function output_vector = fun_process_vector(input_vector)
%% Description:
% This function receives a 18-rows vector, rpunds the components smaller
% than 1e-10 to 0 (in abs value) and converts the rest of the
% components to scientific notation with 3 significant digits

if length(input_vector) ~= 18
    error('Input vector should have 18 rows.');
end
output_vector = zeros(size(input_vector));

for i = 1:length(input_vector)
    real_part = real(input_vector(i));
    imag_part = imag(input_vector(i));

    if abs(real_part) < 1e-10
        real_part = 0;
    else
        real_part = round(real_part, 3, 'significant');
    end

    if abs(imag_part) < 1e-10
        imag_part = 0;
    else
        imag_part = round(imag_part, 3, 'significant');
    end
    output_vector(i) = real_part + 1i * imag_part;
end

end