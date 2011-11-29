%/**
% Растяжение массива на N точек
%@param in_a - входной массив
%@param N - размер выходного массива 
%*/
function out_a = resize_array(in_a, N)

S = length(in_a);

if N <= S % Размер выходного массива должен быть больше или равен входного
    if N == S
        out_a = in_a;
    else
        out_a = NaN;
        return;
    end
else
    
    out_a = nan(1, N); % Выделение памяти под выходной массив
    for n = 1:N
        index_in = ceil(n/N*S);
        out_a(n) = in_a(index_in);
    end
    
end

end

