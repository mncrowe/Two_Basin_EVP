function n = pad_zeros(n,d)
% Outputs a string with extra leading spaces replaced with zeros, e.g. 13 -> '00013'
% - n: number
% - d: number of digits

n = num2str(n);
if length(n)<d
    m = d-length(n);
    for im = 1:m
        n = ['0' n];
    end

end