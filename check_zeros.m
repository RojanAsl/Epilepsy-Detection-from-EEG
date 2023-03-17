% MATRIX --> Matrix that has to get cleaned
% MINIMUM -> normalization of the matrix numbers (5 in the other examples) 
% n = iteration size

function [clean_matrix] = check_zeros(matrix,minimum,n)
%UNTITLED Summary of this function goes here
%   matrix is the original matrix introduced
%  
int_matrix = matrix; %program works here
for e = 1:length(int_matrix)
    f = floor(abs(int_matrix(e)));
    if f <= minimum
        int_matrix(e) = 0;
    else
        int_matrix(e) = f;
    end
end

M = 0;
while (M == 0)
    M = median(int_matrix(1:n));
    if (M == 0)
        int_matrix = int_matrix(n:end);
##        disp('deleting...')
    end
end
%finding the index where to cut
length_int = length(int_matrix);
length_mat = length(matrix); 
cutting_index = length_mat - length_int;
clean_matrix = matrix(cutting_index:end);

end