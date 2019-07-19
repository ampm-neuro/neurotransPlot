function matrix_out = sheets_2d(matrix_in, dim, num_items)
%creates sheets size num_items along dimension dim to add a items to ~dim
%matrix. size(matrix,dim) must be a multiple of num_items.

%check input
if rem(size(matrix_in,dim), num_items) > 0
    error('incorrect size inputs')
end

matrix_out = [];

%iterate and load
idx_min = 1;
for isheet = 1:(size(matrix_in,dim)/num_items)
    
    idx_max = num_items*isheet;
    
    if dim == 1
        matrix_out = [matrix_out matrix_in(idx_min:idx_max, :, :)];
    elseif dim == 2
        matrix_out = [matrix_out; matrix_in(:, idx_min:idx_max)];
    end
    
    
    
    idx_min = num_items*isheet + 1
    
end


end