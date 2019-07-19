function matrix_out = sheets_3d(matrix_in, dim, num_items)
%creates sheets size num_items along dimension dim to add a 3rd dim to 2d
%matrix. size(matrix,dim) must be a multiple of num_items.

%check input
if rem(size(matrix_in,dim), num_items) > 0
    error('incorrect size inputs')
end

%preallocate
if dim == 1
    matrix_out = nan(num_items, size(matrix_in,2), size(matrix_in,1)/num_items);
elseif dim == 2
    matrix_out = nan(size(matrix_in,1), num_items, size(matrix_in,2)/num_items);
else
    error('dim input must be 1 or 2')
end

%iterate and load
idx_min = 1;
for isheet = 1:size(matrix_out,3)

    idx_max = num_items*isheet;
    
    if dim == 1
        matrix_out(:,:,isheet) = matrix_in(idx_min:idx_max, :, :);
    elseif dim == 2
         matrix_out(:,:,isheet) = matrix_in(:, idx_min:idx_max, :);
    end
    
    
    
    idx_min = num_items*isheet + 1;
    
end


end