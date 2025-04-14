function save_image(field,savename,ext,colormap,scale,NaN_col)
% Saves the given field as an image 'savename.ext' using specified colormap
% - field: field(x,y), with x as horizontal axis (column)
% - savename: string, name of saved image, excluding file extension
% - ext: string, file extension
% - colormap: colormap, matrix of size [N 3] describing colour scale
% - scale: vector [v1 v2], fix the ends of the colorbar to these values
% - NaN_col: color of NaN values, [R G B] vector

% ----------------------------------------------------------------------------
% Note: Enter field' for traditional matrix approach with first index as row.
% ----------------------------------------------------------------------------

if nargin < 5; scale = 0; end
if nargin < 6; NaN_col = 0; end
if NaN_col == 0; NaN_col = [0.5 0.5 0.5]; end

if numel(scale) == 1
    field = squeeze(field/max(max(abs(field))))';
else
    field = squeeze(2*(field'-scale(1))/(scale(2)-scale(1))-1);
end

field = field(end:-1:1,:);

if sum(isnan(field),'all') == 0

    imwrite(field*128+128,colormap,[savename '.' ext]);

else
    
    field_rgb = interp1(linspace(0,1,length(colormap)),colormap,field/2+0.5);
    NaN_num = sum(isnan(field),'all');
    field_rgb(isnan(field_rgb)) = reshape(ones(NaN_num,1)*reshape(NaN_col,[1 3]),[],1);

    imwrite(field_rgb,[savename '.' ext]);

end

end

