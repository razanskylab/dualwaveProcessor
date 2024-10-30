function mip_clahe = clahe_to_mip(mip, clipLimit)
%apply CLAHE to MIP for better visualization

mip = single(mip);
oldMax = max(mip(:));
mip = mip ./ max(abs(mip(:)));
mip_clahe = adapthisteq(mip, 'clipLimit', clipLimit);
mip_clahe = mip_clahe ./ max(abs(mip_clahe(:)));
mip_clahe = mip_clahe .* oldMax;

end



