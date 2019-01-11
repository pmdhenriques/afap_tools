function [imoffset] = afaptdms_GetFrLocation_ph(thisid, frid)

imoffset = find(frid == thisid);
% safety feature:
if numel(imoffset)>1
    warndlg('Image duplicated!')
    imoffset = [];
elseif isempty(imoffset)
%     disp('Image not found! Finding nearest buffer...')
    imoffset = findnearest(frid,thisid,1);
    imoffset = imoffset(1);
end
end