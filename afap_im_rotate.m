function [imrot,err] = afap_im_rotate(im,offset,fast)
% Rotates afap im images such that the fish is oriented anterior to the top
%
% Pedro Henriques, Oct 2017
% Modified, Dec 2018

if nargin < 3
    fast = 0;
    if nargin < 2
        offset = 90;
    end
end

err = 0;
eyethr = 0.3;
bodythr = 0.01;
bodythr2 = 0.1;

%%

h = waitbar(0,'Init');
imrot = zeros(size(im),'uint8');
for f = 1:size(im,3)
    waitbar(f/size(im,3),h,'Rotating')
    I = im(:,:,f);
    if fast
        bodybw2 = imbinarize(I,bodythr2);
        s = regionprops(bodybw2,{'Area','Orientation'});
        [~,mix] = max([s.Area]);
        ori = s(mix).Orientation;
        if ori < 0 && ori < 90
            rori = -ori-90;
        elseif ori > 0 && ori < 90
            rori = ori-90;
        end
        imrot(:,:,f) = imrotate(I,rori,'bicubic','crop');
    else
        eyesBW = imfill(bwmorph( ...
            imbinarize(I,eyethr),'open'),'holes');
        
        s = regionprops(eyesBW,'Centroid','Area');
        [~,eix] = sort([s.Area],'descend');
        s = s(eix);
        
        if length(s) >= 3
            s = s(1:3); % Process only the 3 largest objects
            p = combnk(1:3,2);   % Combinations of objects to calculate distances
            eye_c = vertcat(s.Centroid);   % Object centroids
            
            dd = NaN(3,1);   % Object distance vector
            for j = 1:3
                dd(j) = abs(sqrt((diff([eye_c(p(j,1),1),eye_c(p(j,2),1)]).^2) + ...
                    ((diff([eye_c(p(j,1),2),eye_c(p(j,2),2)]).^2))));    % Calculate distances between each object centroid
            end
            [~, ix] = min(dd);  % Minumum distance should correspond to distance between both eyes
            eye_ix = p(ix,:)';   % Eyes ix
            sb_ix = setdiff(p,eye_ix);  % SB ix is the other
            
            eye_xy1 = eye_c(eye_ix(1),:);   % eye 1 centroid
            eye_xy2 = eye_c(eye_ix(2),:);   % eye 2 centroid
            sb_cxy = eye_c(sb_ix,:); % sb centroid
            cp = [(eye_xy1(1)+eye_xy2(1))/2, (eye_xy1(2)+eye_xy2(2))/2];    % Coordenate of middle point between both eyes
            hvec = [cp(1)-sb_cxy(1), cp(2)-sb_cxy(2), 0];   % Vector between sb and eye middle point
            h_ori = atan2d(hvec(2),hvec(1));
            
            imrot(:,:,f) = imrotate(I,h_ori+offset,'bilinear','crop');
        elseif length(s) == 2
            % Use body centroid to calculate headind
            
            eye_xy1 = s(1).Centroid;   % eye 1 centroid
            eye_xy2 = s(2).Centroid;   % eye 2 centroid
            cp = [(eye_xy1(1)+eye_xy2(1))/2, (eye_xy1(2)+eye_xy2(2))/2];
            
            bodybw = imbinarize(I,bodythr);
            s2 = regionprops(bodybw,'Area','Centroid');
            [~,eix] = sort([s2.Area],'descend');
            s2 = s2(eix(1));
            body_cxy = s2.Centroid;
            hvec = [cp(1)-body_cxy(1), cp(2)-body_cxy(2), 0];   % Vector between sb and eye middle point
            h_ori = atan2d(hvec(2),hvec(1));
            
            imrot(:,:,f) = imrotate(I,h_ori+offset,'bicubic','crop');
        else
            disp('Cound''t find orientation')
            err = 1;
        end
    end
end
close(h)
end