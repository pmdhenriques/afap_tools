function O = afap_im_segment_fish(im,eyethr,bodythr,tailanaly)
% Segment a fish from an image frame from the afap pipeline and return
% object properties (eyes, orientation, swim-bladder and tail points).
%
% Should work with other fish images but the proper thresholds need be
% applied.
%
% Pedro Henriques, Oct 2017

properties = {'Area','Centroid','Eccentricity','MajorAxisLength', ...
    'MinorAxisLength','Orientation'};

if nargin < 4
    tailanaly = 0;
    if nargin < 3
        bodythr = 0.01;
        if nargin < 2
            eyethr = 0.3;
        end
    end
end

eyesBW = imfill(bwmorph( ...
    imbinarize(im,eyethr),'open'),'holes');

stats = regionprops(eyesBW,properties);
[~,eix] = sort([stats.Area],'descend');
stats = stats(eix);

bodybw = imbinarize(im,bodythr);
stats_body = regionprops(bodybw,properties);
[~,eix] = sort([stats_body.Area],'descend');
stats_body = stats_body(eix(1));

% Tail parameters
npts = 11; % number of segment points
seglen = 7; % segment length (px)

if length(stats) >= 3
    stats = stats(1:3); % Process only the 3 largest objects
    bcp = combnk(1:3,2);   % Combinations of objects to calculate distances
    eye_c = vertcat(stats.Centroid);   % Object centroids
    
    dd = pdist(eye_c);  % Eucl distance between each object
    [~, ix] = min(dd);  % Minumum distance should correspond to distance between both eyes
    eye_ix = bcp(ix,:)';   % Eyes ix
    sb_ix = setdiff(bcp,eye_ix);  % SB ix is the other
    
    eye_xy1 = eye_c(eye_ix(1),:);   % eye 1 centroid
    eye_xy2 = eye_c(eye_ix(2),:);   % eye 2 centroid
    
    sbcp = eye_c(sb_ix,:); % sb centroid
    cp = [(eye_xy1(1)+eye_xy2(1))/2, (eye_xy1(2)+eye_xy2(2))/2];    % Coordenate of middle point between both eyes
    hvec = [cp(1)-sbcp(1), cp(2)-sbcp(2), 0];   % Vector between sb and eye middle point
    
elseif length(stats) == 2
    % Use body centroid to calculate heading
    
    eye_ix = 1:2;
    eye_xy1 = stats(1).Centroid;   % eye 1 centroid
    eye_xy2 = stats(2).Centroid;   % eye 2 centroid
    cp = [(eye_xy1(1)+eye_xy2(1))/2, (eye_xy1(2)+eye_xy2(2))/2];
    
    sbcp = stats_body.Centroid; % Use body centroid instead
    hvec = [cp(1)-sbcp(1), cp(2)-sbcp(2), 0];
    stats = cat(1,stats,stats_body);
    
else
    disp('Eyes cannot be found');
    O = [];
    return
end

h_ori = atan2d(hvec(2),hvec(1));    % Orientation of sb-eye vector

% Which eye?

eyeu = [eye_xy1(1)-eye_xy2(1), eye_xy1(2)-eye_xy2(2)];  % eye-eye vector
eye_head_ang = atan2d(hvec(1)*eyeu(2)-hvec(2)*eyeu(1), ...
    hvec(1)*eyeu(1)+hvec(2)*eyeu(2));

if eye_head_ang >= 0
    stats(eye_ix(1)).Side = 1;  % Right eye
    stats(eye_ix(2)).Side = 0;  % Left eye
else
    stats(eye_ix(1)).Side = 0;
    stats(eye_ix(2)).Side = 1;
end

% Make sure there are no empty cells so no problem later in concatenating
% structure
ix = find(cellfun(@isempty,{stats.Side}));
for i = ix
    stats(i).Side = nan;
end

%%

if length(stats) >= 3
    lix = [stats.Side] == 0;
    rix = [stats.Side] == 1;
    cl = stats(lix).Centroid;
    cr = stats(rix).Centroid;
    ol = stats(lix).Orientation;
    or = stats(rix).Orientation;
    mal = stats(lix).MajorAxisLength;
    mar = stats(rix).MajorAxisLength;
    
    % Major axis extremes
    extl1 = [cl(1)+(cos(deg2rad(ol)).*mal./2), ...
        cl(2)-(sin(deg2rad(ol)).*mal./2)];
    extl2 = [cl(1)-(cos(deg2rad(ol)).*mal./2), ...
        cl(2)+(sin(deg2rad(ol)).*mal./2)];
    extr1 = [cr(1)+(cos(deg2rad(or)).*mar./2), ...
        cr(2)-(sin(deg2rad(or)).*mar./2)];
    extr2 = [cr(1)-(cos(deg2rad(or)).*mar./2), ...
        cr(2)+(sin(deg2rad(or)).*mar./2)];
    
    extl = [extl1;extl2];
    extr = [extr1;extr2];        
    
    %
    bcp = stats_body.Centroid;
    ddl = abs(sqrt((extl(:,1)-bcp(1)).^2 + (extl(:,2)-bcp(2)).^2));
    ddr = abs(sqrt((extr(:,1)-bcp(1)).^2 + (extr(:,2)-bcp(2)).^2));
    lixx = false(2,1);
    rixx = false(2,1);
    [~, ddlm] = max(ddl);   lixx(ddlm) = 1;
    [~, ddrm] = max(ddr);   rixx(ddrm) = 1;
    
    ul = [extl(lixx,1)-extl(~lixx,1), extl(lixx,2)-extl(~lixx,2), 0];    % Left eye vector
    ur = [extr(rixx,1)-extr(~rixx,1), extr(rixx,2)-extr(~rixx,2), 0];
    
    %
    
    stats(lix).ext_lower = extl(~lixx,:);
    stats(lix).ext_higher = extl(lixx,:);
    stats(lix).eye_vector = ul;
    stats(rix).ext_lower = extr(~rixx,:);
    stats(rix).ext_higher = extr(rixx,:);
    stats(rix).eye_vector = ur;
end

O.stats = stats;
O.stats_body = stats_body;
O.cp = cp;
O.hvec = hvec;
O.h_ori = h_ori;

%% Get tail segments

if tailanaly
    bodydbw = bwdist(~bodybw);
    
    Tpts = NaN(npts,2);
    Tpts(1,:) = sbcp;  % First point is sb centroid
    for i = 1:npts-1
        if i == 1
            cx = seglen*2*cos(0:0.1:2*pi)+Tpts(i,1);  % Draw a circle
            cy = seglen*2*sin(0:0.1:2*pi)+Tpts(i,2);
        else
            cx = seglen*cos(0:0.1:2*pi)+Tpts(i,1);
            cy = seglen*sin(0:0.1:2*pi)+Tpts(i,2);
        end
        [cxix,cyix,cI] = improfile(bodydbw,cx,cy,'bilinear');   % Build intencity profile of the circle perimeter
        cI2 = circshift(cI,round(length(cI)/4));    % shift intensity profile to correct for peaks that are in the extremities
        warning('off','all')
        [pks,locs] = findpeaks(cI,'MinPeakDistance',pi*seglen,'MinPeakHeight',1);  % Find peaks
        [pks2,~] = findpeaks(cI2,'MinPeakDistance',pi*seglen,'MinPeakHeight',1);   % Find peaks in the shifted vector
        warning('on','all')
        
        if length(pks) < 2 && length(pks2) >= 2     % Add a peak in the extremetry
            pks = pks2; locs = [1; locs];
        end
        if length(pks) > 2
            % Use only two largest peaks
            [pks, ix] = sort(pks,'descend');
            pks = pks(1:2);
            locs = locs(ix(1:2));
        end
        if length(pks) == 2
            if i == 1
                dd = abs(sqrt((cxix(locs)-cp(1)).^2 + (cyix(locs)-cp(2)).^2));  % Second point is defined by the maximum distance to the eye-eye center point
            else
                dd = NaN(length(pks),1);    % Next iterations look for points that maximize the distance to the previous point
                for j = 1:length(pks)
                    dd(j) = abs(sqrt((cxix(locs(j))-Tpts(i-1,1))^2 + (cyix(locs(j))-Tpts(i-1,2))^2));
                end
            end
            [~, ix] = max(dd);
            Tpts(i+1,:) = [cxix(locs(ix)),cyix(locs(ix))];
        else
            % Terminate for loop if cant find more than one
            % point
            break
        end
    end
    
    O.tailpts = Tpts;
end

end