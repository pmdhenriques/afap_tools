function [] = afap_btf(datadir,overwrite)
% Adds bout features to aq structure.
% Mostly based on gcbtf.m script.
%
% Pedro Henrques

if nargin < 2
    overwrite = 0;
elseif nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\NI project\Ablations', ...
        'Select fish to process');
end

h = waitbar(0,'Loading');
load(fullfile(datadir,'aq'));
if isfield(aq,'btf') && overwrite
    disp('Removing previous btf field')
    aq = rmfield(aq,'btf');
end

%% Initialize variables

tailtb = 0.001; % 1 kHz timebase for standard tail analysis
tailfiltlength = 7; % 7 ms boxcar filter of (1 kHz) tail data

if aq.settings.analysisversion >= 180119
    stimtypcol = 2; stimlatcol = 3; stimstcol = 4; stimedcol = 5;
else
    stimtypcol = 1; stimlatcol = 2; stimstcol = 3; stimedcol = 4;
end


% summary data variables

theta1 = []; % angl of 1st tail undulation
vel1 = []; % peak velocity first tail undulation
period1 = []; % duration of 1st tail undulation
trough1 = []; % amplitude of first trough
theta2 = [];
vel2 = [];
period2 = [];
pktime1 = []; % time of first peak (ms)

%% Find bouts

tt = aq.timebase;
tbi = 0:tailtb:(tt(end)-tailtb);
alltail = aq.cumtail;
alltail_median = nanmedian(alltail);
cumi = interp1(tt, alltail, tbi);
cumi = interpolate_NaNs(cumi);
cumif = filtfilt(ones(1, tailfiltlength)./tailfiltlength, 1, cumi);
cumif_m = cumif - alltail_median;
cumV = [0 diff(cumif)].*(1/tailtb); % degrees/sec

% tail vigor (rectified cumV, filtered 40ms)
tailvig_40 = filtfilt(ones(1,round(0.040./tailtb))./round(0.040./tailtb), ...
    1, abs(cumV-nanmean(cumV)));

% segmentangles full
seg_F = nan(size(aq.segmentangles_B, 1), numel(tbi));
for s = 1:size(aq.segmentangles_B, 1)
    thisseg = interp1(tt, aq.segmentangles_B(s, :), tbi);
    thisseg = interpolate_NaNs(thisseg);
    thisseg = filtfilt(ones(1, tailfiltlength)./tailfiltlength, 1, thisseg);
    seg_F(s, :) = thisseg;
end

% tail subsegments (median subtracted)
seg_R_m = nanmedian(sum(aq.segmentangles_B(2:3,:)));
seg_R = interp1(tt, sum(aq.segmentangles_B(2:3,:), 1) - seg_R_m, tbi);
seg_R(isnan(seg_R)) = 0;
seg_R = filtfilt(ones(1, tailfiltlength)./tailfiltlength, 1, seg_R);
seg_M_m = nanmedian(sum(aq.segmentangles_B(3:5,:)));
seg_M = interp1(tt, sum(aq.segmentangles_B(3:5,:), 1) - seg_M_m, tbi);
seg_M(isnan(seg_M)) = 0;
seg_M = filtfilt(ones(1, tailfiltlength)./tailfiltlength, 1, seg_M);
seg_C_m = nanmedian(sum(aq.segmentangles_B(5:7,:)));
seg_C = interp1(tt, sum(aq.segmentangles_B(5:7,:)) - seg_C_m, tbi);
seg_C(isnan(seg_C)) = 0;
seg_C = filtfilt(ones(1, tailfiltlength)./tailfiltlength, 1, seg_C);

tailvig_40_m = tailvig_40 - median(tailvig_40);
taildata = tailvig_40_m;

thrONhi = 1000;
thrONlo = 500;
thrOFFhi = 1000;
thrOFFlo = 500;
thrSEPfloor = 500;          % taildata must drop to this level before a new bout can be found
thrminbout = 35;            % min duration (ms) of the >ONhi period of the bout
thrmininterbout = 50;       % min interbout (ms)
boutz = findbouts(taildata, thrONhi, thrONlo, thrOFFhi, thrOFFlo, thrSEPfloor, thrminbout, thrmininterbout, tailtb);

% kill small bouts (min length 61 ms) and bouts with < 20 deg position
if ~isempty(boutz)
    boutz(:,3) = boutz(:,2) - boutz(:,1);
    boutz(boutz(:,3)<(0.061./tailtb), :) = nan;
    for b = 1:size(boutz,1)
        if ~isnan(boutz(b,1)) && ~any(abs(cumif_m(boutz(b,1):boutz(b,2))) > 20)
            boutz(b,:) = nan;
        end
    end
    boutz = unique(boutz, 'rows');
    boutz(isnan(boutz(:,1)) | isnan(boutz(:,2)), :) = [];
end

%%  Extract features

if ~isempty(boutz)
    boutz(:, 3:6) = 0;
    
    k = 1;
    for b = 1:size(boutz, 1)
        waitbar(b/size(boutz,1), h, 'Processing bouts...')
        boutossi = cumif_m(boutz(b,1):boutz(b,2));
        try
            boutdata = boutanal_v200(boutossi, tailtb, false);
            if ~isfield(boutdata, 'boutinfo')
                % abort this non-bout
                continue
            end
        catch
            continue
        end
        aq.btf(k).itst = findnearest(tbi(boutz(b,1)),tt,-1);
        aq.btf(k).ited = findnearest(tbi(boutz(b,2)),tt,-1);
        aq.btf(k).st = tbi(boutz(b,1));
        aq.btf(k).ed = tbi(boutz(b,2));
        aq.btf(k).seg = seg_F(:, boutz(b,1):boutz(b,2));
        
        % duration (s)
        aq.btf(k).dur = aq.btf(k).ed-aq.btf(k).st;
        
        % laterality
        aq.btf(k).intcum60ms = sum(cumif_m(boutz(b,1):boutz(b,1)+(0.060./tailtb)));
        latsign = sign(aq.btf(k).intcum60ms);
        
        % vigor
        if boutz(b,1)+(0.120./tailtb) < length(cumif_m)
            aq.btf(k).vig40 = sum(tailvig_40(boutz(b,1):boutz(b,1)+(0.120./tailtb)));
        else
            aq.btf(k).vig40 = sum(tailvig_40(boutz(b,1):end));
        end
        
        % morphology asymmetry index
        if boutz(b,1)+(0.120./tailtb) < length(cumif_m)
            ossismall = cumif_m(boutz(b,1):boutz(b,1)+(0.120./tailtb));
        else
            ossismall = cumif_m(boutz(b,1):end);
        end
        aq.btf(k).morphAI = ( sum(sign(ossismall)==latsign) - sum(sign(ossismall)~=latsign) )./numel(ossismall);
        aq.btf(k).morphAI2 = sum(ossismall)./numel(ossismall);
        
        % bendwise analysis
        boutinfo = boutdata.boutinfo;
        aq.btf(k).boutinfo = boutinfo;
        aq.btf(k).max_angl = boutdata.max_angl;
        aq.btf(k).mean_angl = nanmean(boutossi);
        aq.btf(k).max_vel = boutdata.max_vel;
        aq.btf(k).mean_vel = nanmean(abs(boutossi));
        aq.btf(k).pre_peak = boutdata.pre_peak;
        aq.btf(k).ifreq = boutdata.ifreq;
        aq.btf(k).mean_TBF = boutdata.mean_TBF;
        aq.btf(k).max_TBF = boutdata.max_TBF;
        
        % for boutinfo summary (to conform with gc structure)
        
        theta1(k,1) = boutinfo(1,4);
        vel1(k,1) = boutinfo(1,5);
        period1(k,1) = (boutinfo(1,2)-boutinfo(1,1))*tailtb;
        trough1(k,1) = boutinfo(1,6);
        pktime1(k,1) = boutinfo(1,3)*tailtb;
        if size(boutinfo,1)>1
            theta2(k,1) = boutinfo(2,4);
            vel2(k,1) = boutinfo(2,5);
            period2(k,1) = (boutinfo(2,2)-boutinfo(2,1))*tailtb;
        else
            theta2(k,1) = nan;
            vel2(k,1) = nan;
            period2(k,1) = nan;
        end
        
        % Visual stimulus present
        if ~isempty(aq.visstim)
            stimix = find(aq.btf(k).itst >= aq.visstim(:,stimstcol) & ...
                aq.btf(k).itst <= aq.visstim(:,stimedcol));
            if any(stimix)
                aq.btf(k).stim = aq.visstim(stimix(1),stimtypcol);
                aq.btf(k).stimdir = aq.visstim(stimix(1),stimlatcol);
            else
                aq.btf(k).stim = 0;
                aq.btf(k).stimdir = 0;
            end
        end
        
        % Convergence bout
        
        % Is bout start within any hunting routine?
        convix = find(aq.btf(k).itst >= aq.huntep(:,1) & ...
            aq.btf(k).itst <= aq.huntep(:,2));
        if ~isempty(convix)
            aq.btf(k).Convbout = convix;
        else
            % Is within 100ms of start of any routine?
            [cminstlat,convix] = min(abs(aq.huntep(:,1)-aq.btf(k).itst));
            if ~isempty(convix) && cminstlat <= aq.settings.camerarate*100/1000
                aq.btf(k).Convbout = convix;
            else
                % Is within 100ms of end of any routine
                [cminedlat,convix] = min(abs(aq.huntep(:,2)-aq.btf(k).itst));
                if ~isempty(convix) && cminedlat <= aq.settings.camerarate*100/1000
                    aq.btf(k).Convbout = convix;
                else
                    aq.btf(k).Convbout = 0;
                end
            end
        end
        
        % Looming bout
        loomix = find(aq.btf(k).itst >= aq.loomz(:,1) & ...
            aq.btf(k).itst <= aq.loomz(:,2));
        if ~isempty(loomix)
            aq.btf(k).Loombout = loomix(1);
        else
            aq.btf(k).Loombout = 0;
        end
        
        if ~isempty(aq.visstim)
            % Bout number within stimulus
            currstim = any(aq.btf(k).stim) | aq.btf(k).Convbout;
            if currstim && k ~= 1
                aq.btf(k).stimbout = aq.btf(k-1).stimbout+1;
            elseif currstim && k == 1
                aq.btf(k).stimbout = 1;
            else
                aq.btf(k).stimbout = 0;
            end
        end
        
        % fourier spectra
        for s = 1:4
            switch s
                case 1
                    boutossi = cumif_m(boutz(b,1):boutz(b,2)); % cumif section, median(whole epoch) subtracted
                case 2
                    boutossi = seg_R(boutz(b,1):boutz(b,2)); % cumif section, median(whole epoch) subtracted
                case 3
                    boutossi = seg_M(boutz(b,1):boutz(b,2)); % cumif section, median(whole epoch) subtracted
                case 4
                    boutossi = seg_C(boutz(b,1):boutz(b,2)); % cumif section, median(whole epoch) subtracted
            end
            boutossi = boutossi(1:2.*floor(numel(boutossi)./2)); % crop to even length
            % Fourier analysis
            fN = numel(boutossi);
            fFs = 1./tailtb;
            % fft of cumifs. mean-subtracted to avoid 0 Hz 'offset' components
            fZ = fft((boutossi-mean(boutossi))');
            fR = real(fZ);
            fI = imag(fZ);
            fP = (sqrt(fR.^2 + fI.^2)).^2;
            % find element corresponding to 100 Hz
            f100i = ceil(100*fN/fFs)+1;
            f100 = fP(1:f100i);
            % interpolate for consistency across bouts, upto 80Hz
            f100 = interp1(fFs.*[0:f100i-1]./fN, f100, [0:1:80]);
            aq.btf(k).sangl(s,:) = boutossi;
            aq.btf(k).fftspec(s,:) = f100;
            aq.btf(k).p5(s,1) = sum(f100(5:7)); % 4-(5)-6 Hz
            aq.btf(k).p30(s,1) = sum(f100(29:33)); % 28-(30)-32 Hz
        end
        % morphological asymmetries
        % fraction of cumulative tail position in Rostral, Middle and Caudal sections, at 1st peak:
        aq.btf(k).fcRMC1 = aq.btf(k).sangl(2:4,boutinfo(1,3))./aq.btf(k).sangl(1,boutinfo(1,3));
        
        %%  Fish max velocity, acceleraton and total displacement
        
        st = aq.btf(k).itst;
        fin = aq.btf(k).ited;
        aq.btf(k).max_fish_vel = max(abs(aq.fvelV(st:fin)));
        aq.btf(k).max_fish_acc = max(abs(diff(aq.fvelV(st:fin))));
        aq.btf(k).delta_ori = aq.fvelV(fin)-aq.fvelV(st);
        C = aq.datainterp(st:fin,2:3);
        mix = isnan(C(:,1));
        if any(mix)
            if sum(~mix) > size(C,1)/2
                C = interp1(find(~mix),C(~mix,:),1:size(C,1),'linear','extrap');
                aq.btf(k).t_disp = sum(sqrt((diff(C(:,1)).^2) + ...
                    (diff(C(:,2)).^2)),'omitnan')/aq.settings.CAMERApxpermm;
            else
                aq.btf(k).t_disp = nan;
            end
        else
            aq.btf(k).t_disp = sum(sqrt((diff(C(:,1)).^2) + ...
                (diff(C(:,2)).^2)),'omitnan')/aq.settings.CAMERApxpermm;
        end
        
        k = k+1;
    end
end

%% Bout classification

% ------------------ GET SUMMARY DATA -----------------------------

dur = cat(1, aq.btf.dur);
vig = cat(1, aq.btf.vig40);
lat = cat(1, aq.btf.intcum60ms);
morphAI = cat(1, aq.btf.morphAI);
morphAI2 = cat(1, aq.btf.morphAI2);
max_angl = cat(1, aq.btf.max_angl);
max_vel = cat(1, aq.btf.max_vel);
pre_peak = cat(1, aq.btf.pre_peak);

theta1 = theta1(:);
vel1 = vel1(:);
period1 = period1(:);
trough1 = trough1(:);
theta2 = theta2(:);
vel2 = vel2(:);
period2 = period2(:);

Convbout = cat(1, aq.btf.Convbout);
Loombout = cat(1, aq.btf.Loombout);

aq.btf_summary.dur = dur;
aq.btf_summary.vig = vig;
aq.btf_summary.lat = lat;
aq.btf_summary.morphAI = morphAI;
aq.btf_summary.morphAI2 = morphAI2;
aq.btf_summary.max_angl = max_angl;
aq.btf_summary.max_vel = max_vel;
aq.btf_summary.pre_peak = pre_peak;

% boutinfo stuff
aq.btf_summary.theta1 = theta1;
aq.btf_summary.vel1 = vel1;
aq.btf_summary.period1 = period1;
aq.btf_summary.trough1 = trough1;
aq.btf_summary.pktime1 = pktime1;
aq.btf_summary.theta2 = theta2;
aq.btf_summary.vel2 = vel2;
aq.btf_summary.period2 = period2;

% context markers
aq.btf_summary.Convbout = Convbout;
aq.btf_summary.Loombout = Loombout;


% ---------------- KINEMATIC THRESHOLD CLASSIFICATION ------------------
% ------------------NOT TESTED FOR AFAP DATA YET!!! --------------------

% v2:
aq.bouttype_classification = 'gcbft_manual/kinematic_v2';

tpratio = trough1./theta1;
FS = ( abs(theta1)<45 & morphAI<0.5 & abs(max_angl) < 100 );
RTL = ( theta1<-100 & tpratio<0 & morphAI<0.5 & abs(theta1-max_angl)<25 );
RTR = ( theta1>100 & tpratio<0 & morphAI<0.5 & abs(theta1-max_angl)<25 );
JL = ( theta1<-60 & tpratio>0.05 & morphAI>0.5 & abs(theta1-max_angl)<25 & abs(max_vel)<1.65e4 );
JR = ( theta1>60 & tpratio>0.05 & morphAI>0.5 & abs(theta1-max_angl)<25 & abs(max_vel)<1.65e4 );
OLBL = ( theta1>-110 & theta1<0 & morphAI>0.2 & max_vel<-2.4e4 ) | ( theta1<-180 & max_vel<0);
OLBR = ( theta1<110 & theta1>0 & morphAI>0.2 & max_vel>2.4e4 ) | ( theta1>180 & max_vel>0 );


% Stepwise classification
bouttype = nan(length(aq.btf),1);
bouttype(FS)    = 1;
bouttype(RTL)   = 2;
bouttype(RTR)   = 3;
bouttype(JL)    = 4;
bouttype(JR)    = 5;
bouttype(OLBL)  = 6;
bouttype(OLBR)  = 7;

aq.btf_summary.bouttype = bouttype;
aq.bouttype_num = [1:7];
aq.bouttype_name = {
    'FS';
    'RTL';
    'RTR';
    'JL';
    'JR';
    'OLBL';
    'OLBR';
    };


%% Check for inmiddle condition
for i = 1:size(aq.btf,2)
    if any(aq.inmiddle(aq.btf(i).itst:aq.btf(i).ited) == 0)
        aq.btf(i).inmiddle = 0;
    else
        aq.btf(i).inmiddle = 1;
    end
end

%%
waitbar(1,h,'Saving...');
save(fullfile(datadir,'aq'),'aq','-v7.3');
close(h)
end