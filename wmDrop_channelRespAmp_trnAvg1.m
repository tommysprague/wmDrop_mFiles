function wmDrop_channelRespAmp_trnAvg1(subj,VOIs,hex_size)
% TCS 2/6/15 - adapted form wmDrop_channelRespAmp1.m
% - this one uses wmDrop_extractSignal_avg.m, which averages over input TRs
%   (here, [3 4], like wm_delay (CB)



% structure I want:
%
% wmDrop_channelRespAmp1
% - within session, loads training and testing data, sends this to
%   computeChanResp fcn
% - runs extractSignal_tst/trn on respective mr files
%   - use trial-level betas for training
%   - use full event-related timecourse for testing
% - 


% TODO: remove first components which make bilateral MR's, put that in
% another file - this function should just be about loading data (bilateral
% mr) and computing over it. 

%
% updated 2/2/2016 by TCS to look in wmDrop_mrStructs instead of subj
% directory for mr files
%

if nargin < 1
    
    subj = {'AI81','AI82','AI83','AP81','AP82','AP83','AR81','AR82','AR83','AS81','AS82','AS83','AL81','AL82','AL83','BC81','BC82','BC83'};
    
end


if nargin < 2
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop','iPCS','SMA','DLPFC'};

end

if nargin < 3
    % the size of the hexagonal grid - we used a grid of 37 channels (7)
    % for Sprague, Ester & Serences, 2016 - but feel free to try other
    % values (5 will reduce size of grid; 9 will increase to 61, 11 will increase to 91). 
    hex_size = 7;
end

trnPrefix = 'wmMapHex';


filt_scale = 1.1; % FWHM: spacing ratio need to multiply by 1/rad2fwhm(1) to convert to cos7 size constant


root = load_root;

if ~exist([root 'wmDrop_chanResp'],'dir')
    fprintf('creating directory: %s\n',[root 'wmDrop_chanResp']);
    mkdir([root 'wmDrop_chanResp']);
end

for ss = 1:size(subj,2)   % loop over subjects


    
    for vv=1:length(VOIs)

        %% load trialData
        fn = sprintf('%swmDrop_trialData/%s_%s_wmDrop_trialData.mat',root,subj{ss},VOIs{vv});
        fprintf('loading imaging data from: %s\n',fn);
        load(fn); clear fn;
        
        % tr_coord, tn_coord are in CARTESIAN coordinates (+,+ for quadrant
        % 1), ang_offset is in radians, was added to theta of stimLocsX/Y in
        % behavioral script BEFORE flipping Y, so when working in screen
        % coords, -ang_offset was added (so use +ang_offset to fix)
        
        %% generate design matrix for training data
        %
        % can do this just once, before loading each ROI's data for
        % each subj (design matrix is the same for all ROIs)
        
        if vv == 1
            
            fprintf('generating design matrix\n');
            
            % load this subject's stimulus parameter file
            fn = sprintf('%s%s_%s_params.mat',root,subj{ss},trnPrefix);
            fprintf('loading stimulus parameters from : %s\n',fn);
            load(fn); clear fn;
            
            % stimulus was part of a make_hexagon(7) grid, made up of
            % equilateral triangles with wmDropHex_params.stepSize sides, so
            % the maximum distance from fixation the stimulus reached was:
            %   max_dist = 3 * wmMapHex_params.stepSize + wmMapHex_params.radDeg
            
            [rfX, rfY] = make_hex(hex_size); % 11 = 91 channels - normalized grid, max value (x,y) = 1/-1, 9 = 61 channels
            max_dist = 3 * wmMapHex_params.stepSize + wmMapHex_params.radDeg;
            rfX = rfX*max_dist*(2/sqrt(3));  rfY = rfY*max_dist*(2/sqrt(3));
            
            basis_sep = rfX(2)-rfX(1); % how far, in DVA, the basis functions are spaced apart
            rfSize = basis_sep * filt_scale;  % size of basis functions, in FWHM units
            
            
            
            trnX = wmDrop_makeX1(stimLocs,wmMapHex_params.radDeg,rfX,rfY,rfSize,501);
            trnX = trnX/max(trnX(:));
            
            
            clear wmMapHex_params max_dist;
            
        end
        

        
        
        
        %% compute channel weights using only training set (wmDropHex)
        
        fprintf('computing channel weights\n');
        w = trnX\b_trn;
        
        
        %% use design matrix to compute channel responses using testing data
        
        fprintf('computing channel responses\n');
        
        chan_resp = inv(w*w')*w*tst';
 
        

        chan_resp = chan_resp - nanmean(chan_resp(:)) + nanmean(tst(:));



        
        
        
        
        
        %% save channel responses, etc
        
        fn2s = sprintf('%swmDrop_chanResp/%s_%s_hex%i_channelResp_trnAvg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('saving to %s...\n\n',fn2s);
        save(fn2s,'chan_resp','w','rfX','rfY','rfSize','conds','tst_scans','tpts','tr_coord','tn_coord','targ_idx','ang_offset');
        clear fn2s w conds tst_scans tpts tr_coord tn_coord targ_idx ang_offset;
        
    end % end VOI loop
    clear rfX rfY rfSize ;

end % end subj loop

end % end fcn




