function wmDrop_computeReconstructions_rotateCoreg1(subj, VOIs,hex_size)
% computes reconstructions along an anulus around fixation spanning 2.9-4.1
% dva. Used for all "representational fidelity" analyses (Fig. 2C, Figs
% 5-6) and 1D reconstructions (Figs 5-6).

if nargin < 1
    subj = {'AI81','AI82','AI83','AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
end

if nargin < 2 
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop','iPCS','SMA','DLPFC'}; 

end

if nargin < 3
    hex_size = 7;
end

root = load_root;


chan_str = 'trnAvg1'; % 

ecc_to_align = 3.5; % align to this point
ecc_rad = 0.6;  % distance from ecc_to_align within whcih targets could ever appear

res_r = 12; % in x, y  % or 101
res_t = 220;   % each step shoudl be close to 1.6372 deg

gridt = linspace(-pi+(pi/res_t),pi,res_t);
gridr = linspace(ecc_to_align-ecc_rad,ecc_to_align+ecc_rad,res_r);
[gridt,gridr] = meshgrid(gridt,gridr);

[gridx,gridy] = pol2cart(gridt,gridr);

gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);



for ss = 1:length(subj)
    
    for vv = 1:length(VOIs)
        
        % load data from /matData/
        % subj_ROI_channelRespAmp.mat
        chan_fn = sprintf('%swmDrop_chanResp/%s_%s_hex%i_channelResp_%s.mat',root,subj{ss},VOIs{vv},hex_size,chan_str);
        
        fprintf('loading %s...\n',chan_fn);
        chan = load(chan_fn);
        
        [rfT, rfR] = cart2pol(chan.rfX,chan.rfY);
        
        % cartesian: x/y
        t_coord_cart{1} = chan.tr_coord;
        t_coord_cart{2} = chan.tn_coord;
        
        % theta/r
        [t,r] = cart2pol(chan.tr_coord(:,1),chan.tr_coord(:,2));
        t_coord_pol{1} = [t r]; clear t r;
        [t,r] = cart2pol(chan.tn_coord(:,1),chan.tn_coord(:,2));
        t_coord_pol{2} = [t r]; clear t r;
        
        % create a trial_num index so that we don't need to keep generating
        % basis sets
        tu = unique(chan.tpts);
        n_trials = size(chan.conds,1)/length(tu);
        
        trial_num = repmat((1:n_trials)',length(tu),1);
        
        recons_vec = cell(2,1);
        
        for ii = 1:2   % TR, TN
            recons_vec{ii} = nan(size(chan.chan_resp,2),numel(gridx));
            for tt = 1:n_trials
                
                % adjust rfX, rfY to correct polar angle, ecc = 3.5
                this_rfT = rfT-t_coord_pol{ii}(tt,1);
                
                this_rfR = rfR;
                
                [this_rfX, this_rfY] = pol2cart(this_rfT,this_rfR);

                
                
                % build basis
                basis_set = build_basis_pts(this_rfX,this_rfY,chan.rfSize/rad2fwhm(1),gridx,gridy);
                
                % get all timepoints from that trial
                thisidx = trial_num==tt;
                thisdata = chan.chan_resp(:,thisidx);
                
                % reconstruct 
                recons_vec{ii}(thisidx,:) = thisdata'*basis_set';
                
            end
            
        end
        
        conds = chan.conds;
        tpts = chan.tpts;
        
        % save file
        fn2s = sprintf('%swmDrop_recons/%s_%s_hex%i_%s_rotate_coreg1.mat',root,subj{ss},VOIs{vv},hex_size,chan_str);
        fprintf('saving to %s...\n',fn2s);
        save(fn2s,'recons_vec','conds','tpts','t_coord_cart','chan_fn','res_r','res_t','gridx','gridy','ecc_rad','ecc_to_align');
        
        clear chan;
        
    end
    
end


return