% wmDrop_concatenateBehav_scanner.m
%
% combines all behavioral files into a single file per session, saves into
% wmDrop_behav
%
% TCS 5/4/2015

function wmDrop_concatenateBehav_scanner(subj)

root = load_root;%'/usr/local/serenceslab/tommy/wmDrop/';

if nargin < 1
    subj = {'AI81','AI82','AI83','AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
   
end

for ss = 1:length(subj)

    exp_name = 'wmDrop1dScanner';
    
    
    % for NC, original directory structure:
    %fn2s = sprintf('%s%s/%s_Behav/%s_%s_*.mat',root,subj{ss},subj{ss},subj{ss},exp_name);
    
    fn2s = sprintf('%swmDrop_behav/%s_Behav/%s_%s_*.mat',root,subj{ss},subj{ss},exp_name);
    
    fprintf('searching for %s...\n',fn2s);
    thisf = dir(fn2s);
    fprintf('%i files found\n',length(thisf));
    
    startidx = 1;
    
    for ff = 1:length(thisf)
        
        fn2load = sprintf('%swmDrop_behav/%s_Behav/%s',root,subj{ss},thisf(ff).name);
        fprintf('Loading %s...\n',fn2load);
        load(fn2load);
        
        if ff == 1
            
            % we want to save:
            % - all_conds
            % - all_resp
            % - all_respErr
            % - all_TRcoord
            % - all_TNcoord
            % - all_respCoord
            
            nblank = size(p.conditions,1)*length(thisf);
            
            all_conds = nan(nblank,size(p.conditions,2));
            all_resp = nan(nblank,1);
            all_respErr = nan(nblank,1);
            all_TRcoord = nan(nblank,2);
            all_TNcoord = nan(nblank,2);
            all_respCoord = nan(nblank,1);
            
        end
        
        thisidx = startidx:(startidx+size(p.conditions,1)-1);
        
        all_conds(thisidx,:) = p.conditions;
        all_resp(thisidx) = p.wmResp/p.ppd;   % already locked to center
        all_respErr(thisidx) = abs(p.wmRecallDistanceDeg.');
        % need to recenter these
        all_TRcoord(thisidx,:) = (p.TRcoord - repmat([1024 768]/2,size(p.conditions,1),1))/p.ppd;
        all_TNcoord(thisidx,:) = (p.TNcoord - repmat([1024 768]/2,size(p.conditions,1),1))/p.ppd;
        all_respCoord(thisidx) = p.respCoord;
        
        startidx = thisidx(end)+1;
        clear p;
        
    end
    
    
    fnsave = sprintf('%swmDrop_behav/%s_%s_behav.mat',root,subj{ss},exp_name);
    fprintf('saving combined behavior to %s...\n\n',fnsave);
    save(fnsave,'all_conds','all_resp','all_respErr','all_TRcoord','all_TNcoord','all_respCoord');
end



return