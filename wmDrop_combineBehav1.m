% wmDrop_combineBehav1.m
%
% combines behavioral data (recall error) across all runs for each subj
%
%
% TCS 1/12/2015

function wmDrop_combineBehav1(subj)

root = load_root;%'/usr/local/serenceslab/tommy/wmDrop/';

if ~exist([root 'wmDrop_combinedBehav'],'dir');
    fprintf('Creating directory: %s\n',[root 'wmDrop_combinedBehav']);
    mkdir([root 'wmDrop_combinedBehav']);
end

if nargin < 1
    subj = {'AI81','AI82','AI83','AL81','AL82','AL83','AP81','AP82','AP83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
end

for ss = 1:length(subj)
    
    recall_error = [];
    conditions = [];
    % on NC:
    %fs = sprintf('%s%s/%s_Behav/%s_wmDrop1dScanner*.mat',root,subj{ss},subj{ss},subj{ss});
    
    % for OSF:
    fs = sprintf('%swmDrop_behav/%s_Behav/%s_wmDrop1dScanner*.mat',root,subj{ss},subj{ss});
    
    fprintf('searching for %s..\n',fs);
    allf = dir(fs);
    
    for ff = 1:length(allf)
        
        thisf = sprintf('%swmDrop_behav/%s_Behav/%s',root,subj{ss},allf(ff).name);
        fprintf('loading %s..\n',thisf);
        load(thisf);
        conditions = [conditions;p.conditions];
        recall_error = [recall_error;p.wmRecallDistanceDeg'];
        
        clear p;
    end
    
    fn2s = sprintf('%swmDrop_combinedBehav/%s_wmDrop_combinedBehav.mat',root,subj{ss});
    fprintf('saving %s data to %s...\n',subj{ss},fn2s);
    save(fn2s,'recall_error','conditions','allf'); clear recall_error allf;
    
    
end


return