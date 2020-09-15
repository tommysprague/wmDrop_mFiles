% wmDrop_allAnalyses1.m
%
% TCS 10/27/14 - this is where we'll put the full workflow of analyses.
% plotting will be a separate function, but it will require that all these
% are run


%addpath(genpath([load_root 'mFiles']));


subj = {'AI81','AI82','AI83','AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop','iPCS','DLPFC','SMA'};

hex_size = 7; % only free parameter for IEM analysis - used to generate basis set
delay_tpts = {[3 4],[7 8]}; % TRs after trial start (0 = volume nearest targets) for each delay (see Fig. S5 for comparison of different choices)

all_tpts = 0:9; % for time-series analyses

%reconstruction_resolution = 51;
reconstruction_resolution = 101;

%% behavioral analyses

% first concatenate behavioral files
wmDrop_concatenateBehav_scanner(subj); % for behavioral analyses



wmDrop_combineBehav1(subj);  % for analyses which split by median recall error


% now run the plotting script (generates Fig. 1C and Fig. S1)
% statistics will spit out to command line

wmDrop_plotBehav2(subj);

%% plot HRFs (S2)
%
% load wmDrop_trialData files and plot mean HRFs across all trials, voxels
% within an ROI (Fig. S2); compare mean delay period activation between
% delay periods for each condition

wmDrop_plotHRFs_ERA1(subj,VOIs,delay_tpts);


%% compute channel responses

wmDrop_channelRespAmp_trnAvg1(subj,VOIs,hex_size);



%% compute reconstructions
% NOTE: this step creates >1000 GB of files using parameters as used for
% Sprague, Ester & Serences, 2016. To allow for reasonable computations on
% typical computing setups, we're reducing the reconstruction resolution to
% 51 x 51 pixels (from 101 x 101). This *will* produce slightly different
% results than plotted in our report. Should your computing systems be
% adequate, or should you adjust the plotting scripts accordingly, you can
% adjust the "reconstruction_resolution" parameter above. (in case you
% still receive a "Requested array exceeds maximum array size"
% warning/error, try a smaller resolution)
%
% for 51 x 51: ~14 GB for all ROIs, subj

% this computes reconstructions as in Figs 3-4 (averaged within target
% separation conditions, Fig. 1B), entirely for visualization
%
% (each command takes ~15-20 mins on new iMac)
wmDrop_computeReconstructions_trnAvg_posCoreg1(subj,VOIs,hex_size,reconstruction_resolution);


% for all 'representational fidelity' analyses (Figs 5-6), we use reconstructions
% computed along an anulus from 2.9-4.1 dva from fixation. I call these
% "rotated" reconstructions. For all analyses presented in the text, we
% averaged over all eccentricities with a polar angle to generate 1d
% reconstructions over polar angle. Generate those reconstructions as
% below. This requires ~12 GB storage for all ROIs, subj
wmDrop_computeReconstructions_rotateCoreg1(subj,VOIs,hex_size)

% we also need to resample across trials for plotting 1-d reconstructions
% with CIs (Fig. 6)
% TCS 9/12/2020 - this correctly resamples
wmDrop_resampleReconstructions_rotateCoreg1(subj,VOIs,delay_tpts,hex_size); % TCS - added delay_tpts arg


% this computes reconstructions used for surface fitting (Figs. 7-8, S8)
wmDrop_computeReconstructions_trnAvg_exactCoreg1(subj,VOIs,hex_size,reconstruction_resolution);




%% plot reconstructions (Figs. 3-4)

% Fig. 3 (this will plot all ROIs, only {'V1','V3A','IPS0','IPS2','sPCS','SuperWMDrop'}
% used for publication)
wmDrop_plotReconstructions_posCoreg_thruTime1(subj,VOIs,0:9,hex_size);


% Fig. 4: (this will plot all ROIs; only
% {'V1','V3A','IPS0','IPS2','sPCS','SuperWMDrop'} for publication
wmDrop_plotReconstructions_trnAvg_posCoreg1(subj,VOIs,delay_tpts,hex_size);


% Fig. 5:
wmDrop_plotReconstructions_rotateCoreg1_thruTime(subj,VOIs,all_tpts,hex_size);

% Fig. 6:
% note - we use the # of subj, not the subj ID's
wmDrop_plotReconstructions_rotateCoreg1_resample(6,VOIs,delay_tpts,hex_size);


%% plot representational fidelity (Figs. 5-6, S5-S7)
% NOTE: to achieve exact plotting results as shown in SES 2016, do not use
% the same subject input order as above. Calling without arguments will use
% the matching subject order with corresponding randomization seed. Feel
% free to try different seeds, etc, but only the included seed will produce
% identical results. (this is because randomization happens on a list of
% trials loaded in order; if the trials are loaded in a different order, a
% different set of trials will be grouped during each randomization
% iteration).
% TCS 9/12/2020: seems likely this is only using n = 3 for delay-period stats
wmDrop_vectorMean_rotateCoreg1_CR;


%% plot/compare target response amplitude
% NOTE: to achieve the exact plotting results as reported, must use subject
% order listed in this function (call w/ no input arguments). 

% plot Fig. S4 (and an extra figure showing each target's response)
wmDrop_plotTargetResp_resample1;


%% compute surface fits (Figs. 7-8, S8)
% As above, to achieve the exact plotting results, run without arguments so
% that original subject order is used in conjunction with random seed. We
% also include the best-fit and resampled reconstructions as computed by
% us if you prefer to forego this step (~20 GB or so of files). 

% for Fig. 7 (all trials within each condition)
wmDrop_fitReconstructions_exactCoreg1_CR;

% for Fig. 8, S8 (split by median within condition; session)
wmDrop_fitRecons_exactCoreg_byRecallError1_CR;

%% plot surface fits

% for Fig. 7 (note we use n subjects)
wmDrop_plotFits_trnAvg_exactCoreg1_resample(6,VOIs,delay_tpts,hex_size)

% for Fig. 8 ({'SuperWMDrop'}) and Fig. S8
% {V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'};
% Fig. 8:
wmDrop_plotFits_trnAvg_exactCoreg_byRecallError1_resamp(6,{'SuperWMDrop'},delay_tpts,hex_size);

% Fig. S8:
wmDrop_plotFits_trnAvg_exactCoreg_byRecallError1_resamp(6,{'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'},delay_tpts,hex_size);


