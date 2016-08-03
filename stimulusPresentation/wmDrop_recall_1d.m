

function wmDrop_recall_1d
% adapted from wmDrop_recall.m
%
% now uses a recall task - thin bar appears on screen at end of each trial,
% must adjust it to match the given coordinate of target as best as
% possible

% - remember 1, remember 2, drop
%
%
% 8 s delay - delay 1
% 8 s delay - delay 2
%


%% get user input about sub name etc...

% file/directory operations
p.root = pwd;

if ~isdir([p.root '/data/'])
    mkdir([p.root '/data/']);
end

db_probe = 0;
draw_target_frames = 0;
p.show_eye_pos = 0;

% for behavioral testing room, 60 Hz, approx 3.2 cd/m^2 luminance (against
% [128 128 128] - 9.68 cd/m^2, for approx 50% contrast (michelson..)
% mygreen = [0 82 0];
% myred   = [100 0 0];
% myblue  = [30 30 255];
% mypurple= [135 0 135];
% myorange = [118 47 0];
% mygray  = [65 65 65];


% these are all easily discriminable, so let's just use them
% TODO: measure luminance of each both at scanner and in testing room
myred = [255 0 0];
myblue = [50 50 255];
mypurple = [135 0 135];

mywhite = [255 255 255];
myblack = [0 0 0];
mygray = [65 65 65]; % fix this one
mygray2 = [110 110 110];

warning('off','MATLAB:dispatcher:InexactMatch');
%Screen('Preference', 'SkipSyncTests', 1)
%get subject info
prompt = {'Subject Name','Session number', 'Run Number', 'Run Subsample', 'Random seed', 'display (1 = fMRI, 2 = laptop, 3 = behavioral room)','Eye tracking','Max angular offset (deg)'};
%grab a random number and seed the generator (include this as a gui field
%in case we want to repeat an exact sequence)
s = round(sum(100*clock));
%put in some default answers
defAns = {'XX','3','X','X',num2str(s),'1','0','15'};
box = inputdlg(prompt,'Enter Subject Information...', 1, defAns);


if length(box)==length(defAns)
    p.subName=char(box{1});
    p.sessionNum=str2num(box{2});
    p.runNum = eval(box{3});
    p.runSub=eval(box{4});
    p.rndSeed=str2num(box{5});
    p.display=str2num(box{6});   
    p.eyeTracking = str2num(box{7});    
    p.maxAngOffset = str2num(box{8});
    rand('state',p.rndSeed);  %actually seed the random number generator
else
    return
end

%debug_stim_loc = 0;

% note: fMRI overrides fullscr param...
if p.display == 1
    p.fMRI = 1;
    p.fullScr = 1;
    
    p.exptName = 'wmDrop1dScanner';

else
    p.fMRI = 0;
    p.fullScr = 1;

    p.exptName = 'wmDrop1dScannerBehav'; % because this will also be run on first session
end

%p.staircaseStepSize = 0.1;

%p.keys = 66; %keys that correspond to elements in p.dir respectively (B)


%% --------------------begin user define parameters----------------------------%

% Screen (parameter) setup
p.windowed = ~p.fullScr;                     % if 1 then the display will be in a small window, use this for debugging.

% monitor stuff
p.refreshRate = 60;                 % refresh rate is normally 100 but change this to 60 when on laptop!

if p.display == 1  % fMRI
    % 3TW
    %p.vDistCM = 308; % CM
    %p.screenWidthCM = 90; % CM
    
    % 3TE
    p.vDistCM = 375;
    p.screenWidthCM = 120;
    
elseif p.display == 2 % laptop
    p.vDistCM = 60;  % rough guesses (from map2d_3)
    p.screenWidthCM = 32;
    p.percentNull = 0;
else
    % measured 7/30/12 - TCS - in eyetracker room
    p.vDistCM = 62;
    p.screenWidthCM = 51;
    p.percentNull = 0;
end


p.nLoc = 6; % stim will appear at random location on nLoc by nLoc grid (should be even)
p.nSepConds = 3;  % same quadrant, neighboring quadrants across or within hemifield, diagnoally-separated quads
% 1 = 60 polar deg
% 2 = 120 polar deg
% 3 = 180 polar deg
% also 50% p(L/R) shift from target position

p.nMemConds = 3;  % remember 1, 2, or 2 then drop items

p.nt = 18;  % nTrials/3, corresponds to 3 subruns

% TODO: below
%stimulus geometry (in degrees, note that these vars have the phrase 'Deg'
%in them, used by pix2deg below)
p.radDeg = 3.5; % center of point clouds  % TODO
p.jitterRadDeg = 0.6; % dispersion of point clouds (radius) % TODO
% if 8 locations, center of each location is, at closest, 2.1048 (w/ 5.5
% deg radius) away from hemifield boundary - so jitterRad << this is good

p.stepSizeDeg = .06; % pixels per frame of holding button down

p.wmDistanceLimitDeg = p.jitterRadDeg * 1.5;

% font stuff
p.fontSize = 24;
p.fontName = 'ARIAL';
p.textColor = [100, 100, 100];

%Response keys
if ismac
    p.escape = 41;
else
    p.escape = 27;
end
if p.fMRI
    p.keys = [KbName('b'),KbName('y')]; % 1 and 2 keys
else
    p.keys = [KbName('left'), KbName('right')];  % left same, right different
end
p.space = KbName('space');
p.start = KbName('t');

p.backColor     = [128, 128, 128];      % background color

%fixation point properties
p.fixColor      = [180 180 180];%[155, 155, 155];
% p.fixSizeDeg    = .25 * .4523;                  % size of dot
p.fixSizeDeg    = 0.2;   %.1131;
p.appSizeDeg    = p.fixSizeDeg * 4;     % size of apperture surrounding dot
p.fixWinDeg     = p.appSizeDeg * 2;   % TODO: CHECK THIS

p.crossHairSizeDeg = .5;

% target/probe properties (WM) ( both 0.5 - now 0.5 * .4523)
p.targetSizeDeg = 0.15;
p.probeSizeDeg = 0.15;
p.wmTargetColors = [myred;myblue];
p.wmCueColors = [myred; myblue; mypurple; mygray];
p.wmProbeColors = [mywhite;myblack];

p.wmRecallColor = mygray2;

% target properties (fix/checks)

% trial setup

p.repetitions = 1;
if p.windowed % DEBUG MODE
    p.repetitions = 1;
end
% this is the TOTAL number of conditions...
p.nTrials = p.repetitions * p.nSepConds * p.nMemConds * p.nLoc;    % length p.cond, for now, should always be 1

% stimulus timing (in seconds - convert to vid frames, if necessary, later)
p.wmDelay1Dur    = 8 * ones(p.nt,1); % keep this constant
p.wmDelay2Dur    = 8 * ones(p.nt,1); % keep this constant (in scanner)
p.wmTargetDur   = 0.5 * ones(p.nt,1); % s
%p.wmProbeDur    = 2.0 * ones(p.nt,1); % s (made this longer, very hard task)

%p.delayAfterResponse = 0.5;
%p.flickerFreq   = 6; % Hz

p.responseWindow = 3.0 * ones(p.nt,1);  % maximum time to make response
p.nTargs = 2;

% time between stimulus presentations, constant
p.minITI = 3.0; % seconds
p.maxITI = 6.0; % seconds - linearly space between these (below)
p.ITI           = linspace(p.minITI,p.maxITI,p.nt)';  % should we randomize this behaviorally?
p.ITI           = p.ITI(randperm(length(p.ITI)));  % shuffle ITIs
p.trialDur      = p.wmTargetDur + p.wmDelay1Dur + p.wmDelay2Dur + p.responseWindow + p.ITI;  % note - removed p.responseWindow from here - that should be part of ITI...
if p.display == 1
    p.nTRsWait      = 0;  % no dummy scans
    p.passiveDur    = 10; % in sec, how long to wait at end of a block (and beginning?)
else
    p.nTRsWait = 0;
    p.passiveDur = 1;
end

p.startWait     = 2;  % after 5 dummy TRs waiting
p.TR            = 2.25;
p.expDur        = sum(p.trialDur) + p.passiveDur + p.startWait + p.nTRsWait*p.TR;

disp(sprintf('SCAN DURATION: %i',p.expDur));
disp(sprintf('# TRs: %i',ceil(p.expDur/p.TR)));  % XXX TRs

ListenChar(2);

%% --------------------Screen properties----------------------------%

p.LUT = 0:255;
% correct all the colors
p.fixColor = p.LUT(p.fixColor)';
%p.targetFixColor = p.LUT(round(p.targetFixColor))';
p.textColor = p.LUT(p.textColor)';
p.backColor = p.LUT(p.backColor)';

%Start setting up the display
AssertOpenGL; % bail if current version of PTB does not use

% Open a screen
Screen('Preference','VBLTimestampingMode',-1);  % for the moment, must disable high-precision timer on Win apps

% figure out how many screens we have and pick the last one in the list
s=max(Screen('Screens'));
p.black = BlackIndex(s);
p.white = WhiteIndex(s);

if p.windowed
    Screen('Preference', 'SkipSyncTests', 1);
    [w p.sRect]=Screen('OpenWindow', s, p.backColor, [50,50,800,600]);
    HideCursor;
else
    [w p.sRect]=Screen('OpenWindow', s, p.backColor);
    HideCursor;
end

p.center = [p.sRect(3)/2 p.sRect(4)/2];

% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% test the refesh properties of the display
p.fps=Screen('FrameRate',w);          % frames per second
p.ifi=Screen('GetFlipInterval', w);   % inter-frame-time
if p.fps==0                           % if fps does not register, then set the fps based on ifi
    p.fps=1/p.ifi;
end
%p.flickerFrames = 1/p.flickerFreq*p.fps;

% make sure the refreshrate is ok
if abs(p.fps-p.refreshRate)>5
    Screen('CloseAll');
    disp('CHANGE YOUR REFRESH RATE')
    ListenChar(0);
    clear all;
    return;
end


% if running the real experiment (not debugging), then hide cursor and set
% priority high
if ~p.windowed
    HideCursor; % Hide the mouse cursor
    % set the priority up way high to discourage interruptions
    Priority(MaxPriority(w));
end

% convert from degrees to pixel units
p = deg2pix(p);

% now we have p.usedScreenSizePix, p.radPix

% compute and store the center of the screen: p.sRect contains the upper
% left coordinates (x,y) and the lower right coordinates (x,y)
center = [(p.sRect(3) - p.sRect(1))/2, (p.sRect(4) - p.sRect(2))/2];

xcenter = center(1); ycenter = center(2);

p.xOffsetPix = 0;
%left = [center(1)-p.xOffsetPix, center(2)];

Screen('TextSize', w, p.fontSize);
Screen('TextStyle', w, 1);
Screen('TextFont',w, p.fontName);
Screen('TextColor', w, p.textColor);

% file/directory operations
%p.root = 'C:\Users\RA\Dropbox\Documents\UCSD\Serences\spatialWM';

if ~isdir([p.root '/data/'])
    mkdir([p.root '/data/']);
end
%p.subNameAndDate = [p.subName '_' datestr(now,30)];


    %fName=[p.root, '/data/', p.exptName,'_', p.subNameAndDate, '_sess', num2str(p.sessionNum), '_blkNum', num2str(b), '.mat'];
    fName = sprintf('%s/data/%s_%s_sess%02.f_run%02.f_subsamp%02.f.mat',p.root, p.subName,p.exptName,p.sessionNum,p.runNum,p.runSub);
    if exist(fName,'file')
        Screen('CloseAll');
        msgbox('File name already exists, please specify another', 'modal')
        ListenChar(0);
        return;
    end
    p.dateTime = datestr(now,30);
    
    % setup eyetracker before every block - this does calibration!
    if p.eyeTracking  == 1 % SMI
        
        addpath('Z:/ETTools/');
        
        try
            Screen('CloseAll');
            [pSampleData pEventData] = setup_SMI(p.backColor(1));
            disp('setup finished');
            
            %  expt/subj/expt_subj_sess_block_date.idf
            p.et_fname = sprintf('%s/%s/%s_%s_sess%02.f_run%01.f_subrun_%01.f_%s.idf',p.exptName,p.subName, p.subName,p.exptName,p.sessionNum,p.runNum,p.runSub, datestr(now,30));
            
            % open Screen again...
            if p.windowed
                Screen('Preference', 'SkipSyncTests', 1);
                [w]=Screen('OpenWindow', s, p.backColor, [50,50,800,600]);
                HideCursor;
            else
                [w]=Screen('OpenWindow', s, p.backColor);
                HideCursor;
            end
            
            verify_calib_SMI(w,pSampleData,pEventData)
            Screen('TextSize', w, p.fontSize);
            Screen('TextStyle', w, 1);
            Screen('TextFont',w, p.fontName);
            Screen('TextColor', w, p.textColor);

        catch
            disp('quit during setup...');
            p.eyeTracking = 0;
        end
        
    elseif p.eyeTracking == 2 % ASL
        
        disp('setting up eyetracking...');
        try
            write_parallel(0); % initial communication w/ eyetracker?
        catch thisError
            disp('EyeTracking error...continuing with p.eyeTracking = 0;');
            disp(thisError);
            p.eyeTracking = 0;
            p.eyeTrackingStopped = 'at setup';
        end
    end
    
    p.allAngOffset = deg2rad(linspace(-p.maxAngOffset,p.maxAngOffset,p.nt))';  % separately randomized from conditions, this ensures full coverage of visual field
    p.allAngOffset = p.allAngOffset(randperm(p.nt));
    
    % now do stimLocsT/R and stimLocsX/Y for each trial (n_trials x n_loc)
    stimLocsT = repmat(linspace(2*pi/p.nLoc,2*pi,p.nLoc)-2*pi/p.nLoc,p.nt,1) + repmat( p.allAngOffset, 1, p.nLoc );% deg2rad(p.angOffset);
    stimLocsR = p.radPix*ones(p.nt,p.nLoc);
    
    stimLocsX = center(1) + stimLocsR.*cos(stimLocsT);
    stimLocsY = center(2) - stimLocsR.*sin(stimLocsT);

    p.conditions = nan(p.nTrials,4);  % nulls are nans at the end...
    cnt = 1;
    %locCnt = 0;
    for k = 1:p.nSepConds
        for j = 1:p.nMemConds       
            for i = 1:p.nLoc
                for ii = 1:p.repetitions
                    p.conditions(cnt,[1 2 3 ]) = [i j k ];  % location, memory condition (0, 1 or 2), separation condition (1-3), probe distance
                    cnt = cnt+1;
                end
            end
        end
    end
    clear k j i ii cnt;
    
    
    p.conditions(:,4) = -1+2*round(rand(size(p.conditions,1),1)); % 1 = CCW, -1 = CCW
    
    
    % shuffle trial order
   % p.nt = 24;  (defined this above - where ITIs are defined)
    if p.runSub==1
        p.rndInd = randperm(p.nTrials);
        p.trInd = p.rndInd(1:p.nt);
    else
        
        oldfName = sprintf('%s/data/%s_%s_sess%02.f_run%02.f_subsamp%02.f.mat',p.root, p.subName,p.exptName,p.sessionNum,p.runNum,1);
        tmp = load(oldfName);
        p.trInd = tmp.p.rndInd(p.nt*(p.runSub-1)+1:p.nt*(p.runSub));
        p.rndInd = tmp.p.rndInd; % save al of them..
        clear tmp;
    end
    % below accomplishes only running nt trials
    p.conditions = p.conditions(p.trInd,:);
    %p.rt =          nan(1, p.nTrials);       % store the rt on each trial
    %p.resp =        zeros(p.nTrials, 2);     % store the response
    
    TNloc = mod(p.conditions(:,1) + p.conditions(:,4) .* p.conditions(:,3),p.nLoc);
    
    TNloc(TNloc==0) = p.nLoc;
    
    % 180 deg is symmetric, don't need +/- 1 here
    p.conditions(p.conditions(:,3)==3,5)=0;
    
    p.TRloc = p.conditions(:,1);
    p.TNloc = TNloc;
    
    
    p.respCoord = round(rand(p.nt,1))+1;  % 1 for left/right (x, vertical line), 2 for up/down (y, horiz line)

    
    p.trialStart =  nan(p.nt,1);    % equal to targetStart
    p.trialEnd =  nan(p.nt,1);
    p.stimStart =  nan(p.nt,1);
    p.stimEnd =  nan(p.nt,1);
    p.delay1Start = nan(p.nt,1);
    p.delay2Start = nan(p.nt,1);
    
    p.wmResp = nan(p.nt,1);
    p.responseTime = nan(p.nt,1);
    
    p.respCoordTrace = nan(p.nt,p.fps*p.responseWindow(1));
    
    % TODO: if in scanner
    textWait = 'Waiting for start of scanner';
    
    tCenterMo = [center(1)-RectWidth(Screen('TextBounds', w, textWait))/2 center(2)/2];
    
    Screen('DrawText', w, textWait, tCenterMo(1), tCenterMo(2), p.textColor);
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
    %after all initialization is done, sit and wait for scanner synch (or
    %space bar)
%     resp=0;
%     
%     while 1
%         [resp, timeStamp] = checkForResp([p.space,p.start],p.escape);
%         if resp==p.space || resp==p.start || resp==p.escape
%             break;
%         end
%         if resp==-1; ListenChar(0); return; end;
%     end
%     disp('starting block');
    %cumTime = GetSecs;
    
    % start eyetracking...
    if p.eyeTracking == 1
        try
            start_recording_SMI();
        catch
            p.eyeTracking = 0;
        end

    elseif p.eyeTracking == 2
        try
            write_parallel(192); % start recording
            disp('recording started');
        catch thisError
            disp('Error during recording initiation');
            p.eyeTracking = 0;
            p.eyeTrackingStopped = 'After 192';
            disp(thisError);
        end
    end
    
    FlushEvents;
    WaitSecs(0.25);
    GetChar;
    
    
    
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
    p.startExp = GetSecs;
    
    while GetSecs < (p.startExp + p.startWait)
        [resp, timeStamp] = checkForResp([],p.escape);
        if resp == p.escape
            Screen('CloseAll');
            clear mex;
            return;
        end
        
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        
        %continue;
    end
    
    fprintf('\n');
    cumTime = GetSecs;
    
    save(fName,'p');
    
    %% here is the start of the trial loop
    for t=1:p.nt
        
        % REMEMBERED TARGET - this is always around location in p.conditions(trial,1)
        xLocTR = stimLocsX(t,p.conditions(t,1)); % remembered target - this is easy, from
        yLocTR = stimLocsY(t,p.conditions(t,1));
        
        xLocTR = xLocTR + rand(1) * p.jitterRadPix * cos(rand(1) * 2*pi); % + random radius 0:p.jitterRadPix * cos(random angle (0:2*pi))
        yLocTR = yLocTR + rand(1) * p.jitterRadPix * sin(rand(1) * 2*pi); % - random radius 0:p.jitterRadPix * sin(random angle (0:2*pi))
        
        % 1 = 60 ang deg 
        % 2 = 120 ang deg
        % 3 = 180 ang deg
  
        % let's put this up in the top...
%         p.conditions(t,5) = -1+2*round(rand(1)); % 1 = CCW, -1 = CCW
%         
%         TNloc = mod(p.conditions(t,1) + p.conditions(t,5) * p.conditions(t,3),6);
%         if TNloc == 0
%             TNloc = 6;
%         end
%         
%         if p.conditions(3) == 3
%             p.conditions(t,5) = 0;
%         end
%         
%         p.TRLoc(t) = p.conditions(t,1);
%         p.TNLoc(t) = TNloc;
        
        xLocTN = stimLocsX(t,p.TNloc(t));
        yLocTN = stimLocsY(t,p.TNloc(t));
        
        ang = rand(1); rad = rand(1);
        xLocTN = xLocTN + rad * p.jitterRadPix * cos(ang * 2*pi); % + random radius 0:p.jitterRadPix * cos(random angle (0:2*pi))
        yLocTN = yLocTN - rad * p.jitterRadPix * sin(ang * 2*pi); % - random radius 0:p.jitterRadPix * sin(random angle (0:2*pi))
        clear ang rad;
        
        p.TRcoord(t,:) = [xLocTR yLocTR];
        p.TNcoord(t,:) = [xLocTN yLocTN];
        
        p.targetDistancePix(t) = sqrt((xLocTN-xLocTR)^2 + (yLocTN-yLocTR)^2);
        p.targetDistanceDeg(t) = p.targetDistancePix(t)/p.ppd;
        
        % if we use up/down or left/right, we'll need to restrict the angle
        % a bit - for now writing for "which dot?" task
        ang = rand(1) * 2*pi;
        p.probeCoord(t,:) = p.TRcoord(t,:) + p.conditions(t,4) * [cos(ang), -1*sin(ang)];  % rad * rand angle... p.conditions,4 <-- is in degrees, so mult by ppd
        clear ang;
        
        % mark trial
        if p.eyeTracking == 1
            % setup data structure
            et(t).gaze = -10000*ones(round(p.trialDur(t)*100),2);
            et(t).fix_dist = -10000*ones(round(p.trialDur(t)*100),1);
            et(t).bad_frames = 0;
            et(t).epoch = zeros(round(p.trialDur(t)*100),1);
            et(t).t = -1*ones(round(p.trialDur(t)*100),1);
            et_cnt = 1;
            startTime = GetSecs;
            
            try
                send_trialNum_SMI(t);
            catch
                p.eyeTracking = 0;
                disp(sprintf('error writing trial num %02.f',trial_eachrun));
            end

        elseif p.eyeTracking == 2
            try
                write_parallel(192+t); % event ID will be trial number for now
                disp(sprintf('Start of trial %i',t));
            catch thisError
                disp(thisError);
                p.eyeTracking = 0;
                p.eyeTrackingStopped = sprintf('Start of trial %i (%i)',t,192+t);
            end
        end
        
        % randomly shuffle color correspondence...
        
        p.targetColorMapping(t,:) = randperm(2);  % p.wmTargetColor(colorMapping(1),:) is color of TR, 2 = TN
        p.probeColorMapping(t,:)  = randperm(2);  % p.wmTargetColor(colorMapping(1),:) is color of TR, 2 = TN
        
        p.trialStart(t) = GetSecs;
        
        %% SHOW TARGET
        p.stimStart(t) = GetSecs;
        while GetSecs <= p.trialStart(t) + p.wmTargetDur(t)
            
            
            if p.eyeTracking == 1

                try
                    [gx gy ex ey] = get_current_gaze_SMI(pSampleData,pEventData);

                    if p.show_eye_pos
                        Screen('DrawDots',w,[gx gy],5,[255 0 0]);
                    end

                    et(t).gaze(et_cnt,:) = [gx gy];
                    fix_dist = sqrt((gx-xcenter)^2 + (gy-ycenter)^2);
                    et(t).fix_dist(et_cnt) = fix_dist;

                    et(t).epoch(et_cnt) = 1; % TARGETS

                    et(t).t(et_cnt) = GetSecs - startTime;

                    if fix_dist > p.fixWinPix % radius
                        et(t).bad_frames = et(t).bad_frames + 1;
                    end

                    et_cnt = et_cnt+1;
                catch
                    disp('eyetracking broke (pre-cue)...');
                    lasterror.stack.file
                    lasterror.stack.name
                    lasterror.stack.line

                    which initIViewXAPI
                    which get_current_gaze_SMI
                    which iViewX.dll

                    p.eyeTracking = 0;
                end
            end
           
            if draw_target_frames == 1
                
                myc = hsv(p.nLoc);
                
                for i = 1:p.nLoc
                    stimRect = CenterRectOnPoint([0 0 2*p.jitterRadPix 2*p.jitterRadPix],stimLocsX(t,i),stimLocsY(t,i));
                    Screen('FrameOval',w,255*myc(i,:),stimRect);
                end
            end
            
            % TARGETs
            % TR - get its color from p.targetColorMapping(t,1)
            Screen('DrawDots',w,[p.TRcoord(t,:)],p.targetSizePix,p.wmTargetColors(p.targetColorMapping(t,1),:),[],1);
            % TN
            Screen('DrawDots',w,[p.TNcoord(t,:)],p.targetSizePix,p.wmTargetColors(p.targetColorMapping(t,2),:),[],1);

            % fixation
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0);

            Screen('DrawingFinished',w);
            Screen('Flip',w);
            
            [resp, timeStamp] = checkForResp(p.escape, p.escape);
            
            if resp == -1
                
                if p.eyeTracking == 1
                    stop_recording_SMI();
                    save_eyedata_SMI(p.et_fname);
                    end_session_SMI();
                end
                
                ListenChar(0); return;
                
            end
        end
        
        p.delay1Start(t) = GetSecs;   % start a clock to get the stim onset time
        p.stimEnd(t) = GetSecs;
        %% DELAY PERIOD 1
        while GetSecs<=(p.trialStart(t) + p.wmTargetDur(t) + p.wmDelay1Dur(t)) % if we want multiple exposure durations, add that here
            
            
            % draw fixation == memory cue color
            if p.conditions(t,2) == 1
                Screen('DrawDots', w, [0,0], p.fixSizePix, p.wmCueColors(p.targetColorMapping(t,1),:), center, 0);
            elseif p.conditions(t,2) == 2
                Screen('DrawDots', w, [0,0], p.fixSizePix, p.wmCueColors(3,:), center, 0);
            else  % 0, hopefully...
                Screen('DrawDots', w, [0,0], p.fixSizePix, p.wmCueColors(3,:), center, 0);
            end
            
            if db_probe
                Screen('DrawDots',w,[p.TRcoord(t,:)],p.targetSizePix,p.wmTargetColors(p.targetColorMapping(t,1),:),[],1);
            end
            
            % like this for eyetracking later....
            
            [kbResp, timeStamp] = checkForResp(p.keys, p.escape);
            if kbResp == -1
                if p.eyeTracking == 1
                    stop_recording_SMI();
                    save_eyedata_SMI(p.et_fname);
                    end_session_SMI();
                end
                ListenChar(0); return;
            end
            
            % eyetrack...
            if p.eyeTracking == 1

                try
                    [gx gy ex ey] = get_current_gaze_SMI(pSampleData,pEventData);

                    if p.show_eye_pos
                        Screen('DrawDots',w,[gx gy],5,[255 0 0]);
                    end

                    et(t).gaze(et_cnt,:) = [gx gy];
                    fix_dist = sqrt((gx-xcenter)^2 + (gy-ycenter)^2);
                    et(t).fix_dist(et_cnt) = fix_dist;

                    et(t).epoch(et_cnt) = 2; % DELAY PERIOD

                    et(t).t(et_cnt) = GetSecs - startTime;

                    if fix_dist > p.fixWinPix
                        et(t).bad_frames = et(t).bad_frames + 1;
                    end

                    et_cnt = et_cnt+1;
                catch
                    disp('eyetracking broke (delay)...');
                    %lasterror.stack.file
                    %lasterror.stack.name
                    %lasterror.stack.line

                    %which initIViewXAPI
                    %which get_current_gaze_SMI
                    %which iViewX.dll

              %      p.eyeTracking = 0;
                end
            end

            
            Screen('DrawingFinished',w);
            Screen('Flip',w);
            %frmCnt = frmCnt + 1;
        end
        
        p.delay1End(t) = GetSecs;
        
        
        %% SECOND DELAY (post-retro-cue)
        
        while GetSecs<=(p.trialStart(t) + p.wmTargetDur + p.wmDelay1Dur + p.wmDelay2Dur) % if we want multiple exposure durations, add that here
            
            % TODO: maybe make this a ring around fixation?
            % draw fixation == midtrial retro cue color
            if p.conditions(t,2) == 3 % fixation is the cued color
                Screen('DrawDots', w, [0,0], p.fixSizePix, p.wmCueColors(p.targetColorMapping(t,1),:), center, 0);
            else % R1 or R2hold
                Screen('DrawDots', w, [0,0], p.fixSizePix, p.wmCueColors(4,:), center, 0);
            end
            
            
            if db_probe
                Screen('DrawDots',w,[p.TRcoord(t,:)],p.targetSizePix,p.wmTargetColors(p.targetColorMapping(t,1),:),[],1);
            end
            
            % like this for eyetracking later....
            
            [kbResp, timeStamp] = checkForResp(p.keys, p.escape);
            if kbResp == -1
                
                
                if p.eyeTracking == 1
                    stop_recording_SMI();
                    save_eyedata_SMI(p.et_fname);
                    clear pSampleData pEventData;
                    end_session_SMI();
                end
                
                ListenChar(0); return;
                
                
                
            end
            
            % eyetrack...
            if p.eyeTracking == 1
                
                try
                    [gx gy ex ey] = get_current_gaze_SMI(pSampleData,pEventData);
                    
                    if p.show_eye_pos
                        Screen('DrawDots',w,[gx gy],5,[255 0 0]);
                    end
                    
                    et(t).gaze(et_cnt,:) = [gx gy];
                    fix_dist = sqrt((gx-xcenter)^2 + (gy-ycenter)^2);
                    et(t).fix_dist(et_cnt) = fix_dist;
                    % TODO: update epochs!
                    et(t).epoch(et_cnt) = 2; % DELAY PERIOD
                    
                    et(t).t(et_cnt) = GetSecs - startTime;
                    
                    if fix_dist > p.fixWinPix
                        et(t).bad_frames = et(t).bad_frames + 1;
                    end
                    
                    et_cnt = et_cnt+1;
                catch
                    disp('eyetracking broke (delay)...');
                    %lasterror.stack.file
                    %lasterror.stack.name
                    %lasterror.stack.line
                    
                    %which initIViewXAPI
                    %which get_current_gaze_SMI
                    %which iViewX.dll
                    
                    %      p.eyeTracking = 0;
                end
            end
            
            Screen('DrawingFinished',w);
            Screen('Flip',w);
            
        end
        
        p.delay2End(t) = GetSecs;
        startResponse = p.delay2End(t);
        
        %% PRESENT RECALL STIMULUS
        % fix will be color of response cue (below)
        %
        % p.respCoord = 1 for Left/Right, 2 for Up/Down
        
        if p.respCoord(t) == 1   % if left/right
            respLineCoords = [0 -p.sRect(3)/2; 0 p.sRect(3)/2]';
            to_add = [1 0; 1 0]';
        else
            respLineCoords = [-p.sRect(4)/2 0; p.sRect(4)/2 0]';
            to_add = [0 1; 0 1]';
        end
        
        fidx = 1;
        while GetSecs <= (p.trialStart(t) + p.wmTargetDur(t) + p.wmDelay1Dur(t) + p.wmDelay2Dur(t) + p.responseWindow(t))
            
            Screen('DrawLines',w,respLineCoords,1,p.wmRecallColor,center);

            % draw cue
            
            if p.conditions(t,2) == 0 % r0
                Screen('DrawDots',w,[0 0], p.fixSizePix, p.wmCueColors(4,:),center,0); % black cue
            else
                Screen('DrawDots',w,[0 0], p.fixSizePix, p.wmCueColors(p.targetColorMapping(t,1),:),center,0);
            end

            [kbResp, timeStamp] = checkForResp(p.keys, p.escape);
            
            if p.eyeTracking == 1

                try
                    [gx gy ex ey] = get_current_gaze_SMI(pSampleData,pEventData);

                    if p.show_eye_pos
                        Screen('DrawDots',w,[gx gy],5,[255 0 0]);
                    end

                    et(t).gaze(et_cnt,:) = [gx gy];
                    fix_dist = sqrt((gx-xcenter)^2 + (gy-ycenter)^2);
                    et(t).fix_dist(et_cnt) = fix_dist;

                    et(t).epoch(et_cnt) = 3; % RESPONSE PERIOD

                    et(t).t(et_cnt) = GetSecs - startTime;

                    if fix_dist > p.fixWinPix
                        et(t).bad_frames = et(t).bad_frames + 1;
                    end

                    et_cnt = et_cnt+1;
                catch
                    disp('eyetracking broke (response)...');
                    p.eyeTracking = 0;
                end
            end
            
            if kbResp && find(p.keys==kbResp)

                which_resp = find(p.keys==kbResp);
                if which_resp == 1 % go left or up, so negative
                    respLineCoords = respLineCoords - to_add * p.stepSizePix;
                elseif which_resp == 2
                    respLineCoords = respLineCoords + to_add * p.stepSizePix;
                end

            end

            
            if kbResp == -1               
                if p.eyeTracking == 1
                    stop_recording_SMI();
                    save_eyedata_SMI(p.et_fname);
                    end_session_SMI();
                end
                ListenChar(0); return;
            end
            
            p.respCoordTrace(t,fidx) = respLineCoords(p.respCoord(t),1);
            fidx = fidx + 1;
            
            Screen('DrawingFinished',w);
            Screen('Flip',w);

        end
        
        clear fidx;
        
        Screen('DrawDots',w,[0 0], p.fixSizePix, p.wmCueColors(p.targetColorMapping(t,1),:),center,0);
        
        Screen('DrawingFinished',w);
        Screen('Flip',w);
        

        p.wmResp(t) = respLineCoords(p.respCoord(t),1);
        p.wmRecallDistanceDeg(t) = ((p.wmResp(t)) - (p.TRcoord(t,p.respCoord(t)) - center(p.respCoord(t)) ))/p.ppd;
        
        
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %draw fixation point
        Screen('Flip',w);
        
        
        %% ITI
        while GetSecs <= p.trialStart(t)  + p.wmTargetDur(t) + p.wmDelay1Dur(t) + p.wmDelay2Dur(t) +  p.responseWindow(t) + p.ITI(t)
            
            
            [resp, timeStamp] = checkForResp(p.keys, p.escape);
            if resp==-1; 
                if p.eyeTracking == 1
                    stop_recording_SMI();
                    save_eyedata_SMI(p.et_fname);
                    end_session_SMI();
                end
                ListenChar(0); return; 
            end;

            if p.eyeTracking == 1

                try
                    [gx gy ex ey] = get_current_gaze_SMI(pSampleData,pEventData);

                    if p.show_eye_pos
                        Screen('DrawDots',w,[gx gy],5,[255 0 0]);
                    end

                    et(t).gaze(et_cnt,:) = [gx gy];
                    fix_dist = sqrt((gx-xcenter)^2 + (gy-ycenter)^2);
                    et(t).fix_dist(et_cnt) = fix_dist;

                    et(t).epoch(et_cnt) = 5; % ITI

                    et(t).t(et_cnt) = GetSecs - startTime;

                    if fix_dist > p.fixWinPix
                        et(t).bad_frames = et(t).bad_frames + 1;
                    end

                    et_cnt = et_cnt+1;
                catch
                    disp('eyetracking broke (response)...');
                    lasterror.stack.file
                    lasterror.stack.name
                    lasterror.stack.line

                    which initIViewXAPI
                    which get_current_gaze_SMI
                    which iViewX.dll

                    p.eyeTracking = 0;
                end

            end
            
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %draw fixation point
            Screen('Flip',w);

        end

%         if p.probeColorMapping(t,1)== p.wmResp(t) && ~isnan(p.wmResp(t))
%             p.correct(t) = 1;
%         elseif p.probeColorMapping(t,2)== p.wmResp(t) && ~isnan(p.wmResp(t))
%             p.correct(t) = 0;
%         end
%         
        p.trialEnd(t) = GetSecs;
        
        if exist('et','var')
            save([fName], 'p','et')

        else
            save([fName], 'p');
        end
        
        
        
    end
    %% end trial loop
    
    %   10s passive fixation at end of block
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    startPassive = GetSecs;
    while GetSecs <= startPassive + p.passiveDur
        [resp, timeStamp] = checkForResp(p.escape, p.escape);
        if resp==-1; ListenChar(0); return; end;
    end
    
    p.endExp = GetSecs;
    
    
    %save trial data from this block
    if exist('et','var')
        save([fName], 'p','et')
        
    else
        save([fName], 'p');
    end
    
    recall_err = nan(2,1);
    %uc = unique(p.conditions(~isnan(p.wmResp),2));
    for ii = 1:2
        recall_err(ii) = mean(abs(p.wmRecallDistanceDeg(p.respCoord==ii)));
    end
    
    str1 = sprintf('Block %i.%i complete',p.runNum,p.runSub);
    tCenter1 = [center(1)-RectWidth(Screen('TextBounds', w, str1))/2 center(2)/2];
    
    %str2 = sprintf('Accuracy: %2.1f %% / %2.1f %% / %2.1f %%',100*nanmean(p.correct(p.conditions(:,2)==1)),100*nanmean(p.correct(p.conditions(:,2)==2)),100*nanmean(p.correct(p.conditions(:,2)==3)));
    %str2 = sprintf('Accuracy: %2.1f%% / %2.1f%% / %2.1f%%',100*behav_acc(1),100*behav_acc(2),100*behav_acc(3));
    str2 = sprintf('Recall error (h/v): %.02f / %.02f',recall_err(1),recall_err(2));
    fprintf('%s\n',str2); % to command window so we don't lose it
    tCenter2 = [center(1)-RectWidth(Screen('TextBounds', w, str2))/2 center(2)/2];
    
    %str3 = sprintf('No response: %i / %i / %i',sum(isnan(p.wmResp(p.conditions(:,2)==1))),sum(isnan(p.wmResp(p.conditions(:,2)==2))),sum(isnan(p.wmResp(p.conditions(:,2)==3))));
    %fprintf('%s\n',str3); % to command window so we don't lose it
    %tCenter3 = [center(1)-RectWidth(Screen('TextBounds', w, str3))/2 center(2)/2];

    
    Screen('DrawText', w, str1, tCenter1(1), tCenter1(2)-100, p.textColor);
    Screen('DrawText', w, str2, tCenter2(1), tCenter2(2), p.textColor);
    %Screen('DrawText', w, str3, tCenter3(1), tCenter3(2)+100, p.textColor);
    
    % put up a message to wait for a space bar press.
    Screen('Flip', w);
    
 
    while resp~=p.space
        [resp, timeStamp] = checkForResp(p.space, p.escape);
        if resp==-1; ListenChar(0); return; end;
    end
    
    % stop recording and save
    if p.eyeTracking == 1
        try
            stop_recording_SMI();
            save_eyedata_SMI(et_fname);
            end_session_SMI()
        catch
            p.eyeTracking = 0;
        end
    elseif p.eyeTracking == 2
        try
            write_parallel(0); % close file, end recording
        catch thisError
            disp(thisError);
            p.eyeTracking = 0;
            p.eyeTrackingStopped = sprintf('End of block %i',p.runNum);
        end
    end
    
try
ListenChar(0);
end
Screen('CloseAll');

return


%--------------------------------------------------------------------------
function p = deg2pix(p)
% converts degrees visual angle to pixel units before rendering
% with PTB. Needs p.screenWidthCM and p.vDistCM
% js - 10.2007

% figure out pixels per degree, p.sRect(1) is x coord for upper left of
% screen, and p.sRect(3) is x coord for lower right corner of screen
p.ppd = pi * (p.sRect(3)-p.sRect(1)) / atan(p.screenWidthCM/p.vDistCM/2) / 360;

% get name of each field in p
s = fieldnames(p);

% convert all fields with the word 'Deg' from degrees visual angle to
% pixels, then store in a renmaed field name
for i=1:length(s)
    ind = strfind(s{i}, 'Deg');
    if ind
        curVal = getfield(p,s{i});
        tmp = char(s{i});
        newfn = [tmp(1:ind-1), 'Pix'];
        p = setfield(p,newfn,curVal*p.ppd);
    end
end
