% 372 s + 5 TRs - (171 TRs @ 2250 ms)

function wmDrop_localizer %manual

%
% this will be used for masking of spatial responses in early regions

addpath('../parallel_tools');


%% get user input about sub name etc...

mygreen = [0 220 0];
myred   = [255 0 0];
mygray = [65 65 65]; % fix this one

warning('off','MATLAB:dispatcher:InexactMatch');
%Screen('Preference', 'SkipSyncTests', 1)
%get subject info
prompt = {'Subject Name', 'Block number', 'seed','Display','Target/probe distance (dva)','Eye tracking'};
%grab a random number and seed the generator (include this as a gui field
%in case we want to repeat an exact sequence)
s = round(sum(100*clock));
%put in some default answers
defAns = {'XX','1',num2str(s),'1','0.8','0'};
box = inputdlg(prompt,'Enter Subject Information...', 1, defAns);

p.exptName = 'wmDrop_localizer';

if length(box)==length(defAns)
    p.subName=char(box{1});
    p.sessionNum=3;%str2num(box{2});
    %p.cond = str2num(box{2});
    %p.staggered = str2num(box{3});
    p.nBlocks=eval(box{2});
    p.rndSeed=str2num(box{3});
    p.display=str2num(box{4});
    p.wmDifficultyDeg = str2num(box{5})
    p.eyeTracking = str2num(box{6});
    rand('state',p.rndSeed);  %actually seed the random number generator
else
    return
end

%p.keys = 66; %keys that correspond to elements in p.dir respectively (B)

if p.eyeTracking
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

%% --------------------begin user define parameters----------------------------%

% Screen (parameter) setup
%p.windowed = ~p.fMRI;                     % if 1 then the display will be in a small window, use this for debugging.
p.windowed = 0;
% monitor stuff
p.refreshRate = 60;                 % refresh rate is normally 100 but change this to 60 when on laptop!

if p.display == 1 % projector in scanner
    % 3TE
    p.fMRI = 1;
    p.vDistCM = 375;
    p.screenWidthCM = 120;
    p.keys = [KbName('b'),KbName('y')]; % 1 and 2 keys
elseif p.display == 3 % behavioral eyetracker room
    p.vDistCM = 62;
    p.screenWidthCM = 51;
    p.percentNull = 0;
    p.fMRI = 0;
    p.keys = [KbName('left'),KbName('right')]; % 1 and 2 keys
end

%p.nLoc = 6; % stim will appear at random location on nLoc by nLoc grid (should be even)

%p.radDeg = p.usedScreenSizeDeg/(p.nLoc+1);

%p.stimSizeDeg = p.usedScreenSizeDeg/2 + p.radDeg/4; % due to stagger each way...

%stimulus geometry (in degrees, note that these vars have the phrase 'Deg'
%in them, used by pix2deg below)
%p.radDeg = p.usedScreenSizeDeg/(p.nLoc+1);
p.sfDeg  = 1.5;
% p.targetProbeMinSepDistance = pi/12; % radians (WM)
% p.targetProbeMaxSepDistance = pi/8;
% p.targetProbeMinRad = 0.2; % innerRad + targProbeMinRad * stim_width (out - in)
% p.targetProbeMaxRad = 0.8; % 
p.wmProbeWindowAng  = 20; % up to plus/minus 20 deg from horiz/vertical offset from target
p.wmTargetWindowAng = 60; % 60 deg polar angle up or down from horizontal meridian
p.wmTargetWindowEccDeg = [2.5 4.5]; % dva from fixation
% 
p.probeCueSizeDeg = 0.5; % the bar that indicates the judgment to make
p.wmRespCueColor = mygray;


p.nAng = 60;  %number of sub-wedges in a circle - RADIAL checkerboard
p.nRad = 15;  %number of sub-rings in a circle

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

p.space = KbName('space');
p.start = KbName('t');

p.backColor     = [128, 128, 128];      % background color

%fixation point properties
p.fixColor      = [180, 180, 180];
p.fixSizeDeg    = .2;                  % size of dot

% target/probe properties (WM)
p.targetSizeDeg = 0.35;
p.probeSizeDeg = 0.35;
p.wmTargetColor = myred;
p.wmProbeColor  = mygreen;

% stimulus properties (IN goes from fixAppRad to inStimRad, OUT goes from
% inStimRad to outStimRad)
p.fixAppSizeDeg    = p.fixSizeDeg * 4;     % size of apperture surrounding dot

p.inStimRadiusDeg = p.fixAppSizeDeg;
p.outStimRadiusDeg = 6; % dva - 3 * p.stepSize + p.radDeg (from wmDrop_hexMap.m) - for now using 1.5 and 1, respectively
% target properties (fix/checks)
p.stimContrast = 1;

% trial setup
p.percentNull = 0.2;
p.repetitions = 8; % FIX THIS!!!! 
if p.windowed % DEBUG MODE
    p.repetitions = 1;
end
p.nCond = 2; % L, R
p.nValidTrials = p.nCond * p.repetitions;    % length p.cond, for now, should always be 1
p.numNullTrials = ((1/(1-p.percentNull)) - 1) * p.nValidTrials;
p.nTrials = p.nValidTrials + p.numNullTrials; % no stim


p.TR = 2.25;
p.nTRsWait = 0;

% stimulus timing (in seconds - convert to vid frames, if necessary, later)
p.stimExposeDur = 10; % keep this constant
p.flickerFreq   = 6; % Hz
p.flickerPer  = 1/p.flickerFreq; % s
p.wmTargetDur   = 0.5; % s
p.wmProbeDur    = 0.75; % s
p.wmTargProbeSepDur = 3.0; % s
% p.wmTargDurPer = p.wmTargetDur * p.flickerFreq;
% p.wmProbeDurPer = p.wmProbeDur * p.flickerFreq;
% p.wmTargProbeSepPer = p.wmTargProbeSepDur * p.flickerFreq;

p.wmResponseWindow = 1.5;
p.nTargs = 1;

p.wmTargWindow = [1.0 p.stimExposeDur-1-p.wmTargProbeSepDur-p.wmProbeDur-p.wmTargetDur-p.wmResponseWindow]; % period over which contrast can change
%p.minTargPer = p.wmTaskWindow(1)*p.flickerFreq + 1; % these are period index (which period can we start target?)
%p.maxTargPer = (p.wmTaskWindow(2)-p.wmTargProbeSepDur-p.wmTargetDur)*p.flickerFreq + 1; % +1 is because want cycle *starting at* this point in time
%p.minTargSepDur = 1.5; % separation b/w targ/probe pair "trials"
%p.minTargSepPer = floor(p.wmTargProbeSepDur*p.flickerFreq); % number of periods (12 = 2 s)
% time between stimulus presentations, constant
p.ITI           = linspace(3,5,p.nTrials); p.ITI = p.ITI(randperm(p.nTrials));%3.0 * ones(p.nTrials,1);
p.trialDur      = p.stimExposeDur + p.ITI;
p.startWait     = 2; % may change...
p.passiveDur    = 10; % in sec, how long to wait at end of a block (and beginning?)
p.expDur        = sum(p.trialDur) + p.passiveDur + p.startWait + p.TR*p.nTRsWait;

disp(sprintf('SCAN DURATION: %i',p.expDur));
disp(sprintf('# TRs: %i',ceil(p.expDur/p.TR)));

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
else
    [w p.sRect]=Screen('OpenWindow', s, p.backColor);
    HideCursor;
end

p.aspectRatio = p.sRect(4)/p.sRect(3);

% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% test the refesh properties of the display
p.fps=Screen('FrameRate',w);          % frames per second
p.ifi=Screen('GetFlipInterval', w);   % inter-frame-time
if p.fps==0                           % if fps does not register, then set the fps based on ifi
    p.fps=1/p.ifi;
end
p.flickerFrames = 1/p.flickerFreq*p.fps;
p.targetFrames = p.wmTargetDur * p.fps;
p.probeFrames = p.wmProbeDur * p.fps;

% make sure the refreshrate is ok
if abs(p.fps-p.refreshRate)>5
    Screen('CloseAll');
    disp('CHANGE YOUR REFRESH RATE')
    ListenChar(0);
    clear all;
    return;
end

% convert relevant timing stuff to vid frames for stim presentation
p.stimExpose = round((p.stimExposeDur*1000)./(1000/p.refreshRate));

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

p.xOffsetPix = 0;
left = [center(1)-p.xOffsetPix, center(2)];

Screen('TextSize', w, p.fontSize);
Screen('TextStyle', w, 1);
Screen('TextFont',w, p.fontName);
Screen('TextColor', w, p.textColor);

% file/directory operations
p.root = pwd;

if ~isdir([p.root '/Subject Data/'])
    mkdir([p.root '/Subject Data/']);
end
p.subNameAndDate = [p.subName '_' datestr(now,30)];

%% start a block loop
for b=p.nBlocks

    fName=[p.root, '/data/',  p.subName,'_',p.exptName, '_sess', num2str(p.sessionNum), '_blkNum', num2str(b), '.mat'];
    if exist(fName,'file')
        Screen('CloseAll');
        msgbox('File name already exists, please specify another', 'modal')
        ListenChar(0);
        return;
    end
    
    % open eyetracker file
    if p.eyeTracking
        try
            write_parallel(64); % open file             
            disp('ET started');
        catch thisError
            disp('Error during file opening');
            p.eyeTracking = 0;
            p.eyeTrackingStopped = 'After 64';
            disp(thisError);
        end
    end
    
    p.respCoord = round(rand(p.nValidTrials,p.nTargs))+1; % LR or UPDOWN
    p.corrResp  = round(rand(p.nValidTrials,p.nTargs))+1; % L/Up or R/Down
    
    % choose which trials are same/different target/probe location
    %p.targetProbeMatch = [zeros(floor(p.nValidTrials*p.nTargs/2),1); ones(ceil(p.nValidTrials*p.nTargs/2),1)]; % if floor/ceil used reliably together, should always come out same... ( sum = p.nTrials)
    
    p.targetLocR = p.wmTargetWindowEccPix(1) + rand(p.nValidTrials,p.nTargs) * (p.wmTargetWindowEccPix(2)-p.wmTargetWindowEccPix(1));
    p.targetLocT = rand(p.nValidTrials,p.nTargs) * 2 * p.wmTargetWindowAng - p.wmTargetWindowAng;
        % p.wmTargetWindowAng = 60; % 60 deg polar angle up or down from horizontal meridian
%    p.wmTargetWindowEccDeg = [1.5 4]; % dva from fixation
    
    % r,theta
%     p.targetLocR =  p.targetProbeMinRad + (p.targetProbeMaxRad - p.targetProbeMinRad)*rand(p.nValidTrials*p.nTargs,1); % pick radii
%     p.targetLocT = [linspace(p.targetProbeMinSepDistance, pi-p.targetProbeMinSepDistance,floor(p.nValidTrials*p.nTargs/2))';... % hemifield..
%         linspace(p.targetProbeMinSepDistance, pi-p.targetProbeMinSepDistance,ceil(p.nValidTrials*p.nTargs/2))'];


    p.probeRelLocR = p.wmDifficultyPix * ones(p.nValidTrials,p.nTargs);
    p.probeRelLocT = rand(p.nValidTrials,p.nTargs) * 2 * p.wmProbeWindowAng - p.wmProbeWindowAng; % degrees
    p.probeRelLocT(p.respCoord==2) = p.probeRelLocT(p.respCoord==2)-90; % UP/DOWN
    p.probeRelLocT(p.corrResp==1) = p.probeRelLocT(p.corrResp==1)+180; % UP or LEFT
    
    
    % randomize trial order (w/ ITI)
    % 1 = LEFT, 2 = RIGHT
    
    p.stim = [ones(1,p.repetitions) 2*ones(1,p.repetitions)]';
    p.stim = Shuffle(p.stim);
    
    % should:
    % - generate random sequence of 1/2
    % - 2 s ITI
    % - insert null trials randomly, as before
    
    p.flickerSequ = repmat([ones(1,p.flickerFrames/2) 2*ones(1,p.flickerFrames/2)],1,p.stimExpose/p.flickerFrames);
    %p.flickerSequ(2,:) = p.flickerSequ(1,:) + 2;
    %p.stimSequ = repmat(p.flickerSequ,p.nTrials,1);
    p.wmTargOnset = rand(p.nValidTrials,p.nTargs)*(p.wmTargWindow(2)-p.wmTargWindow(1)) + p.wmTargWindow(1);
    p.wmTargOnsetFr = floor(p.wmTargOnset/(p.flickerFrames*p.ifi/2));
    p.wmTargOnset = p.wmTargOnsetFr*p.ifi*p.flickerFrames/2;
    
    
    

    
    % NULL   
    tmp = Shuffle(1:p.nValidTrials-1);
    p.nullIdx = sort(tmp(1:p.numNullTrials),'descend'); % insert null trials *after* these
    
    % insert null trials
    for i = 1:length(p.nullIdx)
        p.stim(p.nullIdx(i)+2:end+1) = p.stim(p.nullIdx(i)+1:end);
        p.stim(p.nullIdx(i)+1) = NaN;
%        p.wmSequ(p.nullIdx(i)+2:end+1,:) = p.wmSequ(p.nullIdx(i)+1:end,:);
%        p.wmSequ(p.nullIdx(i)+1,:) = NaN;
%        p.targetProbeMatch(p.nullIdx(i)+2:end+1,:) = p.targetProbeMatch(p.nullIdx(i)+1:end,:);
        p.targetLocR(p.nullIdx(i)+2:end+1,:) = p.targetLocR(p.nullIdx(i)+1:end,:);
        p.targetLocT(p.nullIdx(i)+2:end+1,:) = p.targetLocT(p.nullIdx(i)+1:end,:);
        p.probeRelLocR(p.nullIdx(i)+2:end+1,:) = p.probeRelLocR(p.nullIdx(i)+1:end,:);
        p.probeRelLocT(p.nullIdx(i)+2:end+1,:) = p.probeRelLocT(p.nullIdx(i)+1:end,:);
        p.respCoord(p.nullIdx(i)+2:end+1,:) = p.respCoord(p.nullIdx(i)+1:end,:);
        p.corrResp(p.nullIdx(i)+2:end+1,:) = p.corrResp(p.nullIdx(i)+1:end,:);
        p.wmTargOnset(p.nullIdx(i)+2:end+1,:) = p.wmTargOnset(p.nullIdx(i)+1:end,:);
        
        p.targetProbeMatch(p.nullIdx(i)+1,:) = NaN;
        p.targetLocR(p.nullIdx(i)+1,:) = NaN;
        p.targetLocT(p.nullIdx(i)+1,:) = NaN;
        p.probeRelLocR(p.nullIdx(i)+1,:) = NaN;
        p.probeRelLocT(p.nullIdx(i)+1,:) = NaN;
        p.respCoord(p.nullIdx(i)+1,:) = NaN;
        p.corrResp(p.nullIdx(i)+1,:) = NaN;
        p.wmTargOnset(p.nullIdx(i)+1,:) = NaN;
    end
    
    
    p.null = isnan(p.stim);
    
    
    p.stimStart =  nan(p.nTrials,1);
    p.stimEnd =  nan(p.nTrials,1);
    p.wmRespT  =  nan(p.nTrials,p.stimExpose);
    p.wmResp   = nan(p.nTrials,1);
    p.responseTime = nan(p.nTrials,1);
    
    p.probeStart = nan(p.nTrials,1);
    % generate checkerboards we use...
    %ch = make_checkerboard_rect(p.stimSizePix,2*p.stimSizePix,p.sfPix,p.stimContrast);
    
    %c_in = make_checkerboard(p.inStimRadiusPix,p.sfInPix,p.stimContrast);
    ch = make_radial_checkerboard(p);
    
    stim(1)=Screen('MakeTexture', w, ch{1});
    stim(2)=Screen('MakeTexture', w, ch{2});

    
    %c_in = make_checkerboard(p.inStimRadiusPix,p.sfInPix,p.stimContrast);
    %ch = make_radial_checkerboard(p);
    
    stim(1)=Screen('MakeTexture', w, ch{1});
    stim(2)=Screen('MakeTexture', w, ch{2});
    %c_out = make_checkerboard(p.outStimRadiusPix,p.sfOutPix,p.stimContrast);
    
    %[x y] = meshgrid(linspace(-p.outStimRadiusPix,p.outStimRadiusPix,size(ch{1},1)),linspace(-p.outStimRadiusPix,p.outStimRadiusPix,size(ch{1},1)));
    %ch_r = sqrt(x.^2 + y.^2);
    %ch{1}(ch_r > p.inStimRadiusPix) = p.backColor(1);
    %ch{2}(ch_r > p.inStimRadiusPix) = p.backColor(1);
    %stim(3)=Screen('MakeTexture', w, ch{1});
    %stim(4)=Screen('MakeTexture', w, ch{2});
    
    % put up a message to wait for a space bar press
    textMotion = 'Compare probe position to remembered position';

    tCenterMo = [center(1)-RectWidth(Screen('TextBounds', w, textMotion))/2 center(2)/2];
    
    Screen('DrawText', w, textMotion, tCenterMo(1), tCenterMo(2), p.textColor);
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, left, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
    %after all initialization is done, sit and wait for scanner synch (or
    %space bar)
    FlushEvents;
    WaitSecs(0.25);
    GetChar;
    

    % start eyetracking...
    if p.eyeTracking
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
    

    
    

    
    
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, left, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
    p.startExp = GetSecs;
    
    while GetSecs < p.startExp + p.startWait
        [resp, timeStamp] = checkForResp([],p.escape);
        if resp == p.escape
            Screen('CloseAll');
            clear mex;
            return;
        end
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, left, 0); %change fixation point
        Screen('DrawingFinished', w);
        Screen('Flip', w);
    end
    
    cumTime = GetSecs;
    
    rTcnt = 1;
    
    save(fName,'p');
    
    %% here is the start of the trial loop
    for t=1:p.nTrials
        
        p.trialStart(t) = GetSecs;
        
        disp(sprintf('STIM = %i',p.stim(t)));
        
        if p.eyeTracking        
            try
                write_parallel(192+t); % event ID will be trial number for now
                disp(sprintf('Start of trial %i',t));
            catch thisError
                disp(thisError);
                p.eyeTracking = 0;
                p.eyeTrackingStopped = sprintf('Start of trial %i (%i)',t,192+t);
            end
        end

        txtRect = Screen('Rect',stim(1));
        stimRect = CenterRect(txtRect,Screen('Rect',w));
        maskOval = CenterRect([0 0 2*p.fixAppSizePix 2*p.fixAppSizePix],Screen('Rect',w));
        %thisRadTarg = p.fixAppSizePix + (p.stimSizePix - p.fixAppSizePix) * p.targetLocR(t,:);
        %thisRadProbe = p.fixAppSizePix + (p.stimSizePix - p.fixAppSizePix) * p.probeLocR(t,:);
        
        targetLocX =  p.targetLocR(t).*cosd(p.targetLocT(t));
        targetLocY = -p.targetLocR(t).*sind(p.targetLocT(t));
        
        
        
        rct = Screen('Rect',w);
        if p.stim(t) == 1 % LEFT
            maskRect = [rct(3)/2 0 rct(3) rct(4)];
            targetLocX = -targetLocX;
            %stimRect = stimRect - [txtRect(3)/2 0 txtRect(3)/2 0];
        else              % RIGHT
            
            maskRect = [0 0 rct(3)/2 rct(4)];
            %stimRect = stimRect + [txtRect(3)/2 0 txtRect(3)/2 0];
        end
        
        if p.respCoord(t) == 1 % L/R
            bar_coord = [-p.probeCueSizePix p.probeCueSizePix; 0 0];
        else
            bar_coord = [0 0; -p.probeCueSizePix p.probeCueSizePix];
        end
        
        p.targetCoord(t,:) = [targetLocX targetLocY];
        
        probeLocX = targetLocX + p.probeRelLocR(t)*cosd(p.probeRelLocT(t));
        probeLocY = targetLocY - p.probeRelLocR(t)*sind(p.probeRelLocT(t));
        
        p.probeCoord(t,:) = [probeLocX probeLocY];
        
%         
%         if p.stim(t) == 1 || p.stim(t) == 3 % inner
%             txtRect  = CenterRect([0 0 2*p.inStimRadiusPix 2*p.inStimRadiusPix],Screen('Rect',stim(1))); % use center of radial checkerboard
%             stimRect = CenterRect(txtRect,Screen('Rect',w));
%             maskOval = CenterRect([0 0 2*p.fixAppRadiusPix 2*p.fixAppRadiusPix],Screen('Rect',w));
%             thisRadTarg = p.fixAppRadiusPix + (p.inStimRadiusPix - p.fixAppRadiusPix) * p.targetLocR(t,:);
%             thisRadProbe = p.fixAppRadiusPix + (p.inStimRadiusPix - p.fixAppRadiusPix) * p.probeLocR(t,:);
%         else  % outer
%             txtRect  = Screen('Rect',stim(1)); % use whole radial checkerboard
%             stimRect = CenterRect(txtRect,Screen('Rect',w));
%             maskOval = CenterRect([0 0 2*p.inStimRadiusPix 2*p.inStimRadiusPix],Screen('Rect',w));
%             thisRadTarg = p.inStimRadiusPix + (p.outStimRadiusPix - p.inStimRadiusPix) * p.targetLocR(t,:);
%             thisRadProbe = p.inStimRadiusPix + (p.outStimRadiusPix - p.inStimRadiusPix) * p.probeLocR(t,:);
%         end
%         
%         rct = Screen('Rect',w);
%         if p.stim(t) == 2
%             %maskRect = [0 0 rct(3)/2 rct(4)];
%             thisAngTarg = p.targetLocT(t,:) - pi/2; % if right side, rotate "top" hemi-circle to right
%             thisAngProbe = p.probeLocT(t,:) - pi/2;
%         else
%             %maskRect = [rct(3)/2 0 rct(3) rct(4)];
%             thisAngTarg = p.targetLocT(t,:) + pi/2; % otherwise, rotate it left
%             thisAngProbe = p.probeLocT(t,:) + pi/2;
%         end
%         clear rct;
        
%         % get x, y screen coords for targ/probe x & ys
%         targetLocX = thisRadTarg.*cos(thisAngTarg);
%         targetLocY = thisRadTarg.*sin(thisAngTarg);
%         probeLocX  = thisRadProbe.*cos(thisAngProbe);
%         probeLocY  = thisRadProbe.*sin(thisAngProbe);
%         
        frmCnt=1;
        p.stimStart(t) = GetSecs;   % start a clock to get the stim onset time
        
        % STIMULUS
        while frmCnt<=p.stimExpose % if we want multiple exposure durations, add that here
            
            if ~p.null(t)
                % stim: 1 = L_in, 2 = L_out, 3 = R_in, 4 = R_out

                Screen('DrawTexture',w, stim(p.flickerSequ(1,frmCnt)),txtRect,stimRect);
                Screen('FillOval',w,p.backColor,maskOval);
                Screen('FillRect',w,p.backColor,maskRect);
                

                tpt = GetSecs-p.stimStart(t);
                
                % present target
                
                if tpt >= p.wmTargOnset(t) && tpt <= (p.wmTargOnset(t)+p.wmTargetDur)
                    Screen('DrawDots',w,[targetLocX targetLocY],...
                            p.targetSizePix,p.wmTargetColor,left,1);
                end
                
                % prseent probe
                if tpt >= (p.wmTargOnset(t)+p.wmTargetDur+p.wmTargProbeSepDur) && tpt <= (p.wmTargOnset(t)+p.wmTargetDur+p.wmTargProbeSepDur+p.wmProbeDur)
                        Screen('DrawDots',w,[probeLocX probeLocY],...
                            p.probeSizePix,p.wmProbeColor,left,1);
                    % also drawlines....
                    Screen('DrawLines',w,bar_coord,3,p.wmRespCueColor,left);
                    
                    % start response clock
                    if isnan(p.probeStart(t))
                        p.probeStart(t) = GetSecs;
                    end
                    
                end

                
                Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, left, 0);
                
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                Screen('Flip', w);

                % check response...
                [resp, timeStamp] = checkForResp(p.keys, p.escape); % checks both buttons...
                if resp==-1; ListenChar(0); return; end;
                if resp && find(p.keys==resp)
                    p.wmRespT(t, frmCnt) = find(p.keys==resp);
                    if isnan(p.wmResp(t)) && tpt >=(p.wmTargOnset(t)+p.wmTargetDur+p.wmTargProbeSepDur) && tpt <= (p.wmTargOnset(t)+p.wmTargetDur+p.wmTargProbeSepDur+p.wmProbeDur+p.wmResponseWindow)
                        p.wmResp(t) = find(p.keys==resp);
                        p.responseTime(t) = GetSecs - p.probeStart(t);
                    end
                end
                
            end
            frmCnt = frmCnt + 1;
        end
        
        p.stimEnd(t) = GetSecs;
                
        
        % clear out screen
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, left, 0); %draw fixation point
        Screen('Flip',w);
        
        while GetSecs <= p.trialStart(t) + p.trialDur(t)
            [resp, timeStamp] = checkForResp(p.keys, p.escape);
            if resp==-1; ListenChar(0); return; end;
        end
 
        cumTime = cumTime + p.trialDur(t);
        rTcnt = rTcnt + 1;  % increment real trial counter.

        p.trialEnd(t) = GetSecs;
        
        save(fName, 'p');
        
    end
    %% end trial loop

    %   10s passive fixation at end of block
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, left, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    while GetSecs <= p.trialEnd(end) + p.passiveDur
        [resp, timeStamp] = checkForResp(p.escape, p.escape);
        if resp==-1; ListenChar(0); return; end;
    end
    
    p.endExp = GetSecs;
    %if p.cond == 3
    %    p.accuracy = p.wmCorrect/(p.wmCorrect+p.wmIncorrect)*100;
    %else
    %    p.accuracy = p.hits/sum(p.validTarget)*100;
    %    p.meanCorrectRT = nanmean(p.rt(find(p.validTarget)));
    %end
    
    %save trial data from this block
    save(fName, 'p');
   
    str1 = sprintf('Block %i complete',b);
    tCenter1 = [center(1)-RectWidth(Screen('TextBounds', w, str1))/2 center(2)/2];
    
    p.correct = p.wmResp == p.corrResp;
    acc = nanmean(p.wmResp(~p.null) == p.corrResp(~p.null));
    p.accuracy = acc;
    str2 = sprintf('Accuracy: %.02f%%',acc*100);
    fprintf('%s\n',str2);
    tCenter2 = [center(1)-RectWidth(Screen('TextBounds', w, str2))/2 center(2)/2];

    Screen('DrawText', w, str1, tCenter1(1), tCenter1(2)-100, p.textColor);
    Screen('DrawText', w, str2, tCenter2(1), tCenter2(2), p.textColor);


    % put up a message to wait for a space bar press.
    Screen('Flip', w);
    while resp~=p.space
        [resp, timeStamp] = checkForResp(p.space, p.escape);
        if resp==-1; ListenChar(0); return; end;
    end
    
    % stop recording and save
    if p.eyeTracking
        try
            write_parallel(0); % close file, end recording
        catch thisError
            disp(thisError);
            p.eyeTracking = 0;
            p.eyeTrackingStopped = sprintf('End of block %i',b);
        end
    end
    
end
%% end block loop

ListenChar(0);
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


function s = make_radial_checkerboard(p)
%make cartesian matries
[x,y] = meshgrid( linspace(-1,1,2*p.outStimRadiusPix), linspace(-1,1,2*p.outStimRadiusPix));
%[x,y] = meshgrid( linspace(-p.outStimRadiusDeg,p.outStimRadiusDeg,2*p.outStimRadiusDeg), linspace(-p.outStimRadiusDeg,p.outStimRadiusDeg,2*p.outStimRadiusDeg));

%make polar angle matrices
r = sqrt(x.^2+y.^2);
a = atan2(y,x)/(2*pi)+.5;  % 0<=a<=1

%make the checkerboard
radCheck = sign(sin(p.nRad*pi*r).* sin(p.nAng*pi*a));

%make it round
radCheck(r>1)  =0;

s{1} = (radCheck+1)*127+1;
s{2} = ((radCheck*-1) + 1)*127 + 1;

return

