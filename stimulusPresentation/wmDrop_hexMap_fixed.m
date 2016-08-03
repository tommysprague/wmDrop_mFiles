% 295 s - 131 TRs @ 2.25 sec/TR (per attn/spatial condition)

function wmDrop_hexMap_fixed 

% wmDrop Localizer - adapted from map2dAttn_flicker4.m
%
% instead of grid of mapping positions, now present stimuli along a
% hexagonal grid which is offset a different amount each run
%


% NOTE: this is the "corrected" version of the mapping task, where accuracy
% is computed correctly (see Supplemental Experimental Procedures). Only
% change from wmDrop_hexMap.m is to fix accuracy computation.
% wmDrop_hexMap.m was used for most subjects in Sprague, Ester & Serences,
% 2016 (BC was run w/ this script, but accuracy was held to ~89% to
% maintain consistency w/ previous subjects)
e

%% get user input about sub name etc...

% mygreen = [0 165 0];
% myred   = [255 0 0];

myred = [255 0 0];
myblue = [50 50 255];
mypurple = [135 0 135];

mywhite = [255 255 255];
myblack = [0 0 0];
mygray = [65 65 65]; % fix this one



warning('off','MATLAB:dispatcher:InexactMatch');
%Screen('Preference', 'SkipSyncTests', 1)
%get subject info
prompt = {'Subject Name','Ang offset? (degrees, -30 to 30)', 'Block number', 'Random seed', 'fMRI','Eye tracking','WM Distance (dva)'};
%grab a random number and seed the generator (include this as a gui field
%in case we want to repeat an exact sequence)
s = round(sum(100*clock));
%put in some default answers
defAns = {'XX','XXX','1',num2str(s),'1','0','0.6'};
box = inputdlg(prompt,'Enter Subject Information...', 1, defAns);

p.exptName = 'wmDrop_hexMap';

if length(box)==length(defAns)
    p.subName=char(box{1});
    p.sessionNum=3;%str2num(box{2});
    %p.cond = str2num(box{2});
    p.angOffset = str2num(box{2});
    p.nBlocks=eval(box{3});
    p.rndSeed=str2num(box{4});
    p.fMRI=str2num(box{5});
    p.eyeTracking = str2num(box{6});
    p.wmDifficultyDeg = str2num(box{7}); % WM difficulty is how far, in radians from center of stim, targ & probe are apart (multiples of pi/2)
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
p.windowed = ~p.fMRI;                     % if 1 then the display will be in a small window, use this for debugging.

% monitor stuff
p.refreshRate = 60;                 % refresh rate is normally 100 but change this to 60 when on laptop!
if p.fMRI
    p.vDistCM = 375; % CM
    p.screenWidthCM = 120; % CM
else % behavioral room w/ eyetracker
    p.vDistCM = 62;
    p.screenWidthCM = 51;
end

% ideally, let's do 2 hex steps = wmDrop's target eccentricity (max ecc of
% target center = 3*steps, then stimRad on top oft that

% maybe 35 if we drop center point, which we should...
p.nLoc = 36; % instead of grid, let's just do number of positions, then have x, y coords for this number of positions

%stimulus geometry (in degrees, note that these vars have the phrase 'Deg'
%in them, used by pix2deg below)

p.wmEcc = 3.5; % use this to construct geometry?
p.stepSize = 1.625;% p.wmEcc/2;

% for now, let's do rad = stepsize
p.radDeg = 3.25/3;%p.wmEcc/3;

%p.radDeg = p.usedScreenSizeDeg/(p.nLoc+2);
%p.sfDeg  = 1.5 * .4532;
p.sfDeg = .6785;  % TODO: check this out...

p.targetProbeMaxRadDeg = p.radDeg/2; % ...and this or SMALLER TCS 11/13/12

p.wmProbeWindowAng = 20; % plus/minus 20 degrees from target in offset direction (left/right, up/down)

p.probeCueSizeDeg = 0.5; % the bar that indicates the judgment to make
p.wmRespCueColor = mygray;

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

% TODO: also allow left/right arrow...
p.keys = [KbName('b'),KbName('y')]; % 1 and 2 keys
p.space = KbName('space');
p.start = KbName('t');

p.backColor     = [128, 128, 128];      % background color

%fixation point properties
p.fixColor      = [180 180 180];%[155, 155, 155];
% p.fixSizeDeg    = .25 * .4523;                  % size of dot
p.fixSizeDeg    = .2;
p.appSizeDeg    = p.fixSizeDeg * 4;     % size of apperture surrounding dot

% target/probe properties (WM)
% target/probe properties (WM) ( both 0.5 - now 0.5 * .4523)
p.targetSizeDeg = 0.2;
p.probeSizeDeg = 0.2;

% TODO - get this sorted out a bit like wmDrop.m
p.wmTargetColor = myred;
p.wmProbeColors  = [myblack];

% target properties (fix/checks)
p.stimContrast = 1;
%p.contrastChange = 0.5; % will change by this*100% difference from bg color
%p.targetContrast = p.stimContrast - p.stimContrastChange;
%p.targetFixColor = p.backColor + (abs(p.fixColor - p.backColor)) * (1-p.fixContrastChange);

% trial setup
p.percentNull = 0.2;
p.repetitions = 1; 
if p.windowed % DEBUG MODE
    p.repetitions = 1;
end
p.nValidTrials = p.repetitions * p.nLoc;    % length p.cond, for now, should always be 1
p.numNullTrials = 6;%((1/(1-p.percentNull)) - 1) * p.nTrials; % TODO: check this...
p.nTrials = p.nValidTrials + p.numNullTrials; % no stim

% stimulus timing (in seconds - convert to vid frames, if necessary, later)
p.stimExposeDur = 3; % keep this constant
p.wmTargetDur   = 0.5; % s
p.wmProbeDur    = 0.75; % s
p.flickerFreq   = 6; % Hz
p.flickerPer  = 1/p.flickerFreq; % s
%p.contrastChangeWindow = [0.5 p.stimExposeDur-0.75]; % period over which contrast can change
%p.minTargFrame = p.contrastChangeWindow(1)*p.flickerFreq + 1; % these are period index (which period can we start target?)
%p.maxTargFrame = p.contrastChangeWindow(2)*p.flickerFreq + 1; % +1 is because want cycle *starting at* this point in time
p.responseWindow = 1.0;
%p.nTargs = 1;
%p.minTargSep = 1; % number of periods
% time between stimulus presentations, constant

p.ITIs = linspace(2,4.5,p.nTrials)';
p.ITIs = p.ITIs(randperm(length(p.ITIs)));
%p.ITI           = 2.0 * ones(p.nTrials,1);
p.trialDur      = p.wmTargetDur + p.stimExposeDur + p.wmProbeDur + p.ITIs;
p.nTRsWait      = 0;
p.startWait     = 2;
p.passiveDur    = 10; % in sec, how long to wait at end of a block (and beginning?)
p.TR            = 2.25;
p.expDur        = sum(p.trialDur) + p.passiveDur + p.startWait + p.nTRsWait*p.TR;

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

% make sure the refreshrate is ok
if abs(p.fps-p.refreshRate)>5
    Screen('CloseAll');
    disp('CHANGE YOUR REFRESH RATE')
    ListenChar(0);
    clear all;
    return;
end

%% HEXAGONAL GRID - IMPORTANT! for now, just manually generating

% 35 points (don't include 0,0)
xgrid = [-1.5 -0.5 0.5 1.5  -2 -1 0 1 2  -2.5 -1.5 -0.5 0.5 1.5 2.5  -3 -2 -1 1 2 3  -2.5 -1.5 -0.5 0.5 1.5 2.5  -2 -1 0 1 2  -1.5 -0.5 0.5 1.5];
ygrid = [-3 -3 -3 -3   -2 -2 -2 -2 -2  -1 -1 -1 -1 -1 -1  0 0 0 0 0 0  1 1 1 1 1 1  2 2 2 2 2  3 3 3 3]*(sqrt(3)/2);

% cartesian coords

[thgrid, rgrid] = cart2pol(xgrid,ygrid); % CHECK THIS FOR up/down
thgrid = thgrid + deg2rad(p.angOffset);
[xgrid_adj, ygrid_adj] = pol2cart(thgrid,rgrid);

% SCREEN COORDINATES!!!!

p.xGridDeg =  xgrid_adj*p.stepSize;
p.yGridDeg = -ygrid_adj*p.stepSize;

clear xgrid ygrid thgrid rgrid xgrid_adj ygrid_adj;

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

    
    fName = sprintf('%s/data/%s_%s_sess%02.f_run%02.f.mat',p.root, p.subName,p.exptName,p.sessionNum,b);
    %fName=[p.root, '/Subject Data/', p.exptName, num2str(p.cond),'_', p.subNameAndDate, '_sess', num2str(p.sessionNum), '_blkNum', num2str(b), '.mat'];
    if exist(fName,'file')
        Screen('CloseAll');
        msgbox('File name already exists, please specify another', 'modal')
        ListenChar(0);
        return;
    end
    
    % when eyetracking updated, use: p.et_fname = sprintf('%s/%s/%s_%s_sess%02.f_run%01.f_%s.idf',p.exptName,p.subName, p.subName,p.exptName,p.sessionNum,p.runNum, datestr(now,30));
    
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
 
    
    p.stimLocs = (1:length(p.xGridDeg))';
    p.stimLocsDeg = [p.xGridDeg(p.stimLocs)' p.yGridDeg(p.stimLocs)'];    
    p.stimLocsPix = p.stimLocsDeg * p.ppd;

    % mark null trials
    p.stimLocs(end+1:p.nTrials) = NaN;
    %p.stimLocsY(end+1:p.nTrials) = NaN;
    p.null = zeros(p.nTrials,1); p.null(isnan(p.stimLocs)) = 1;
    
    
    p.stimLocsDeg(isnan(p.stimLocs),:) = NaN;
    p.stimLocsPix(isnan(p.stimLocs),:) = NaN;
    

    % uses same behavioral scheme as wmDrop - and will adjust accordingly
    % if we need to do up/down left/right task
    nRealTrials = p.nTrials - sum(p.null);
    %p.probeColorMapping(:,1) = 1+[zeros(floor(nRealTrials/2),1); ones(ceil(nRealTrials/2),1)]; % if floor/ceil used reliably together, should always come out same... ( sum = p.nTrials)
    %p.probeColorMapping(:,2) = 1+[ones(ceil(nRealTrials/2),1); zeros(floor(nRealTrials/2),1); ]; % if floor/ceil used reliably together, should always come out same... ( sum = p.nTrials)
    
    
    % let's just do these as x, y rather than r, theta
    %tang = rand(p.nTrials,1); trad = rand(p.nTrials,1)*2*pi;
    trad = rand(nRealTrials,1); tang = rand(nRealTrials,1)*2*pi;
    p.targetLoc =  (repmat(trad*p.targetProbeMaxRadPix,1,2) .* [cos( tang ) sin( tang ) ]);% + p.stimLocsPix; % pick radii
    clear tang trad;
    %p.targetLoc(:,2) = [linspace(0,2*pi,floor(p.nTrials/2))'; linspace(0,2*pi,ceil(p.nTrials/2))'];
    
    
    
    
    %% FIRST, shuffle taret/probe positions for all real trials

    p.rndIndWM = randperm(nRealTrials); % different random indices for working memory targets (so no wm/location correspondence)
    %p.probeColorMapping = p.probeColorMapping(p.rndIndWM,:);
    p.targetLoc = p.targetLoc(p.rndIndWM,:);  % both in pixels
    %p.probeLoc  = p.probeLoc(p.rndIndWM,:);    
    
    % now, add p.stimLocsPix, p.targLoc
    p.targetLoc = p.targetLoc + p.stimLocsPix(1:nRealTrials,:);
    %p.probeLoc  = p.targetLoc + p.probeLoc;

    % now, add in the null trials
    p.targetLoc(end+1:p.nTrials,:) = NaN;
%    p.probeLoc(end+1:p.nTrials,:) = NaN;
    %p.probeColorMapping(end+1:p.nTrials,:) = NaN;
    
    %% SHUFFLE trial order (EVERYTHING, including WM stuff)
    
    p.rndInd = randperm(p.nValidTrials);
    tmp1 = randperm(p.nValidTrials-2); tmp2 = 2:p.nValidTrials-1;
    null_after = sort(tmp2(tmp1(1:p.numNullTrials)));
    null_after = null_after + (0:(length(null_after)-1));
    for ii = 1:length(null_after)
        p.rndInd = [p.rndInd(1:null_after(ii)) p.nValidTrials+ii p.rndInd((null_after(ii)+1):end)];
    end
    
    
    %p.rndInd = randperm(p.nTrials);
    p.stimLocs = p.stimLocs(p.rndInd);
    p.stimLocsDeg = p.stimLocsDeg(p.rndInd,:);
    p.stimLocsPix = p.stimLocsPix(p.rndInd,:);
    
    p.null = p.null(p.rndInd);
    
    p.targetLoc = p.targetLoc(p.rndInd,:);
%    p.probeLoc  = p.probeLoc(p.rndInd,:);
   % p.probeColorMapping = p.probeColorMapping(p.rndInd,:);
    
 
    % this should be true now
%    p.targetLoc(p.null,:) = NaN;
%    p.probeLoc(p.null,:)  = NaN;
%    p.probeColorMapping(p.null,:) = NaN;
    
%    p.targProbeVector = p.probeLoc - p.targetLoc;
        
    %allocate some arrays for storing the subject response
    p.responseTime =          nan(1, p.nTrials);       % store the rt on each trial
    %p.resp =        zeros(p.nTrials, p.stimExpose);     % store the response
    p.wmResp =      zeros(p.nTrials,1);
    p.wmCorrect =   0;
    
    p.corrResp = round(rand(p.nTrials,1))+1;
    p.respCoord = round(rand(p.nTrials,1)) + 1;
    
    p.corrResp(p.null==1) = NaN;
    p.respCoord(p.null==1) = NaN;
    
    p.trialStart =  nan(p.nTrials,1);    % equal to targetStart
    p.trialEnd =  nan(p.nTrials,1);
    p.stimStart =  nan(p.nTrials,1);
    p.stimEnd =  nan(p.nTrials,1);
    
    % generate checkerboards we use...
    c = make_checkerboard(p.radPix,p.sfPix,p.stimContrast);
    stim(1)=Screen('MakeTexture', w, c{1});
    stim(2)=Screen('MakeTexture', w, p.backColor(1)*ones(size(c{2})));
    stim(3)=Screen('MakeTexture', w, c{2});
    
    % pick the stimulus sequence for every trial (the exact grating to be shown)
    for i=1:p.nTrials
        %p.flickerSequ = repmat([ones(1,p.flickerFrames/2) 2*ones(1,p.flickerFrames/2)],1,p.stimExpose/p.flickerFrames);
        %p.stimSequ(i,:)=p.flickerSequ;
        p.flickerSequ = repmat([ones(1,p.flickerFrames/2) 2*ones(1,p.flickerFrames/2) 3*ones(1,p.flickerFrames/2) 2*ones(1,p.flickerFrames/2)],1,(0.5)*p.stimExpose/p.flickerFrames);
        %p.stimDimSequ(1,:) = zeros(1,size(p.flickerSequ,2));
        % mark the tarket spots with low contrast stims
    end
    
    % use the stim sequences generated above, shuffled wrt trial order, for
    % fix dimming sequences
%    p.fixDimSequ = p.stimDimSequ; % 1 if fix dims, 0 if not
    %p.rndInd2 = randperm(size(p.fixDimSequ,1));
%    p.fixDimSequ = p.fixDimSequ(p.rndInd2,:);
    
    % put up a message to wait for a space bar press
%    elseif p.cond == 3
    textMotion = 'Compare probe position to remembered position(s)';
%        p.validTarget = zeros(size(p.dimFix));
    
    tCenterMo = [center(1)-RectWidth(Screen('TextBounds', w, textMotion))/2 center(2)/2];
    
    Screen('DrawText', w, textMotion, tCenterMo(1), tCenterMo(2), p.textColor);
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
    %after all initialization is done, sit and wait for scanner synch (or
    %space bar)
    resp=0;
%     disp('checking for response');
%     while 1
%         [resp, timeStamp] = checkForResp([p.space,p.start],p.escape);
%         if resp==p.space || resp==p.start || resp==p.escape
%             break;
%         end
%         if resp==-1; ListenChar(0); return; end;
%     end
    
    FlushEvents;
    GetChar;
    
    disp('starting block');
    cumTime = GetSecs;
    
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
    
    %waited = 1; % we've already waited for 1 TR here...
    

    FlushEvents;
    
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    cumTime = GetSecs;
    p.startExp = GetSecs;%cumTime;
    
    while GetSecs < p.startExp + p.startWait
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
    
    cumTime = GetSecs;
    
    rTcnt = 1;
    
    save(fName,'p');
    
    %% here is the start of the trial loop
    for t=1:p.nTrials

        xLoc = p.stimLocsPix(t,1)   + center(1);
        yLoc = p.stimLocsPix(t,2)   + center(2);
        
        targXLoc = p.targetLoc(t,1);% + center(1);
        targYLoc = p.targetLoc(t,2);% + center(2);

        if ~p.null(t)
            pang = deg2rad(rand(1) * 2 * p.wmProbeWindowAng - p.wmProbeWindowAng);
            % in pixels, (x,y) coords
            prad = p.wmDifficultyPix;
            
            if p.respCoord(t) == 1 % LEFT/RIGHT
                if p.corrResp(t) == 1 % LEFT
                    pang = pang+pi;
                end
                
            elseif p.respCoord(t) == 2 % UP/DOWN
                pang = pang + pi/2;
                if p.corrResp(t) == 2 % DOWN
                    pang = pang+pi;
                end
            end
            
            p.probeLoc(t,:) =  (repmat(prad,1,2) .* [cos(pang) -sin(pang)]) + p.targetLoc(t,:);% + p.targetLoc;
            clear pang prad;
            
            
            probeXLoc = p.probeLoc(t,1);% + center(1);
            probeYLoc = p.probeLoc(t,2);% + center(2);
        else
            probeXLoc = NaN;
            probeYLoc = NaN;
            p.probeLoc(t,:) = [NaN NaN];
        end
        
        if ~p.fMRI
            disp(sprintf('dim fix = %i, dim stim = %i',p.dimFix(t),p.dimStim(t)));
        end
        
        p.trialStart(t) = GetSecs;
        
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
        
        % now just use CenterRectOnPoint
        stimRect = CenterRectOnPoint([0 0 2*p.radPix 2*p.radPix],xLoc, yLoc);
        

        
        while GetSecs <= p.trialStart(t) + p.wmTargetDur
            [resp, timeStamp] = checkForResp(p.escape, p.escape);
            if resp==-1; ListenChar(0); return; end;
            
            % TARGET
            if ~p.null(t)
                Screen('DrawDots',w,[targXLoc targYLoc],p.targetSizePix,p.wmTargetColor,center,1);
            end
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0);
            Screen('DrawingFinished',w);
            Screen('Flip',w);
        end
        
        frmCnt=1;
        p.stimStart(t) = GetSecs;   % start a clock to get the stim onset time
        
        % STIMULUS
        while frmCnt<=p.stimExpose % if we want multiple exposure durations, add that here
                       
            if ~p.null(t)
                
                % checkerboard
                Screen('DrawTexture',w,stim(p.flickerSequ(1,frmCnt)),Screen('Rect',stim(p.flickerSequ(1,frmCnt))),stimRect);
                
                % apperture around fixation
                Screen('DrawDots', w, [0 0], p.appSizePix, p.backColor, center, 1);
                
%                 if resp && find(p.keys==resp)
%                     p.resp(t, frmCnt) = find(p.keys==resp);
%                 end
                % TODO - gather responses longer...

            end
            
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0);
            Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            Screen('Flip', w);
            
            % check response...
            [resp, timeStamp] = checkForResp(p.keys, p.escape); % checks both buttons...
            
            if resp==-1; ListenChar(0); return; end;

            frmCnt = frmCnt + 1;
        end
        
        p.stimEnd(t) = GetSecs;

        if p.respCoord(t) == 1 % horizontal bar
            bar_coord = [-1 1; 0 0]*p.probeCueSizePix;
            
        elseif p.respCoord(t) == 2 % vertical bar
            bar_coord = [0 0; -1 1]*p.probeCueSizePix;
        end
        
        
        % PROBES
        if ~p.null(t)
            %Screen('DrawDots',w,[targXLoc  targYLoc] ,p.probeSizePix,p.wmProbeColors(p.probeColorMapping(t,1),:),center,1);
            Screen('DrawDots',w,[probeXLoc probeYLoc],p.probeSizePix,p.wmProbeColors(1,:),center,1);
            Screen('DrawLines',w,bar_coord,3,p.wmRespCueColor,center);
        end
        

        
        
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); 
        Screen('DrawingFinished',w);
        Screen('Flip',w);
        
        while GetSecs <= cumTime + p.stimExposeDur + p.wmTargetDur + p.wmProbeDur
            [resp, timeStamp] = checkForResp(p.keys, p.escape);
            if resp==-1; ListenChar(0); return; end;
            %if p.cond == 3
                if resp && find(p.keys==resp)
                    p.wmResp(t) = find(p.keys==resp); % (most recent response is that recorded)
                end
            %end
        end
        
        % clear out screen
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %draw fixation point
        Screen('Flip',w);
        
        % write end of eyetracking event (give time for return to fix)
        if p.eyeTracking        
            try
                write_parallel(255); % event ID will be trial number for now
                disp(sprintf('ITI of trial %i',t));
                %clear msg out;
            catch thisError
                disp(thisError);
                p.eyeTracking = 0;
                p.eyeTrackingStopped = sprintf('ITI of trial %i (%i)',t,255);
                % attempt to save
            end
        end
        
        while GetSecs <= cumTime + p.trialDur(t)
            [resp, timeStamp] = checkForResp(p.keys, p.escape);
            if resp==-1; ListenChar(0); return; end;
            %if p.cond == 3
                if resp && find(p.keys==resp)
                    p.wmResp(t) = find(p.keys==resp); % (most recent response is that recorded)
                end
            %end
        end
            
        if p.null(t)
            disp(sprintf('Trial %03i: Null',t));
            p.correct(t) = NaN; % added TCS 1/28/2016, to compute acc correctly
        else 
            %mystr = {'no match','match'};
            if p.corrResp(t) == p.wmResp(t)
                p.correct(t) = 1;
            %if p.targetProbeMatch(t) && p.wmResp(t) == 1 || ~p.targetProbeMatch(t) && p.wmResp(t) == 2
                %disp(sprintf('Trial %03i: Correct (%s)',t,mystr{p.targetProbeMatch(t)+1}));
                %p.wmCorrect = p.wmCorrect + 1;
            else
                p.correct(t) = 0;
                %disp(sprintf('Trial %03i: Incorrect (%s)',t,mystr{p.targetProbeMatch(t)+1}));
                %p.wmIncorrect = p.wmIncorrect + 1;
            end
            clear mystr
        end
        cumTime = cumTime + p.trialDur(t);
        rTcnt = rTcnt + 1;  % increment real trial counter.

        p.trialEnd(t) = GetSecs;
        
        save(fName, 'p');
        
    end
    %% end trial loop

    %   10s passive fixation at end of block
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    while GetSecs <= cumTime + p.passiveDur
        [resp, timeStamp] = checkForResp(p.escape, p.escape);
        if resp==-1; ListenChar(0); return; end;
    end
    
    p.endExp = GetSecs;
    
    p.accuracy = 100*nanmean(p.correct);% p.wmCorrect/(p.wmCorrect+p.wmIncorrect)*100;
    
    %save trial data from this block
    save(fName, 'p');
   
    str1 = sprintf('Block %i complete',b);
    tCenter1 = [center(1)-RectWidth(Screen('TextBounds', w, str1))/2 center(2)/2];
    str2 = sprintf('Accuracy: %.02f%%',p.accuracy);
    tCenter2 = [center(1)-RectWidth(Screen('TextBounds', w, str2))/2 center(2)/2];

    Screen('DrawText', w, str1, tCenter1(1), tCenter1(2)-100, p.textColor);
    Screen('DrawText', w, str2, tCenter2(1), tCenter2(2), p.textColor);

    %if p.cond ~= 3
    %    str3 = sprintf('Average response time: %4i ms',round(p.meanCorrectRT*1000));
    %    tCenter3 = [center(1)-RectWidth(Screen('TextBounds', w, str3))/2 center(2)/2];
    %    Screen('DrawText', w, str3, tCenter3(1), tCenter3(2)+100, p.textColor);
    %end
    
    
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
