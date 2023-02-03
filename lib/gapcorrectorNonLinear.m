function tn = gapcorrectorNonLinear( tk,debug )
    %GAPCORRECTORNONLINEAR Detect and fill gaps in pulse series by interpolation

    % Inputs: tk = Event series [1xN]
    %         degub = Logical. Set true for gap-to-gap visual inspection (Default: false)
    %         baseline = Baseline for threshold computing [1x(N-1)]

    warning('off', 'MATLAB:interp1:NaNstrip')
    
    % Threshold multipliers for upper and lower thresholds
    kupper = 1.5;
    kupperFine = 1/kupper*1.15;
    klower = 1/kupper*0.75;
    
    if nargin < 2
        debug = false;
    end
    
    tk = tk(:);
    dtk = diff(tk);
    
    % Remove false positives
    if nargin < 3
        baseline = computeBaseline(dtk);
    end
    fp = dtk<0.7*baseline;
    tk(find(fp)+1) = [];
    tn = tk;
    dtk = diff(tk);
    dtn = dtk;
    
    % Gaps are detected by deviation from the median in difference series
    if nargin < 3
        baseline = computeBaseline(dtk);
    end   
    gaps = find(dtk>baseline*kupper & dtk>0.5);
    if isempty(gaps), return; end
    thresholdAtGap = baseline(gaps)*kupper;
    
    % Gaps on first and last pulses are not allowed
    while gaps(1)<2
        tn(1) = [];
        dtk(1) = [];
        baseline(1) = [];
        gaps = find(dtk>baseline*kupper);
        thresholdAtGap = baseline(gaps)*kupper;
        if isempty(gaps), return; end
    end
    while gaps(end)>numel(dtk)-1
        tn(end) = [];
        dtk(end) = [];
        baseline(end) = [];
        gaps = find(dtk>baseline*kupper);
        thresholdAtGap = baseline(gaps)*kupper;
        if isempty(gaps), return; end
    end
     
    if debug
        f = set(gcf, 'Position', get(0, 'Screensize'));
    end
    
    nfill = 1; % Start filling with one sample
    while ~isempty(gaps)
        % In each iteration, try to fill with one more sample
        for kk = 1:length(gaps)
            if kk==1 && debug
                subplot(211);
                hold off
                stem(dtn); hold on;
                stem(gaps,dtn(gaps),'r');
                hold on
                plot(baseline*kupper,'k--')
                axis tight
                ylabel('Original RR [s]')
            end
            
            auxtn = nfillgap(tn,gaps,gaps(kk),nfill);
            auxdtn = diff(auxtn);
            
            correct = auxdtn(gaps(kk):gaps(kk)+nfill)<kupperFine*thresholdAtGap(kk);
            limitExceeded = auxdtn(gaps(kk):gaps(kk)+nfill)<klower*thresholdAtGap(kk) | ...
                auxdtn(gaps(kk):gaps(kk)+nfill)<0.5;
            
            if debug
                if limitExceeded
                    debugplots(auxdtn,gaps(kk),kupperFine*thresholdAtGap(kk),klower*thresholdAtGap(kk),nfill,false);
                else
                    debugplots(auxdtn,gaps(kk),kupperFine*thresholdAtGap(kk),klower*thresholdAtGap(kk),nfill,correct);
                end
            end
            
            if limitExceeded
                % Check that lower theshold is not exceeded. Use previous nfill instead
                auxtn = nfillgap(tn,gaps,gaps(kk),nfill-1);
                auxdtn = diff(auxtn);
                if debug
                    debugplots(auxdtn,gaps(kk),kupperFine*thresholdAtGap(kk),klower*thresholdAtGap(kk),nfill-1,true);
                end
                tn = auxtn;
                gaps = gaps+nfill-1;
            elseif correct
                % If correct number of samples, save serie
                tn = auxtn;
                gaps = gaps+nfill;
            end
        end
        
        % Compute gaps for next iteration
        dtn = diff(tn);
        if nargin < 3
            baseline = computeBaseline(dtn);
        end
        gaps = find(dtn>baseline*kupper & dtn>0.5);
        thresholdAtGap = baseline(gaps)*kupper;
        nfill = nfill+1;
    end
    
    if debug
        close(f);
    end

end

function tn = nfillgap(tk,gaps,currentGap,nfill)
    dtk = diff(tk);
    gaps(gaps==currentGap) = [];
    dtk(gaps) = nan;
    gap = dtk(currentGap);
    previousIntervals = dtk(max(1,currentGap-20):currentGap-1);
    posteriorIntervals = dtk(currentGap+1:min(end,currentGap+20));
    npre = numel(previousIntervals);
    npos = numel(posteriorIntervals);
    intervals = interp1([1:npre nfill+npre+2:nfill+npre+npos+1],[previousIntervals; posteriorIntervals],...
        npre+1:npre+nfill+1,'pchip');
    intervals = intervals(1:end-1)*gap/(sum(intervals,'omitnan')); % map intervals to gap
    tn = [tk(1:currentGap); tk(currentGap)+cumsum(intervals)'; tk(currentGap+1:end)];  
end

function baseline = computeBaseline(rr)
    wind = 30;
    if length(rr)<wind, wind = length(rr); end
    mf = medfilt1([flipud(rr(1:wind/2)); rr; flipud(rr(end-wind/2+1:end))],wind-1);
    mf(mf>1.5) = 1.5;
    baseline = mf(wind/2+1:end-wind/2);
end

function debugplots(dtn,gap,upperThreshold,lowerThreshold,nfill,correct)
    subplot(212); hold off;
    stem(dtn); hold on;
    if correct
        stem(gap:gap+nfill,dtn(gap:gap+nfill),'g','LineWidth',1);
    else
        stem(gap:gap+nfill,dtn(gap:gap+nfill),'r','LineWidth',1);
    end
    xlim([max(0,gap-50) min(gap+50,numel(dtn))])
    ylabel('Corrected RR [s]')
    xlabel('Samples');% ylim([0 1.7])
    line(xlim,[upperThreshold upperThreshold],'Color','k');
    line(xlim,[lowerThreshold lowerThreshold],'Color','k');
    pause;
end