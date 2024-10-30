function [h,s,p,edges] = pretty_hist(x,varargin)
  % pretty_hist(x,color,nBins,normalization,doFit)
  % creats pretty histograms by default, kinda inspired by nHist function
  % but a lot less hassel and lot more Joe!
  % plots histograms without edges and uses stair plot to get
  % envelope which is just prettier...
  % color need to be either a 1x3 or 2x3 rgb value,
  % when color is a 2x3, then first row is the line color and second is the face
  % 
  % [h,s,p] = pretty_hist(x,varargin)
  % [h,s,p] = pretty_hist(x,color,nBins,normalization,doFit)
  % [h,s,p] = pretty_hist(x,[1 0 0],100,'probability',false,);

  % h - handle to histogram
  % s - handle to stairs
  % p - handle to fit plot

  % check if we already have an open figure with a valid axis and get hold state
  fig = get(groot,'CurrentFigure');
  if ~isempty(fig) 
    wasHold = ishold;
  end

  alpha = 0.15; % set default alpha, can be overwritten using a color trick...
  nVarargs = length(varargin);
  defaultFaceColor = [0.0, 0.6, 0.0];

  if (nVarargs > 0) && ~isempty(varargin{1}) % color argument was provided
    if ~(size(varargin{1},1)>1) %same color for line and face
      lineColor = varargin{1};
      faceColor = lineColor;
    else
      lineColor = varargin{1}(1,:); % first row line Color
      faceColor = varargin{1}(2,:); % second row line Color
      if faceColor(1)<0
        alpha = 0;
        faceColor = lineColor; % won't be used but needs to be a real color value
      end
    end
  elseif (nVarargs > 0) && isempty(varargin{1}) % color arg was emptpy
    faceColor = defaultFaceColor;
    lineColor = faceColor;
  end

  if(nVarargs == 0) % only default parameters
    faceColor = defaultFaceColor; % nice dark green
    lineColor = faceColor;
    nBins = 50;
    normalization = 'probability';
    doFit = true;
    edges = [];
  elseif (nVarargs == 1)
    nBins = 50;
    normalization = 'probability';
    doFit = true;
    edges = [];
  elseif (nVarargs == 2)
    nBins = varargin{2};
    normalization = 'probability';
    doFit = true;
    edges = [];
  elseif (nVarargs == 3)
    nBins = varargin{2};
    normalization = varargin{3};
    doFit = true;
    edges = [];
  elseif (nVarargs == 4)
    nBins = varargin{2};
    normalization = varargin{3};
    doFit = varargin{4};
    edges = [];
  elseif (nVarargs == 5)
    nBins = varargin{2};
    normalization = varargin{3};
    doFit = varargin{4};
    edges = varargin{5};
  else
    short_warn('Unknown pretty_hist parameters. Trying to use default!');
    faceColor = [0.0, 0.6, 0.0]; % nice dark green
    nBins = 50;
    normalization = 'probability';
    doFit = true;
  end

  if isempty(edges)
    [counts,edges] = histcounts(x,nBins,'Normalization',normalization);
  else
    [counts,edges] = histcounts(x,edges,'Normalization',normalization);
  end 

  % make y axis nicer for Probability plots...
  if strcmp(normalization,'probability')
    counts = counts*100; % show prob. in %
  end
  % in some case matlab creates empty bins, but we no like empty bins!
  % nonZeroBins = length(nonzeros(counts));
  % if nonZeroBins < nBins
  %   nBins = nonZeroBins;
  %   [counts,edges] = histcounts(x,nBins,'Normalization',normalization);
  % end
  h = histogram('BinEdges',edges,'BinCounts',counts); % normal hist
  h.FaceColor = faceColor;
  h.FaceAlpha = alpha;
  h.EdgeColor = 'none';
  ax = gca; % get current axes handle
  
  % now we have a histgram without edges
  % we add the edge on top, because it looks nice...
  hold on;
  fullEdges = edges;
  edges = edges(1:end-1);
  s = stairs(edges,counts); % stairs on top
  s.LineWidth = 1;
  s.Color = lineColor;

  if strcmp(normalization,'probability')
    ylabel(normalization);
  end

  % ger kernel estimate for nBins points 
  isNoNaN = sum(isfinite(x));
  if doFit && isNoNaN
    xEval = linspace(fullEdges(1),fullEdges(end),nBins);
    % [est,pointLoc] = ksdensity(x,xEval,'Bandwidth',5,'BoundaryCorrection','reflection'); 
    [est,pointLoc] = ksdensity(x,xEval,'BoundaryCorrection','reflection'); 
    est = est./max(est(:)).*max(counts(:));
    p = plot(pointLoc,est);
    p.LineWidth = 2;
    p.Color = lineColor;
  elseif doFit && ~isNoNaN % prepare plot anyway for later updates with data?
    xTemp = zeros(1,numel(x));
    xEval = linspace(fullEdges(1),fullEdges(end),nBins);
    [~,pointLoc] = ksdensity(xTemp,xEval,'BoundaryCorrection','reflection'); 
    p = plot(pointLoc,xTemp);
    p.LineWidth = 2;
    p.Color = lineColor;
    p = []; 
  end

  axis tight;
  grid on;
  if wasHold
    hold on;
  else
    hold off;
  end
end