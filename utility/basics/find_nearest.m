function [nearestVal,nearestIdx] = find_nearest(findVals,data)
  % [nearestVal,nearestIdx] = find_nearest(findVals,data)
  % looks for nearerst val and index for findVals on data
  findVals = double(findVals);
  data = double(data);
  % check for data comformity
  if ~isvector(data)
    error('ndims(data) > 2! data must be a vector!')
  end
  if ~isvector(findVals)
    error('ndims(findVals) > 2! findVals must be a scalar or a vector!')
  end

  % handle scalar and vector inputs for findVals
  if numel(findVals) > 1
    nearestIdx = zeros(numel(findVals),1);
    nearestVal = zeros(numel(findVals),1);
    for iVal = 1:numel(findVals)
      %FIXME is there a way to vectorize this?
      [~,nearestIdx(iVal)] = min(abs(data-findVals(iVal)));
      nearestVal(iVal) = data(nearestIdx(iVal));
    end
  else % return scalar answer
    difference = abs(data-findVals);
    [~,nearestIdx] = min(difference); %index of closest value
    nearestVal = data(nearestIdx); %closest value
  end
end