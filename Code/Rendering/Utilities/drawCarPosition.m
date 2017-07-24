function pos = drawCarPosition(city)
% A helper function that defines the allowable car positions for
% each city.
%
% Copyright, Henryk Blasinski 2017

%% Three possible cities.  A name for each would be good

switch city
    case 1
        xRange = {[-25000, -3500],[-3499,3500],[3501 25000]};
        yRange = {[-3500/2,2500/2],[-3500/2,2500/2],[-3500/2,2500/2]};
    case 2
        xRange = {[-14450-3500/2,-14450+3500/2],[-10950-3500/2 28549],[ -28784 -17950+3500/2]};
        yRange = {[-32299, 28619],[-14385-3500/2 -14385+3500/2],[7119+3500/2 7119-3500/2]};
    case 3
        xRange = {[-14450-3500/2,-14450+3500/2]};
        yRange = {[-35580 35785]};
    case 4
        xRange = {[-115-3500/2,-115+3500/2],[-35948 -3615 + 3500/2]};
        yRange = {[-35580 35785],[-47-3500/2 -47+3500/2]};
end


%% Not sure what these are

xIntervals = cellfun(@(x) abs(x(1)-x(2)),xRange);
yIntervals = cellfun(@(x) abs(x(1)-x(2)),yRange);

ranges = cumsum(xIntervals/sum(xIntervals));
xRangeDraw = rand(1,1);

xRangeId = find(xRangeDraw < ranges,1);

xPos = round(rand(1,1)*xIntervals(xRangeId) + xRange{xRangeId}(1));
yPos = round(rand(1,1)*yIntervals(xRangeId) + yRange{xRangeId}(1));

pos = [xPos yPos];

end



