function [x, y] = Hammer_clean(x,y)
%{
Connor Marrs
October 25, 2018
Cleaning up Noise in Hammer Video 
%}

if isempty(x); return; end

%get rid of single point outliers

for n=1:length(x)
    dx = abs(x(n)-x);
    dx(n) = nan;
    mindx(n) = min(dx);
    dy = abs(y(n)-y);
    dy(n) = nan;
    mindy(n) = min(dy);
end
x = x(mindx<=2 & mindy<=2);
y = y(mindx<=2 & mindy<=2);

%get rid of groups far from center
if isempty(x); return; end
for n=1:length(x)
    dis = sqrt((x(n)-x).^2 + (y(n)-y).^2);
    meandis(n) = mean(dis(~isnan(dis)));
end
Cdis = mean(meandis);
Cstd = std(meandis);
ind = meandis<Cdis+2*Cstd;
x = x(ind);
y = y(ind);