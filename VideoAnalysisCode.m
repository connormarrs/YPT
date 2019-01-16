%Connor Marrs
%RCDS YPT 2018
%October 2018
%Hammer Analysis Code

clear;
pth = 'Video Analysis/';
FPS = 299.7;
dt = 1/299.7;

opt = 2;
switch opt
    case 1
        fname = [pth '1x5x10_1.MOV'];
        t0 = 28; %time to start analyzing video- determine by eye
        ycut = 1:110; NY = length(ycut); %zoom in on the relevant part of the frame. (use ginput)
        xcut = 50:275; NX = length(xcut);
    case 2
        fname = [pth '1x5x10_1.MOV'];
        t0 = 30; %time to start analyzing video- determine by eye
        ycut = 1:176; NY = length(ycut); %zoom in on the relevant part of the frame. (use ginput)
        xcut = 50:342; NX = length(xcut);
end
    
vidObj = VideoReader(fname, 'CurrentTime', t0);

if 0 %QC graph - check the x, y cuts and choose the color.
	vidFrame = readFrame(vidObj);
	vidFrame = vidFrame(ycut, xcut, :);
	figure(1); clf;
	%subplot(221); 
    H1 = pcolor(flipud(squeeze(vidFrame(:,:,1))'));  title('red');
	%subplot(222); H2 = pcolor(squeeze(vidFrame(:,:,2))');  title('green');
	%subplot(223); H3 = pcolor(squeeze(vidFrame(:,:,3))');  title('blue');
	set(H1, 'linestyle', 'none');
	set(H2, 'linestyle', 'none');
	set(H3, 'linestyle', 'none');
	return
end

NT = 300;
oldFrame = readFrame(vidObj);
oldFrame = double(squeeze(oldFrame(ycut,xcut, 1)));
%taking nxnx1 '> 2dimensional double array
[pTx, pTy, pBx, pBy] = deal(nan*ones(NT, 2));
[angTx,angTy,angBx,angBy] = deal(nan*ones(NT,1));
t=(1:NT)*dt;

%initialize graph
figure(1);clf;
subplot(121); H1 = pcolor(oldFrame'); hold on;
set(H1, 'linestyle', 'none');
set(gca, 'ydir', 'reverse');
subplot(122); H2 = pcolor(oldFrame'); hold on;
set(H2, 'linestyle', 'none');
set(gca, 'ydir', 'reverse');
HT = plot(nan,nan,'k.');
HB = plot(nan,nan,'r.');
HslpT = plot(nan,nan,'b-');
HslpB = plot(nan,nan, 'r-');

for k=1:NT
	newFrame = readFrame(vidObj);
	newFrame = double(squeeze(newFrame(ycut,xcut, 1)));
	diffFrame = newFrame - oldFrame;
	diffFrame(20:40, 130:end) = 0; %blocking rods in the video

	subplot(121); set(H1, 'cdata', newFrame');
	title(['time step ' num2str(k)]);
	subplot(122); set(H2, 'cdata', diffFrame');

	%find the areas that have changed since the last frame
	[xindT, yindT] = find(diffFrame>25);
	[xindB, yindB] = find(diffFrame<-25);

	%pick only the points along a line predicted by earlier steps
	if length(xindT)>5 %have to have enough points
		[xindT, yindT] = Hammer_clean(xindT, yindT);
		if k>5
			if all(abs(angTx(k-(1:4)))<50) %small angles - fit y(x)
				slp = pTy(k-1, 1) + (pTy(k-1, 1)-pTy(k-2,1));
				int = pTy(k-1, 2) + (pTy(k-1, 2)-pTy(k-2,2));
				yT = slp*xindT + int;
				dis = yindT - yT;
				set(HslpT, 'xdata', xindT, 'ydata', yT, 'color', 'c');
			else %large angles - fit x as a function of y
				slp = pTx(k-1, 1) + (pTx(k-1, 1)-pTx(k-2,1));
				int = pTx(k-1, 2) + (pTx(k-1, 2)-pTx(k-2,2));
				xT = slp*yindT + int;
				dis = xindT - xT;
				set(HslpT, 'xdata', xT, 'ydata', yindT, 'color', 'g');
            end
            
			xindT = xindT(abs(dis)<=2);
			yindT = yindT(abs(dis)<=2);
		end
	end

	%repeat for the other set of points
	if length(xindB)>5 %have to have enough points
		[xindB, yindB] = Hammer_clean(xindB, yindB);
		if k>10
			if all(abs(angBx(k-(1:4)))<50) %small angles - fit y(x)
				slp = pBy(k-1, 1) + (pBy(k-1, 1)-pBy(k-2,1));
				int = pBy(k-1, 2) + (pBy(k-1, 2)-pBy(k-2,2));
				yB = slp*xindB + int;
				dis = yindB - yB;
				set(HslpB, 'xdata', xindB, 'ydata', yB, 'color', 'r');
			else %large angles - fit x as a function of y
				slp = pBx(k-1, 1) + (pBx(k-1, 1)-pBx(k-2,1));
				int = pBx(k-1, 2) + (pBx(k-1, 2)-pBx(k-2,2));
				xB = slp*yindB + int;
				dis = xindB - xB;
				set(HslpB, 'xdata', xB, 'ydata', yindB, 'color', 'm');
			end
			xindB = xindB(abs(dis)<=2);
			yindB = yindB(abs(dis)<=2);
		end
	end

	%find slopes of the lines
    if length(xindT)>5
       [pTy(k,:)] = polyfit(xindT,yindT, 1);
       [pTx(k,:)] = polyfit(yindT,xindT, 1);
       angTy(k) = 180/pi*atan(pTy(k,1));
       angTx(k) = -(90+180/pi*atan(pTx(k,1)));
    end
    
    if length(xindB)>5
       [pBy(k,:)] = polyfit(xindB,yindB, 1);
       [pBx(k,:)] = polyfit(yindB,xindB, 1);
       angBy(k) = 180/pi*atan(pBy(k,1));
       angBx(k) = -(90+180/pi*atan(pBx(k,1)));
    end
    
    %update the graph
    set(HT, 'xdata', xindT, 'ydata', yindT);
    set(HB, 'xdata', xindB, 'ydata', yindB);
    pause(.2);
    oldFrame = newFrame;
end

%Can we
%use geometry?
wid = nan*ones(NT,1);
for k=1:NT
    if abs(angBy(k))<45
        x = 1:20;
        yB = pBy(k,1)*x + pBy(k,2);
        yT = pTy(k,1)*x + pTy(k,2);
        for n=1:length(x);
            tmp(n,:) = sqrt((xB(n)-x).^2 + (yB(n)-yT).^2);
        end
        wid(k) = min(tmp(:));
    end
end

%{
wid = CMsmooth(wid, 3, 'median');
widS = CMsmooth(wid, 6, 'mean');

omega = [0; diff(angTy)];
omegaS = CMsmooth(omega, 6, 'mean');

figure(2);clf;
subplot(221); plot(t, angTx, 'k.', t, angTy, 'b.', t, angBx, 'r.', t, angby, 'm.');
xlabel('time (s)'); ylabel('angle (deg)');
subplot(222); plot(t, wid, 'k.', t, widS, 'b-');
xlabel('time (s)'); ylabel('width (pixels)');
subplot(223); plot(t, omega, 'k.', t, omegaS, 'b-');
ylabel('rate of change of angle');
subplot(224); plot(t, [0; diff(widS)], 'b-');
grid on;
ylabel('rate of change of width');
%}