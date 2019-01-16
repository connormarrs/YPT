% Hammer_Euler
%
% 2019 YPT Hammer problem
% Numerical model of Euler's Equations
% Connor Marrs, RCDS June 2019

clear;

%init_w_array = 1:6; % initial ang speed around z axis 
init_w = 12; % initial ang speed around z axis 
eps = 0.01;   % initial ang speed around x axis (to become r)
M = 1;  % mass of block

% orientation of block
a = 10; % size in x-axis: lying along launcher
b = 5;  % size in y axis: across launcher
c = 1;  % size in vertical z axis

I1 = M/12*(b^2 + c^2); %moment of inertia around x axis
I2 = M/12*(a^2 + c^2); % around y
I3 = M/12*(a^2 + b^2); % around z
Istr = ['I_x = ' num2str(I1,'%6.2e') ': eps      ';...
        'I_y = ' num2str(I2,'%6.2e') ': fast spin';...
        'I_z = ' num2str(I3,'%6.2e') ': eps      ']; %string for display/key
%time array parameters    
NT = 10001; %number of time steps
dt = .001; % length of each time step
t = (0:NT-1)*dt;

%initialize rotational kinematic variables
[w1, w2, w3, a1, a2, a3] = deal(nan(1,NT));
w1(1) = eps;
w2(1) = init_w; %rad/s
w3(1) = eps;

for n=1:NT-1
    %Euler's equations
  a1(n) = -(I3-I2)*w2(n)*w3(n)/I1; % alpha
  a2(n) = -(I1-I3)*w3(n)*w1(n)/I2;
  a3(n) = -(I2-I1)*w1(n)*w2(n)/I3;
  
  %kinematics to find omega at next time
  w1(n+1) = w1(n) + a1(n)*dt; % omega
  w2(n+1) = w2(n) + a2(n)*dt;
  w3(n+1) = w3(n) + a3(n)*dt;
end
%vector sum total of omega
wtot = sqrt(w1.^2 + w2.^2 + w3.^2);

%angular momentum
L1 = I1*w1;
L2 = I2*w2;
L3 = I3*w3;
Ltot = sqrt(L1.^2 + L2.^2 + L3.^2);

q1 = dt*cumsum(w1)/2/pi;
q2 = dt*cumsum(w2)/2/pi;
q3 = dt*cumsum(w3)/2/pi;


figure(1); clf;
ax(1) = subplot(2,2,1); %set(ax(1),'position',[.1  .8 .38 .15]);
ax(2) = subplot(2,2,2); %set(ax(2),'position',[.55 .8 .38 .15]);
ax(3) = subplot(2,2,3); %set(ax(3),'position',[.1  .1 .38 .55]);
ax(4) = subplot(2,2,4); %set(ax(4),'position',[.55 .1 .38 .55]);
subplot(ax(1));
%[H1] = plot_rect3D(a,b,c);

% subplot(ax(2)); axis off
% text(0,.5, Istr);
% title('units of mass and length cancel');

subplot(ax(3));

H3(1) = plot(t,w1,'k'); hold on;
H3(2) = plot(t,w2,'b'); 
H3(3) = plot(t,w3,'r'); 
H3(4) = plot(t,wtot,'g');
legend(H3,{'w_x','w_y','w_z','w_t_o_t'});
xlabel('time (s)'); ylabel('Angular Speed (rad/s)');
grid on;

subplot(ax(4));
H4(1) = plot(t,L1,'k'); hold on;
H4(2) = plot(t,L2,'b');
H4(3) = plot(t,L3,'r');
H4(4) = plot(t,Ltot,'g');
legend(H4,{'L_x','L_y','L_z','L_t_o_t'});
xlabel('time (s)'); ylabel('Angular Momentum');
grid on;


subplot(ax(2));
H2(1) = plot(t,q1,'k'); hold on;
H2(2) = plot(t,q2,'b');
H2(3) = plot(t,q3,'r');
%H2(4) = plot(t,sqrt(q1.^2 + q2.^2 + q3.^2),'g');
legend(H2,{'q_x','q_y','q_z'});%,'q_t_o_t'});
xlabel('time (s)'); ylabel('Angular Position (rot)');
grid on;

return


tmp = max(Ltot)/10;
figure(2); clf;
for k = 1:100:NT
  ang = [q1(k) q2(k) q3(k)];
  cla;
  [Hbox] = plot_rect3D(a,b,c,ang);
  title(['time step ' num2str(k) ' of ' num2str(NT)])
  HL0 = plot3([0 L1(k)]/tmp,[0 L2(k)]/tmp,[0 L3(k)]/tmp,...
    'b','linewidth',2);l
  axis([-1 1 -1 1 -1 1]*max([a b c])*1.5)
  pause
end

return

figure(2); clf;
for k=1:50:2500
  if k>1
    set(Hw0,'linestyle',':','linewidth',1);
    set(HL0,'linestyle',':','linewidth',1);
  end
  Hw0 = plot3([0 w1(k)],[0 w2(k)],[0 w3(k)],'k','linewidth',2); hold on
  HL0 = plot3([0 L1(k)],[0 L2(k)],[0 L3(k)],'b','linewidth',2);
  plot3([0 w1(k)],[0 0],[0 0],'k')
  plot3([0 0],[0 w2(k)],[0 0],'k')
  plot3([0 0],[0 0],[0 w3(k)],'k')
  plot3([0 L1(k)],[0 0],[0 0],'b')
  plot3([0 0],[0 L2(k)],[0 0],'b')
  plot3([0 0],[0 0],[0 L3(k)],'b')
%   set(gca,'xlim',[0 max(L1)]);
%   set(gca,'ylim',[0 max(L2)]);
%   set(gca,'zlim',[0 max(L3)]);
  pause;
end

  
