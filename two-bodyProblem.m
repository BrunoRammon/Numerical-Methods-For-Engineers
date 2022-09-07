e=0.7; m1=1; m2=4;
T=2*pi./(1-e).^1.5; tspan=linspace(0,T,1000);
options=odeset('RelTol',1.e-6);
%%%%% Solve differential equations for x and y using ode45 with arguments tspan and options.
%%%%% Determine x1, y1 and x2, y2
%
%
%
w0 = 0; x0 = -1; y0 = 0; z0 = sqrt(1+e);
[t, wxyz] = ode45(@(t,wxyz)orbit(wxyz),tspan,[w0; x0; y0; z0],options);

x1=m2/(m1+m2)*wxyz(:,2);
y1=m2/(m1+m2)*wxyz(:,3);
x2=-m1/(m1+m2)*wxyz(:,2);
y2=-m1/(m1+m2)*wxyz(:,3);

function d_wxyz_dt = orbit(wxyz)
    % define the differential equation here
    w=wxyz(1); x=wxyz(2); y=wxyz(3);   z=wxyz(4);
    d_wxyz_dt = [-x./((x.^2+y.^2).^1.5); w; z; -y./((x.^2+y.^2).^1.5)];
end
%%%%% graphics: UNCOMMENT TO RUN ON MATLAB ONLINE OR DESKTOP %%%%%%%%%%%%%%
k=0.1;
R1=k*(m1)^(1/3); R2=k*(m2)^(1/3); %radius of masses
theta = linspace(0,2*pi); 
figure; axis equal; hold on; set(gcf,'color','w');
axis off; 
xlim([-2,5]); ylim([-2.5,2.5]);
planet=fill(R1*cos(theta)+x1(1), R1*sin(theta)+y1(1),'b'); 
sun=fill(R2*cos(theta)+x2(1), R2*sin(theta)+y2(1),'r'); 
pause(1);
nperiods=5; %number of periods to plot
for j=1:nperiods
    for i=1:length(t)
        planet.XData=R1*cos(theta)+x1(i); planet.YData=R1*sin(theta)+y1(i); 
        sun.XData=R2*cos(theta)+x2(i); sun.YData=R2*sin(theta)+y2(i); 
        drawnow;
    end
end
%%%%% Write local function for differential equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
