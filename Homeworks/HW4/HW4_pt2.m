%% part 3(a)
theta = 0;
phi = 0;
k = 1;
sigma = 2/k;
sigma_x = sigma;
sigma_y = sigma;
A = 1;


[x, y] = meshgrid(-5:.1:5);
z = Ds(x,y,phi,sigma_x,sigma_y,k);
surf(x,y,z)
xlim([-5,5])
ylim([-5,5])
xlabel('x (degrees)','fontsize',16)
ylabel('y (degrees)','fontsize',16)
title('Gabor Function Representation as 2D Surface','fontsize',16);
width=900;
height=600;
set(gcf,'position',[0,0,width,height])
%% Now Plot the Image
I = image(z,'CDataMapping','scaled');

set(gca,'XTickLabel',-4:5)
set(gca,'yTickLabel',-4:5)
xlabel('x (degrees)','fontsize',16)
ylabel('y (degrees)','fontsize',16)
title('Gabor Function Representation as Colorplot','fontsize',16);
width=900;
height=600;
set(gcf,'position',[0,0,width,height])

%% Part 3(b)
tau = 0:0.5:300;

plot(tau,Dt(tau,1/15),'LineWidth',2)
set(gca,'XDir','reverse')
xlabel('\tau (ms)','fontsize',16)
ylabel('Dt (Hz)','fontsize',16)
title('Temporal Evolution of Spatial Receptive Field','fontsize',16);
width=900;
height=600;
set(gcf,'position',[0,0,width,height])

%% Part 3(c)
theta = 0;
phi = 0;
k = 1;
sigma = 2/k;
sigma_x = sigma;
sigma_y = sigma;
A = 1;


tau = 50:50:300;
alpha = 1/15;
[x, y] = meshgrid(-5:.1:5);
for i=1:numel(tau)
     z = D(x,y,phi,sigma_x,sigma_y,k,tau(i),alpha);
     subplot(2,3,i);  
     contour(x,y,z,'LineWidth',2)
     title(['\tau =',' ',num2str(tau(i)),'ms']);
     xlabel('x (degrees)')
     ylabel('y (degrees)')
     width=900;
     height=600;
     set(gcf,'position',[0,0,width,height])
end
colorbar('eastoutside')

%% Part 3(d)
t = 0;
phi = 0;
A = 1;
w=800;

x = 0:0.1:200;
y = 0:0.1:200;
K = 0.2;
theta = 0.5*pi;

z = s(x,y,t,K,theta,w,A,phi);

R = image(z,'CDataMapping','scaled');
title('Counterphase Grating: t = 0','FontSize',16);
set(gca,'xtick',[])
set(gca,'ytick',[])

width=900;
height=600;
set(gcf,'position',[0,0,width,height])

%% Part 3(d) Continued... **Please run the last block first** Now, let's make a movie!

%t = 0:.1:2;
%theta = 0:.1:2*pi;

t= linspace(0,2,20);
theta = linspace(0,2*pi,20);
for i=1:numel(t)
    z = s(x,y,t(i),K,theta(i),w,A,phi);
    V = image(z,'CDataMapping','scaled');
    title(['Counterphase Grating: t =','',num2str(t(i))],'FontSize',16);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
    pause(.3);


end
%% 3(d) continued
t = 0:.5:2;

for i=1:6
    z = s(x,y,t(i),K,theta,w,A,phi);
    subplot (1,6,i)
    V = image(z,'CDataMapping','scaled');
    title(['t =','',num2str(t(i))],'FontSize',16);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    width=900;
    height=600;
    set(gcf,'position',[0,0,width,height])
    
end




%% Part 3(e) 

%Reproduce Figure A - Vary theta
clear variables
theta = 0;
phi = 0;
k = 1;
sigma = 2/k;
sigma_x = sigma;
sigma_y = sigma;
A = 1;
phi = 0;
K =k;

theta = -1.5:.02:1.5;

for i=1:numel(theta)
    
    fun = @(x,y) integrand(x,y,A,K,theta(i),phi,sigma_x,sigma_y,k);
    Ls_A(i) = quad2d(fun,-100,100,-100,100);
    
end

plot(theta,Ls_A,'LineWidth',2)
title('Orientation dependent Spatial Linear Response ','fontsize',16);
xlabel('\theta (radians)','fontsize',16)
ylabel('L_{s}','fontsize',16)
width=900;
height=600;
set(gcf,'position',[0,0,width,height])
%% Part 3(e) Continued
%Reproduce Figure B - Vary ratio K/k (vary K set k = 1)
clear variables
theta = 0;
phi = 0;
k = 1;
sigma = 2/k;
sigma_x = sigma;
sigma_y = sigma;
A = 1;
phi = 0;

K = 0:0.02:3;
for i=1:numel(K)
    fun = @(x,y) integrand(x,y,A,K(i),theta,phi,sigma_x,sigma_y,k);
    Ls_B(i) = quad2d(fun,-100,100,-100,100);


end

plot(K,Ls_B,'LineWidth',2)
title('Frequency dependent Spatial Linear Response ','fontsize',16);
xlabel('K/k','fontsize',16)
ylabel('L_{s}','fontsize',16)
width=900;
height=600;
set(gcf,'position',[0,0,width,height])

%% Part 3(e) Continued
%Reproduce Figure C - Vary phase phi
clear variables
theta = 0;
k = 1;
sigma = 2/k;
sigma_x = sigma;
sigma_y = sigma;
A = 1;
K = k;

phi = -2.5:0.02:2.5;
for i=1:numel(phi)
    fun = @(x,y) integrand(x,y,A,K,theta,phi(i),sigma_x,sigma_y,k);
    Ls_C(i) = quad2d(fun,-50,50,-50,50);


end

plot(phi,Ls_C,'LineWidth',2)
title('Frequency dependent Spatial Linear Response ','fontsize',16);
xlabel('\phi','fontsize',16)
ylabel('L_{s}','fontsize',16)
width=900;
height=600;
set(gcf,'position',[0,0,width,height])

%% Part 3(f)
clear variables
alpha = 1/15;
dt = 0.2;
tau = 0:dt:1000;
t = 0:dt:1000;
w = 0:.01*.02:2*pi * .02;

for i=1:numel(w)
    func1 = Dt(tau,alpha);
    func2 = cos(w(i)*(t));
    convolution = dt*conv(func2,func1,'same');
    pks=max(convolution);
    amplitude(i) = (pks);
end
plot(w/(2*pi)*1000,amplitude,'linewidth',2)  
title('Model Frequency Response','fontsize',16);
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Amplitude','fontsize',16)
width=900;
height=600;
set(gcf,'position',[0,0,width,height])


%% Functions
function out = Ds(x,y,phi,sigma_x,sigma_y,k)

[out,] = (1/(2*pi*sigma_x*sigma_y))*exp(-(x.^2)/(2*(sigma_x.^2)) - (y.^2)/(2*(sigma_y.^2))).*cos(k*x - phi);

end

function out = Dt(tau,alpha)
[out, ] = alpha*exp(-alpha*tau).*( ((alpha*tau).^5)/(125) - ((alpha*tau).^7)/5040);
end

function out = D(x,y,phi,sigma_x,sigma_y,k,tau,alpha)
[out, ] = Ds(x,y,phi,sigma_x,sigma_y,k).*Dt(tau,alpha);
end

function out = s(x,y,t,K,theta,w,A,phi)
[out, ] = A*cos( (K*x*cos(theta)) + (K*y*sin(theta)) - phi).*cos(w*t);
end

function out = integrand(x,y,A,K,theta,phi,sigma_x,sigma_y,k)
[out, ] = Ds(x,y,phi,sigma_x,sigma_y,k).*(A*cos(K.*x.*cos(theta) + K.*y.*sin(theta) - phi));
end

function out = integrand1(t,tau,w,alpha)
[out, ] = Dt(tau,alpha) .* cos(w(t-tau));
end

