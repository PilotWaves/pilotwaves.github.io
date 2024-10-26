%Ethanol diffusion in tumor
tic

dr = 0.005;
dt = dr;
nu = dt/dr;
mu = dt/dr^2;
T = 2;
X = 1;
N = round(T/dt + 1);
M = round(X/dr + 1);

PVol = 0.1;
Vol = PVol*4*pi/3;
U_thresh = 0.25;
eps = 5;


r = [0:dr:1];

U = zeros(M,N);
A = zeros(M);
B = A;
width = 1/10;
r1 = round(length(r)*width);
height = Vol/(sum((r(1:r1 - 1).^2).*exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2)))*dr*4*pi);
U(1:r1 - 1,1) = exp(1 - 1./(1 - (r(1:r1 - 1)/r(r1)).^2))*(height);

 writerObj = VideoWriter('EthanolPartial.mp4');
 writerObj.FrameRate = round(N/100);
 open(writerObj);

%Plotting Initial Condition

hfig = figure;
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1]);
   
subplot(1,2,1);
   x = r'*cos([0:2*pi/(length(r)-1):2*pi]);
y = r'*sin([0:2*pi/(length(r)-1):2*pi]);
z = U(:,1)*ones(size(r));
col = 1./z;  % This is the color, vary with u in this case.
colormap(jet)
surf(x,y,z,col,...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
axis([-1.5 1.5 -1.5 1.5],'square');
caxis([0 max(U(1,:))/100])
title('Radial Diffusion');
axis off;

subplot(1,2,2);
h = plot(r,U(:,1),'b',-r,fliplr(U(:,1)),'b',[-1,1],[U_thresh U_thresh],'r--');
set(h(1),'Linewidth',3);
set(h(2),'Linewidth',3);
set(h(3),'Linewidth',2);
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
xlabel('Radius');
ylabel('Concentration');
title('Concentration Profile');
%axis([-1 1 0 height]);

   pause
   
    frame = getframe(hfig);
 writeVideo(writerObj,frame);
 pcount = 0;
   
for n = 2:N
    
    %Point at origin and boundary
     A(1,1) = mu + 1; A(1,2) = -mu;
     A(M,M) = mu + 1 + eps*dr*(nu + mu); A(M,M-1) = -mu;
    
   %Points in the middle
   for m = 2:M-1
      A(m,m+1) = -nu/(2*r(m)) - mu/2; A(m,m) = mu + 1;
        A(m,m-1) = nu/(2*r(m)) - mu/2;
   end
   
   
   U(:,n) = (A^(-1))*U(:,n-1);
   
   %Plotting
   
   subplot(1,2,1);
z = U(:,n)*ones(size(r));
col = 1./z;  % This is the color, vary with u in this case.
colormap(jet)
surf(x,y,z,col,...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
axis([-1.5 1.5 -1.5 1.5],'square');
caxis([0 max(U(1,:))/100])
title('Radial Diffusion');
axis off;

subplot(1,2,2);
 h = plot(r,U(:,n),'b',-r,fliplr(U(:,n)),'b',[-1 1],[U_thresh U_thresh],'r--');
set(h(1),'Linewidth',3);
set(h(2),'Linewidth',3);
set(h(3),'Linewidth',2);
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
xlabel('Radius');
ylabel('Concentration');
title('Concentration Profile');
if max(U(:,n)) >= height/10
axis([-1 1 0 height]);
n_crit1 = n;
else if max(U(:,n)) >= U(:,n_crit1)/10
        axis([-1 1 0 max(U(:,n_crit1))]);
        n_crit2 = n;
    else if max(U(:,n)) >= U(:,n_crit2)/10
        axis([-1 1 0 max(U(:,n_crit2))]);
        n_crit3 = n;
        pause
        else
            axis([-1 1 0 max(U(:,n_crit3))]);
        end
    end
end

if U(M,n) > 0.01 && pcount == 0
   pause
   pcount = pcount + 1;
end
   
    frame = getframe(hfig);
 writeVideo(writerObj,frame);
 
      %%% Makes it blink %%%
 if U(M,n) > U_thresh
     pause
         subplot(1,2,1);
z = U(:,n)*ones(size(r));
col = 1./z;  % This is the color, vary with u in this case.
colormap(jet)
surf(x,y,z,col,...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
axis([-1.5 1.5 -1.5 1.5],'square');
caxis([0 max(U(1,:))/100])
title('Radial Diffusion');
axis off;

subplot(1,2,2);
h = plot(r,U(:,n),'w',-r,fliplr(U(:,n)),'w',[-1 1],[U_thresh U_thresh],'w--');
set(h(1),'Linewidth',3);
set(h(2),'Linewidth',3);
set(h(3),'Linewidth',2);
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
xlabel('Radius');
ylabel('Concentration');
title('Concentration Profile');

if max(U(:,n)) >= height/10
axis([-1 1 0 height]);
n_crit1 = n;
else if max(U(:,n)) >= U(:,n_crit1)/10
        axis([-1 1 0 max(U(:,n_crit1))]);
        n_crit2 = n;
    else if max(U(:,n)) >= U(:,n_crit2)/10
        axis([-1 1 0 max(U(:,n_crit2))]);
        n_crit3 = n;

        else
            axis([-1 1 0 max(U(:,n_crit3))]);
        end
    end
end
    frame = getframe(hfig);
 writeVideo(writerObj,frame);
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 
 if abs(U(1,n) - U(M,n)) <= 0.0001 || U(M,n) > U_thresh*1.1 || U(1,n) < U_thresh*0.9
     break;
 end
 
end

 close(writerObj);

toc