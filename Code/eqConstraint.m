% NZ InEquality Constraints 
clc; clear all; 
m = 2; b=-2; e= 2;
t = (0 : 0.5 : 10)'; 
yM = m .* t + b + e .* randn(size(t));
rFncs = @(T) yM - (T(1)*t + T(2)) + - 0.001*log(T(2)-T(1)) + 1e9*(T(2)+T(1)-1).^2 ; 
thetaI = [-5;5];
thetaLMA = lma(rFncs,thetaI,500);
thetaLMA = real(thetaLMA);
J = [t,t.*0+1];
yH = J*thetaLMA; % model output
thetaLSQ = pinv(J'*J)*J'*yM;
yLSQ = J*thetaLSQ;

%% Plot fig 1 

fig1 = figure(1);
clf(fig1);
hold on
scatter(t,yM);
plot(t,yH)
plot(t,yLSQ,'--','color','black')
grid on 

titleStr = ['min r(\theta)^2 = (data - f(\theta,x))^2 , f(\theta,x) = \theta_0 * x + \theta_1','\newlines.t:  \theta: \theta_0<\theta_1  ,  \theta_0 + \theta_1 = 1\newline\theta_{calc} = [',num2str(thetaLMA(1)),',',num2str(thetaLMA(2)),']'];
title(titleStr)
legend('Data','Model','LSQ','location','se')
xlabel('x')

%% 3d Plot Data

yMM = -5:0.1:5;
xMM = -5:0.1:5;
Z = zeros(length(yMM),length(xMM));
Zc1 = 0 .* Z;
for i = 1:length(yMM)
   for j = 1:length(xMM)
       Z(i,j) = norm(yM - (xMM(j)*t + yMM(i)));
       Zc1(i,j) = 5;
   end
end

%% Plot fig 2 

fig2 = figure(2);
clf(fig2);
surf(yMM,xMM,Z)
shading interp
hold on 

M = zeros(4,3);
M(1,:) = [-5,-5,0];
M(2,:) = [5,5,0];
M(3,:) = [5,5,200];
M(4,:) = [-5,-5,200];
M = M';

patch(M(1,:),M(2,:),M(3,:),'w','FaceAlpha',0.7);

M = zeros(4,3);
M(1,:) = [5,-4,0];
M(2,:) = [-4,5,0];
M(3,:) = [-4,5,200];
M(4,:) = [5,-4,200];
M = M';

patch(M(1,:),M(2,:),M(3,:),'m','FaceAlpha',0.5,'EdgeColor','m');


xlabel('\theta_0');
ylabel('\theta_1');
z = norm(rFncs(thetaLMA))+5;
scatter3(thetaLMA(1),thetaLMA(2),z,'r','filled');
scatter3(thetaLSQ(1),thetaLSQ(2),norm(yM - (t*thetaLSQ(1) + thetaLSQ(2)))+3,'g','filled');
view(25.2,55.6);

fig3 = copyobj(gcf,0);
view(2);
title(titleStr)
legend('Cost','Inequality Bound','Equality Bound','Constrainted Optimal','Global Optimal','location','NW')