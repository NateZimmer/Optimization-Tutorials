%NZ

x = (-1.4:0.01:2.5)'; % This is the time samples 
y = 5*sin(log(x)); % This is the true model
noise = 5*randn(length(x),1); % Noise model 
yMeasured = y + noise; % This is the measurement. 
J = [x.^5,x.^4,x.^3,x.^2,x.^1,x.^0]; % Formulate Jacobian 
theta = pinv(J' * J) * J' * yMeasured % Normal Equation
yCalculated = J*theta; % Our predicted model

%% Below Code is to make a pretty plot 

fig = figure(1); % Make a figure
clf(fig) % Erase figure
hold on % Hold plot
error = (yMeasured - y).^2; % Error
scatter(x,yMeasured,[],error); % Plot measured data 
plot(x, yCalculated ,'r--','LineWidth',2) % Plot Results
grid on
title('Poly fit')
set(gca,'FontSize',10,'FontWeight','bold');
set(gcf,'Units','Pixels');
set(gcf, 'Position', [500, 500, 700, 350]);
