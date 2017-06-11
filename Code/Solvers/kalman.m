% NZ 
function [x,P] = kalman(x,y,u,A,B,C,Q,R,P)

%% Predict
xHat = A * x + B * u;
Ph = A * P * (A') + Q;

%% Update 
error = y - C * x; 
S = R + C * P * (C');
K = P * (C') * inv(S);
x = xHat + K*error; 
P = (eye(length(x)) - K * C) * Ph; 

end 