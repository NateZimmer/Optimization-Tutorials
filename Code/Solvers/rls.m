% NZ RLS
function [theta,P] = rls(x,y,theta,lambda,P)

P = P/lambda - (P * (x') * x * P)./(lambda + x * P * x')./lambda;
theta = theta + P*(x')*(y - x * theta);

end