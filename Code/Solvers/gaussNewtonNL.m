function theta = gaussNewton(Rfnc,theta,itLimit)
%Jf(fnc,params)
for i = 1:itLimit;
    J = Jf(Rfnc,theta);
    r = Rfnc(theta);
    g = (J')*r;
    Hinv = pinv((J')*J);
    p = -Hinv*g;
    theta = theta + p; % Update Theta
end

end