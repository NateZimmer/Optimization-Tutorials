function theta = LMS(Rfnc,params,iterations)

alpha = 1;
theta = params;
oldCost = norm(Rfnc(theta));

for i =1:iterations;
    r = Rfnc(theta);
    J = Jf(Rfnc,theta);
    p = -pinv(J'*J + alpha*eye(length(params)))*J'*r;
    newCost = norm(Rfnc(theta+p));
    if(newCost<oldCost)
        theta = theta+p;  
        oldCost = newCost;
        alpha =0.1*alpha;
    else
        alpha = 10*alpha;
    end
end

end

