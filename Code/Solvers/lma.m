function theta = lma(Rfnc,params,iterations)

alpha = 1;
theta = params;
oldCost = norm(Rfnc(theta));

for i =1:iterations;
    r = Rfnc(theta);
    J = Jf(Rfnc,theta);
    p = -pinv(J'*J + alpha*eye(length(params)))*J'*r;
    rNew = Rfnc(theta+p);
    newCost = norm(rNew);
    if(newCost<oldCost && norm(imag(rNew))<1e-9)
        theta = theta+p;  
        oldCost = newCost;
        alpha =0.1*alpha;
    else
        alpha = 10*alpha;
        if(alpha>1e100)
            break; % Avoids infinity breaking SVD
        end
    end
end

end