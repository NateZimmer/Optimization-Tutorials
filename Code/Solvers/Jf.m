function J=Jf(fnc,params)

eps = 1e-8;
x1 = fnc(params);
m = size(x1,1);
n = length(params);
J = zeros(m,n);

for i = 1:n
    paramCpy = params; 
    paramCpy(i)= paramCpy(i) + eps;
    J(:,i) = (fnc(paramCpy) - x1)/eps;
end

end