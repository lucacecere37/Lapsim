function [ out ] = quadroots(p)
%implements quadratic equation
%   Detailed explanation goes here
    a = p(1);
    b = p(2);
    c = p(3);
    r1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    r2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
    
    out = [r1;r2];

end

