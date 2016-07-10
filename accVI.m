function [CCV, I] = accVI(ePower, OCV, packR, regen)

%usable electrical power = (Vnom - I*R)*I
% ePower = OCV*I - I^2 *R

%rearrange as a*I^2 + b*I + c = 0;

a = -packR;
b = OCV;

if(regen == true)
    c = ePower;
else
    c = -ePower;
end

I = min(roots([ a b c]));
CCV = (OCV - I*packR);

end