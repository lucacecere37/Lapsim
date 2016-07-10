function [FmF, FmR] = calcPowerLimitFm(car,v,currentFxF,currentFxR, maxP, regen)

%equation is of the following form:
%total electrical power = total mechanical power + losses
%ePowerTotal = mechPowerTotal + motor copper losses + motor iron losses
%ePowerTotal = 2*[(FmF + FmR)*v + ((FmF*r)/(grF * KtF))^2 * RmF + ((FmR*r)/(grR * KtR))^2 * RmR ...
% + BEMFF*IoF + BEMFR*IoR]

%rewrite in the following form:

% a*(FmF)^2 + b*FmF + c*(FmR)^2 + d*FmR = e


%maximum electrical power through energy meter in W
%maxP = car.acc.powerLimiter * 1000;

%backEMF 
backEMFF = v/car.tire.radius * car.ptF.gr * convert('rad/s','rpm') ...
    /car.ptF.motor.Kv;
backEMFR = v/car.tire.radius * car.ptR.gr * convert('rad/s','rpm') ...
    /car.ptR.motor.Kv;

if(~regen)
    %calculation for e
    e = maxP - 2*backEMFF*car.ptF.motor.Io - 2*backEMFR*car.ptR.motor.Io;

    %calculation for a
    a = 2*(car.tire.radius/(car.ptF.gr*car.ptF.motor.Kt))^2 * car.ptF.motor.Rm;

    %calculation for c
    c = 2*(car.tire.radius/(car.ptR.gr*car.ptR.motor.Kt))^2 * car.ptR.motor.Rm;

    %calculation for b
    b = 2*v;

    %calculation for d
    d = 2*v;
else
    e = maxP + 2*backEMFF*car.ptF.motor.Io + 2*backEMFR*car.ptR.motor.Io;
    a = -2*(car.tire.radius/(car.ptF.gr*car.ptF.motor.Kt))^2 * car.ptF.motor.Rm;
    c = -2*(car.tire.radius/(car.ptR.gr*car.ptR.motor.Kt))^2 * car.ptR.motor.Rm;
    b = 2*v;
    d = 2*v;
end
    
    
%keep currentFxF (x), solve for FxR (y)
%now x is known, rearrange to solve for y
%ax^2 + bx + cy^2 + dy = e
%cy^2 + dy + [ax^2 + bx - e] = 0;
%replace term in brackets with constant, f

%keeping front, solving for rear
xF = currentFxF;
f = a*xF^2 + b*xF - e;

if(~regen)
    yF = max(roots([c d f]));
else
    yF = min(roots([c d f]));
end

lossesF = 2*calcLoss(xF,car.ptF,car.tire,backEMFF);
lossesR = 2*calcLoss(yF,car.ptR,car.tire,backEMFR);

lossKeepF = lossesF+lossesR;

%keeping rear, solving for front
yR = currentFxR;
f = c*yR^2 + d*yR - e;

if(~regen)
    xR = max(roots([a b f]));
else
    xR = min(roots([a b f]));
end

lossesF = 2*calcLoss(xR,car.ptF,car.tire,backEMFF);
lossesR = 2*calcLoss(yR,car.ptR,car.tire,backEMFR);

lossKeepR = lossesF+lossesR;

%calculate force distribution at optimal level
% FmF = c/a * FmR, plug this in to solve for each

% first term a*(FmF)^2 is now c^2/a * FmR^2
% second term b*FmF is now (b*c)/a * FmR
%add these to original FmR^2 and FmR coeffs

%c^2/a * FmR^2 + c * FmR^2, new coeff g
g = abs(c^2/a) + abs(c);

%(b*c)/a * FmR + d * FmR, new coeff h
h = abs((b*c)/a) + abs(d);

%solve for optimal rear force, yO
yO = max(roots([g h -e]));

%solve for optimal front force, xO
xO = (c/a) * yO;

if(xO <= currentFxF && yO <= currentFxR)
    FmF = xO;
    FmR = yO;
else
    if(lossKeepF >= lossKeepR)
        FmF = xR;
        FmR = yR;
    else
        FmF = xF;
        FmR = yF;
    end
end



end