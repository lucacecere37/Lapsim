function [Fm] = calcFm(pt, tire, instaAccV, instaBackEMF)
    
    %current angular velocity of motor
    w = (instaBackEMF * pt.motor.Kv * convert('rpm','rad/s'));

    %limit imposed due to heat dissipation
    currentLimit = pt.inverter.peakPhaseCurrent;
    clT = currentLimit*pt.motor.Kt;

    %limit imposed by peak motor power
    
    %mechanical power model commented out
%     plT = pt.motor.peakPower*1000/ ...
%         (instaBackEMF * pt.motor.Kv * convert('rpm','rad/s'));
    
    %ePower = mPower + loss
    % mPower + loss - ePower = 0;
    % T*w + (T/Kt)^2 * R + instaBackEMF*Io - powerLimit = 0;
    %rewrite in form a*T^2 + b*T + c = 0;
    a = pt.motor.Rm/(pt.motor.Kt)^2;
    b = w;
    c = instaBackEMF*pt.motor.Io - pt.motor.peakPower*1000;
    
    plT = max(quadroots([ a b c ]));
    
    

    %limit imposed by instantaneous acc voltage
    maxCurrent = (instaAccV - instaBackEMF)/pt.motor.Rm;
    mcT = maxCurrent*pt.motor.Kt;

    torqueMotor = min(clT,min(mcT,plT));

    torqueTire = torqueMotor*pt.gr;
    
    Fm = max(0,torqueTire/tire.radius);

end