function [loss] = calcLoss(Fx,pt,tire,backEMF)
    
    %Fx*tire.radius/(pt.gr*pt.motor.Kt)
    
    tireTorque = Fx*tire.radius;
    
    motorTorque = tireTorque/pt.gr;
    
    I = motorTorque/pt.motor.Kt;
    
    loss = I^2* pt.motor.Rm + backEMF*pt.motor.Io;
    
end