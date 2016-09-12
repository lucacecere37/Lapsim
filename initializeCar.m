function [car] = initializeCar(car)

if(~isfield(car.acc,'SOC'))
    car.acc.SOC = 100;
end

car.ptF.gr = car.ptF.gr * car.grMultiplier;
car.ptR.gr = car.ptR.gr * car.grMultiplier;
car.acc.maxVoltage = (3 + 1.2*(car.acc.SOC/100))*car.acc.sCount; %V
car.acc.packR = car.acc.cellR/car.acc.pCount * car.acc.sCount;
car.ptF.motor.Kt =  1/(car.ptF.motor.Kv * convert('rpm','rad/s'));
car.ptR.motor.Kt =  1/(car.ptR.motor.Kv * convert('rpm','rad/s'));