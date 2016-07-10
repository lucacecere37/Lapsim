function [car] = initializeCar(car)

car.ptF.gr = car.ptF.gr * car.grMultiplier;
car.ptR.gr = car.ptR.gr * car.grMultiplier;
car.acc.maxVoltage = 4.2*car.acc.sCount; %V
car.acc.packR = car.acc.cellR/car.acc.pCount * car.acc.sCount;
car.ptF.motor.Kt =  1/(car.ptF.motor.Kv * convert('rpm','rad/s'));
car.ptR.motor.Kt =  1/(car.ptR.motor.Kv * convert('rpm','rad/s'));