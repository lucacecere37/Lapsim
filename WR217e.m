function [car] = WR217e()

%% Powertrain Parameters
car.acc.powerLimiter = 80; %kW
car.acc.maxRegen = 20; %kW
car.acc.SOC = 100;

%Energus Pack with HG2 cells
car.acc.cellR = .025;
car.acc.sCount = 30;
car.acc.pCount = 24;

% %Energus Pack with 25R cells
% car.acc.cellR = .023;
% car.acc.sCount = 30;
% car.acc.pCount = 24;
% 
% %Melasta Pack
% car.acc.cellR = .0015;
% car.acc.sCount = 28;
% car.acc.pCount = 5;

car.ptF.motor.Kv = 110; %rpm/V
car.ptR.motor.Kv = 51; %rpm/V

car.ptF.motor.peakPower = 15; %kW
car.ptR.motor.peakPower = 30;
car.ptF.inverter.peakPhaseCurrent = 280; %A
car.ptR.inverter.peakPhaseCurrent = 280; %A
car.ptF.gr = 7.5;
car.ptR.gr = 4.5;
car.ptF.motor.Rm = .0128; %ohms
car.ptR.motor.Rm = .0086; %ohms
car.ptF.motor.Io = 17; %A
car.ptR.motor.Io = 8.5; %A
car.grMultiplier = 1;

%% Chassis Parameters
car.ch.effMass = 280; %kg
car.ch.cgZ = 0.25; %m
car.ch.wb = 1.5748; %m
car.ch.staticWD = 0.5; %less is more frontwards
car.ch.twF = 1.1938; %track width front
car.ch.twR = 1.1430; %track width rear

%% Tire Parameters
car.tire.radius = 0.22352; %m
car.tire.longMu = 1.4;
car.tire.latMu = 1.7;
car.tire.muLoadSenLin = 2.2851e-5; % mu/N 12psi 6in rim
car.tire.muLoadSenSq = -1.20311e-7; %mu/N^2

%% Aero Parameters
car.aero.cl = 2; %coeff of lift/downforce
car.aero.cp = 0.52; %center of pressure, less is more frontwards
car.aero.cd = 1.25; %coeff of drag
car.aero.airDensity = 1.2; %kg/m^3
car.aero.fa = 1.4; %frontal aera, m^2

%% Modifiers
     


end
        