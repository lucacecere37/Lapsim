function [rR, i] = simCurveR(varargin)

    p = inputParser;
    p.addOptional('sim',@isstruct)
    p.addOptional('car', WR217e, @isstruct);
    p.addOptional('plots',false, @islogical);
    p.addOptional('debug',false, @islogical);
    p.addOptional('radius', 30, @isnumeric);
    p.addOptional('theta', 45, @isnumeric);
    p.addOptional('startIndex',2,@isnumeric);
    p.addOptional('r',zeros(100000,62),@ismatrix);
    p.addOptional('hand',1,@isnumeric);
    p.parse(varargin{:});
    plots = p.Results.plots;
    debug = p.Results.debug;
    car = p.Results.car;
    i = p.Results.startIndex;
    startIndex = p.Results.startIndex;
    rR = p.Results.r;
    sim = p.Results.sim;
    radius = p.Results.radius;
    theta = p.Results.theta * convert('deg','rad');
    hand = p.Results.hand;
    runUntilx = radius*theta;
    
    
    %% Assign Variable Columns
    
    assignVarColumns;
    
    %% Run Simulation
    
    if(i <= 2)
        %initial guess for long accel based on mu
        rR(1:2,longAccel) = car.tire.longMu ...
            * sim.g;
        rR(1:2,accV) = car.acc.maxVoltage;
    end
    
    %calculate top speed for each axle, take minimum
    vMaxGuessF = car.ptF.motor.Kv*car.acc.maxVoltage/car.ptF.gr ...
        * convert('rpm','rad/s') * car.tire.radius;
    vMaxGuessR = car.ptR.motor.Kv*car.acc.maxVoltage/car.ptR.gr ...
        * convert('rpm','rad/s') * car.tire.radius;
    try
        vMaxGuessAero = realsqrt((car.ch.effMass*sim.g*car.tire.latMu*radius) ...
             /(car.ch.effMass - .5 * car.aero.cl * car.aero.airDensity * ...
             car.aero.fa * car.tire.latMu*radius));
    catch
        vMaxGuessAero = inf;
    end
    vMaxGuess = min(vMaxGuessAero,min(vMaxGuessF,vMaxGuessR));
    
    %disp('Determining maximum corner speed for this feature...')

    [res, index] = simCurve('sim', sim, 'car', car, 'radius', radius,...
        'calibrateCurve',true, 'vGuess', vMaxGuess);
    
    

    res = res(1:index-1,:);

    vMax = res(end,v);
    
    %disp('Corner speed successfully calculated. Simulating...')

    %keyboard

    %set initial conditions
    if(vMax < rR(i-1, v))
        rR(i-1, v) = vMax;
        rR(i-1, longAccel) = 0;
        rR(i-1, yawRate) = vMax/radius * hand;
    end
    

    while((rR(i-1,x) - rR(startIndex-1,x)) < runUntilx)
        
        %solve for aero forces out here, only depend on v. use previous v
        rR(i,aeroDownF) = (1/2)*car.aero.airDensity*rR(i-1,v)^2 ...
            *car.aero.cl*car.aero.fa*(1-car.aero.cp);
        rR(i,aeroDownR) = (1/2)*car.aero.airDensity*rR(i-1,v)^2 ...
            *car.aero.cl*car.aero.fa*car.aero.cp;
        rR(i,aeroDrag) = (1/2)*car.aero.airDensity*rR(i-1,v)^2 ...
            *car.aero.cd*car.aero.fa;
        
        %calculate back-emfs, also only v-dependent
        rR(i,backEMFF) = rR(i-1,v)/car.tire.radius * car.ptF.gr ...
            * convert('rad/s','rpm') /car.ptF.motor.Kv;
        rR(i,backEMFR) = rR(i-1,v)/car.tire.radius * car.ptR.gr ...
            * convert('rad/s','rpm') /car.ptR.motor.Kv;
        
        %determine instantaneous lateral acceleration
        rR(i,latAccel) = rR(i-1,v)^2 / radius;
        
        %reset iteration and error
        iteration = 1;
        error = 1;
        powerLimiter = false;
        extraLoop = 0;
        
        %solve for longitudinal acceleration at this time step
        while(error > sim.wtErrorThresh || extraLoop <= 1)
        
            if(iteration == 1)
                %calculate longitudinal weight transfer
                rR(i,longWT) = -(car.ch.cgZ/car.ch.wb)*car.ch.effMass*rR(i-1,longAccel);
            else
                rR(i,longWT) = -(car.ch.cgZ/car.ch.wb)*car.ch.effMass*rR(i,longAccel);
            end
            
            %calculate simplified lateral weight transfer
            rR(i,latWTF) = (car.ch.cgZ/car.ch.twF)*car.ch.effMass...
                *(1 - car.ch.staticWD)*rR(i,latAccel);
            rR(i,latWTR) = (car.ch.cgZ/car.ch.twR)*car.ch.effMass...
                *car.ch.staticWD*rR(i,latAccel);
            
            %calculate normal loads on tires
            rR(i,FzFL) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rR(i,aeroDownF) - rR(i,longWT))/2 - rR(i,latWTF)*hand;
            rR(i,FzFR) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rR(i,aeroDownF) - rR(i,longWT))/2 + rR(i,latWTF)*hand;
            rR(i,FzRL) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rR(i,aeroDownR) + rR(i,longWT))/2 - rR(i,latWTR)*hand;
            rR(i,FzRR) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rR(i,aeroDownR) + rR(i,longWT))/2 + rR(i,latWTR)*hand;
            
            %saturate at 0 (no wheelies allowed)
            rR(i,FzFL) = max(0,rR(i,FzFL));
            rR(i,FzFR) = max(0,rR(i,FzFR));
            rR(i,FzRL) = max(0,rR(i,FzRL));
            rR(i,FzRR) = max(0,rR(i,FzRR));
            
            %calculate effective long and lat mu at tires (load sensitivity)
            rR(i,longMuFL) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzFL) ...
                + car.tire.muLoadSenSq * rR(i,FzFL)^2;
            rR(i,longMuFR) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzFR) ...
                + car.tire.muLoadSenSq * rR(i,FzFR)^2;
            rR(i,longMuRL) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzRL) ...
                + car.tire.muLoadSenSq * rR(i,FzRL)^2;
            rR(i,longMuRR) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzRR) ...
                + car.tire.muLoadSenSq * rR(i,FzRR)^2;
            
            rR(i,latMuFL) = car.tire.latMu + car.tire.muLoadSenLin * rR(i,FzFL) ...
                + car.tire.muLoadSenSq * rR(i,FzFL)^2;
            rR(i,latMuFR) = car.tire.latMu + car.tire.muLoadSenLin * rR(i,FzFR) ...
                + car.tire.muLoadSenSq * rR(i,FzFR)^2;
            rR(i,latMuRL) = car.tire.latMu + car.tire.muLoadSenLin * rR(i,FzRL) ...
                + car.tire.muLoadSenSq * rR(i,FzRL)^2;
            rR(i,latMuRR) = car.tire.latMu + car.tire.muLoadSenLin * rR(i,FzRR) ...
                + car.tire.muLoadSenSq * rR(i,FzRR)^2;
            
            %calculate maximum tire forces
            maxFxFL = rR(i,FzFL) * rR(i,longMuFL);
            maxFxFR = rR(i,FzFR) * rR(i,longMuFR);
            maxFxRL = rR(i,FzRL) * rR(i,longMuRL);
            maxFxRR = rR(i,FzRR) * rR(i,longMuRR);
            
            maxFyFL = rR(i,FzFL) * rR(i,latMuFL);
            maxFyFR = rR(i,FzFR) * rR(i,latMuFR);
            maxFyRL = rR(i,FzRL) * rR(i,latMuRL);
            maxFyRR = rR(i,FzRR) * rR(i,latMuRR);
            
            maxTotalFy = (maxFyFL + maxFyFR + maxFyRL + maxFyRR);
            
            pctFyFL = maxFyFL/maxTotalFy;
            pctFyFR = maxFyFR/maxTotalFy;
            pctFyRL = maxFyRL/maxTotalFy;
            pctFyRR = maxFyRR/maxTotalFy;
            
            totalFy = rR(i,latAccel) * car.ch.effMass;
            
            rR(i,FyFL) = pctFyFL*totalFy;
            rR(i,FyFR) = pctFyFR*totalFy;
            rR(i,FyRL) = pctFyRL*totalFy;
            rR(i,FyRR) = pctFyRR*totalFy;
            

            %determine available longitudinal traction based on
            %instantaneous traction ellipse 
            %(how much friction remains for long accel)

            sqrtArgFL = (1 - (rR(i,FyFL)/maxFyFL)^2)*(maxFxFL^2);
            sqrtArgFR = (1 - (rR(i,FyFR)/maxFyFR)^2)*(maxFxFR^2);
            sqrtArgRL = (1 - (rR(i,FyRL)/maxFyRL)^2)*(maxFxRL^2);
            sqrtArgRR = (1 - (rR(i,FyRR)/maxFyRR)^2)*(maxFxRR^2);

            if(sqrtArgFL > 0)
                rR(i,FtFL) = realsqrt(sqrtArgFL);
            else
                rR(i,FtFL) = 0;
            end
            if(sqrtArgFR > 0)
                rR(i,FtFR) = realsqrt(sqrtArgFR);
            else
                rR(i,FtFR) = 0;
            end
            if(sqrtArgRL > 0)
                rR(i,FtRL) = realsqrt(sqrtArgRL);
            else
                rR(i,FtRL) = 0;
            end
            if(sqrtArgRR > 0)
                rR(i,FtRR) = realsqrt(sqrtArgRR);
            else
                rR(i,FtRR) = 0;
            end
            
            if(~powerLimiter)
                if(iteration == 1)
                    %calculate available motor force at tires
                    rR(i,FmFL) = calcFm(car.ptF,car.tire,rR(i-1,accV),rR(i,backEMFF));
                    rR(i,FmFR) = calcFm(car.ptF,car.tire,rR(i-1,accV),rR(i,backEMFF));
                    rR(i,FmRL) = calcFm(car.ptR,car.tire,rR(i-1,accV),rR(i,backEMFR));
                    rR(i,FmRR) = calcFm(car.ptR,car.tire,rR(i-1,accV),rR(i,backEMFR));
                else
                    rR(i,FmFL) = calcFm(car.ptF,car.tire,rR(i,accV),rR(i,backEMFF));
                    rR(i,FmFR) = calcFm(car.ptF,car.tire,rR(i,accV),rR(i,backEMFF));
                    rR(i,FmRL) = calcFm(car.ptR,car.tire,rR(i,accV),rR(i,backEMFR));
                    rR(i,FmRR) = calcFm(car.ptR,car.tire,rR(i,accV),rR(i,backEMFR));
                end
            end
                
            %calculate actual Fx's
            rR(i,FxFL) = rR(i,FtFL);
            rR(i,FxFR) = rR(i,FtFR);
            rR(i,FxRL) = rR(i,FtRL);
            rR(i,FxRR) = rR(i,FtRR);
            
            FxTotal = rR(i,FxFL) + rR(i,FxFR) + rR(i,FxRL) + rR(i,FxRR);
            
            if(rR(i-1,v) >= vMax)
                rR(i,aeroDrag) = FxTotal;
            end
            
            %calculate new long accel
            newLongAccel = (rR(i,FxFL) + rR(i,FxFR) + rR(i,FxRL) + rR(i,FxRR) + rR(i,aeroDrag)) ...
                /car.ch.effMass;
            
            if(iteration == 1)
                %calculate error
                error = abs(newLongAccel - rR(i-1,longAccel));

                %reset rR(i,longAccel) for next iteration- * by factor to converge 
                rR(i,longAccel) = newLongAccel - (newLongAccel - rR(i-1,longAccel))* 0.9;
            else
                %calculate error
                error = abs(newLongAccel - rR(i,longAccel));

                %reset rR(i,longAccel) for next iteration- * by factor to converge 
                rR(i,longAccel) = newLongAccel - (newLongAccel - rR(i,longAccel))*0.9;
            end
            
            
            %calculate motor losses based on tire force/speed
            rR(i,lossFL) = calcLoss(rR(i,FmFL),car.ptF,car.tire,rR(i,backEMFF));
            rR(i,lossFR) = calcLoss(rR(i,FmFR),car.ptF,car.tire,rR(i,backEMFF));
            rR(i,lossRL) = calcLoss(rR(i,FmRL),car.ptR,car.tire,rR(i,backEMFR));
            rR(i,lossRR) = calcLoss(rR(i,FmRR),car.ptR,car.tire,rR(i,backEMFR));

            %calculate mechanical power at each tire
            rR(i,mPowerFL) = rR(i-1,v)*rR(i,FmFL);
            rR(i,mPowerFR) = rR(i-1,v)*rR(i,FmFR);
            rR(i,mPowerRL) = rR(i-1,v)*rR(i,FmRL);
            rR(i,mPowerRR) = rR(i-1,v)*rR(i,FmRR);
            rR(i,mPowerTotal) = rR(i,mPowerFL) + rR(i,mPowerFR) + ...
                rR(i,mPowerRL) + rR(i,mPowerRR);


            %calculate electrical power for each tire
            %ePower = mPower + Losses

            rR(i,ePowerFL) = rR(i,mPowerFL) - rR(i,lossFL);
            rR(i,ePowerFR) = rR(i,mPowerFR) - rR(i,lossFR);        
            rR(i,ePowerRL) = rR(i,mPowerRL) - rR(i,lossRL);        
            rR(i,ePowerRR) = rR(i,mPowerRR) - rR(i,lossRR);
            rR(i,ePowerTotal) = rR(i,ePowerFL) + rR(i,ePowerFR) + rR(i,ePowerRL) ...
                + rR(i,ePowerRR);
            
            %set the regen flag
            regen = true;
            
            %engage power limiter if necessary
            if(rR(i,ePowerTotal) > car.acc.maxRegen*1000)
                powerLimiter = true;
                
                minFxF = min(rR(i,FmFL),rR(i,FmFR));
                minFxR = min(rR(i,FmRL),rR(i,FmRR));
                
                %recalculate motor forces
                [FmF,FmR] = calcPowerLimitFm(car,rR(i-1,v),minFxF,minFxR,...
                    car.acc.maxRegen*1000,regen);
                rR(i,FmFL) = FmF;
                rR(i,FmFR) = FmF;
                rR(i,FmRL) = FmR;
                rR(i,FmRR) = FmR;
            end

            [rR(i,accV),rR(i,accI)] = accVI(rR(i,ePowerTotal),...
                car.acc.maxVoltage, car.acc.packR, regen);


            iteration = iteration + 1;
            
            %check if loop exit conditions have been satisfied
            if(error <= sim.wtErrorThresh)
                extraLoop = extraLoop + 1;
            end

        end

        
        rR(i,a) = rR(i,longAccel);
        rR(i,v) = rR(i-1,v) + rR(i-1,a)*sim.dt;
        rR(i,x) = rR(i-1,x) + rR(i,v)*sim.dt;
        
        rR(i,t) = rR(i-1,t) + sim.dt;
        
        %correct velocity if required
        if(rR(i,v) > vMax)
            rR(i,v) = vMax;
        end
        
        %update yaw rate
        rR(i,yawRate) = rR(i,v)/radius * hand;

        %update heading of vehicle
        rR(i, heading) = rR(i, yawRate)*sim.dt + rR(i-1,heading);

        %update x and y of vehicle in world coordinate system
        rR(i,xCS) = (rR(i,x) - rR(i-1,x))*cos(rR(i,heading)) + rR(i-1,xCS);
        rR(i,yCS) = (rR(i,x) - rR(i-1,x))*sin(rR(i,heading)) + rR(i-1,yCS);


        i = i + 1;
    end
    
    %clean up last row due to overshoot
    fixEnd = zeros(1,length(rR(1,:)));
    for j = 1:length(fixEnd)
        fixEnd(1,j) = interp1(rR(i-2:i-1,x),rR(i-2:i-1,j),runUntilx+rR(startIndex-1,x));
    end
    rR(i-1,:) = fixEnd;

    
