function [rF, i] = simCurve(varargin)

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
    p.addOptional('calibrateCurve',false,@islogical);
    p.addOptional('vGuess', 31, @isnumeric)
    p.parse(varargin{:});
    plots = p.Results.plots;
    debug = p.Results.debug;
    car = p.Results.car;
    i = p.Results.startIndex;
    startIndex = p.Results.startIndex;
    rF = p.Results.r;
    sim = p.Results.sim;
    radius = p.Results.radius;
    theta = p.Results.theta * convert('deg','rad');
    hand = p.Results.hand;
    runUntilx = radius*theta;
    calibrateCurve = p.Results.calibrateCurve;
    vGuess = p.Results.vGuess;
    
    
    %% Assign Variable Columns
    
    assignVarColumns;
    
    %% Run Simulation
    
    if(i <= 2)
        %initial guess for long accel based on mu
        rF(1:2,longAccel) = car.tire.longMu ...
            * sim.g;
        rF(1:2,accV) = car.acc.maxVoltage;
    end
    
    if(~calibrateCurve)
        %calculate top speed for each axle, take minimum
        vMaxGuessF = car.ptF.motor.Kv*car.acc.maxVoltage/car.ptF.gr ...
            * convert('rpm','rad/s') * car.tire.radius;
        vMaxGuessR = car.ptR.motor.Kv*car.acc.maxVoltage/car.ptR.gr ...
            * convert('rpm','rad/s') * car.tire.radius;
        
        %this calculates top curve speed based on mv^2/r and downforce
        try
            vMaxGuessAero = realsqrt((car.ch.effMass*sim.g*car.tire.latMu*radius) ...
                 /(car.ch.effMass - .5 * car.aero.cl * car.aero.airDensity * ...
                 car.aero.fa * car.tire.latMu*radius));
        catch
            vMaxGuessAero = inf;
        end
        
        %this calculates top curve speed based on drag force
        %the resulting v is close to the cubed root of the power limit
        vMaxDrag = nthroot((2*car.acc.powerLimiter)/ ...
            (car.aero.airDensity*car.aero.fa*car.aero.cd),3);
        
        
        vMaxGuess = min(vMaxDrag,min(vMaxGuessAero,min(vMaxGuessF,vMaxGuessR)));
        
        %disp('Determining maximum corner speed for this feature...')
        
        %tic
        
        [res, index] = simCurve('sim', sim, 'car', car, 'radius', radius,...
            'calibrateCurve',true, 'vGuess', vMaxGuess);
        
        %toc
        
        res = res(1:index-1,:);
        
        vMax = res(end,v);
        
        %keyboard
        
        %set initial conditions
        if(vMax < rF(i-1, v))
            rF(i-1, v) = vMax;
            rF(i-1, longAccel) = 0;
            rF(i-1, yawRate) = vMax/radius * hand;
        end
    
        %disp('Corner speed successfully calculated. Simulating...')
        
        whileArg = '(rF(i-1,x) - rF(startIndex-1,x)) < runUntilx';
    else
        
        whileArg = ['abs(rF(i-1,a)) > ' num2str(sim.aErrorThresh)];
        rF(1:2,v) = vGuess;
        rF(1:2,a) = 1;
    end
    

    while(eval(whileArg))
        
        %solve for aero forces out here, only depend on v. use previous v
        rF(i,aeroDownF) = (1/2)*car.aero.airDensity*rF(i-1,v)^2 ...
            *car.aero.cl*car.aero.fa*(1-car.aero.cp);
        rF(i,aeroDownR) = (1/2)*car.aero.airDensity*rF(i-1,v)^2 ...
            *car.aero.cl*car.aero.fa*car.aero.cp;
        rF(i,aeroDrag) = (1/2)*car.aero.airDensity*rF(i-1,v)^2 ...
            *car.aero.cd*car.aero.fa;
        
        %calculate back-emfs, also only v-dependent
        rF(i,backEMFF) = rF(i-1,v)/car.tire.radius * car.ptF.gr ...
            * convert('rad/s','rpm') /car.ptF.motor.Kv;
        rF(i,backEMFR) = rF(i-1,v)/car.tire.radius * car.ptR.gr ...
            * convert('rad/s','rpm') /car.ptR.motor.Kv;
        
        %determine instantaneous lateral acceleration
        rF(i,latAccel) = rF(i-1,v)^2 / radius;
        
        %reset iteration and error
        iteration = 1;
        error = 1;
        error2 = 1;
        powerLimiter = false;
        extraLoop = 0;
        
        %solve for longitudinal acceleration at this time step
        while((error > sim.wtErrorThresh && error2 > sim.vErrorThresh) ...
            || extraLoop <= 1)
        
            if(iteration == 1)
                %calculate longitudinal weight transfer
                rF(i,longWT) = (car.ch.cgZ/car.ch.wb)*car.ch.effMass*rF(i-1,longAccel);
            else
                rF(i,longWT) = (car.ch.cgZ/car.ch.wb)*car.ch.effMass*rF(i,longAccel);
            end
            
            %calculate simplified lateral weight transfer
            rF(i,latWTF) = (car.ch.cgZ/car.ch.twF)*car.ch.effMass...
                *(1 - car.ch.staticWD)*rF(i,latAccel);
            rF(i,latWTR) = (car.ch.cgZ/car.ch.twR)*car.ch.effMass...
                *car.ch.staticWD*rF(i,latAccel);
            
            %calculate normal loads on tires
            rF(i,FzFL) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rF(i,aeroDownF) - rF(i,longWT))/2 - rF(i,latWTF)*hand;
            rF(i,FzFR) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rF(i,aeroDownF) - rF(i,longWT))/2 + rF(i,latWTF)*hand;
            rF(i,FzRL) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rF(i,aeroDownR) + rF(i,longWT))/2 - rF(i,latWTR)*hand;
            rF(i,FzRR) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rF(i,aeroDownR) + rF(i,longWT))/2 + rF(i,latWTR)*hand;
            
            %saturate at 0 (no wheelies allowed)
            rF(i,FzFL) = max(0,rF(i,FzFL));
            rF(i,FzFR) = max(0,rF(i,FzFR));
            rF(i,FzRL) = max(0,rF(i,FzRL));
            rF(i,FzRR) = max(0,rF(i,FzRR));
            
            %calculate effective long and lat mu at tires (load sensitivity)
            rF(i,longMuFL) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzFL) ...
                + car.tire.muLoadSenSq * rF(i,FzFL)^2;
            rF(i,longMuFR) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzFR) ...
                + car.tire.muLoadSenSq * rF(i,FzFR)^2;
            rF(i,longMuRL) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzRL) ...
                + car.tire.muLoadSenSq * rF(i,FzRL)^2;
            rF(i,longMuRR) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzRR) ...
                + car.tire.muLoadSenSq * rF(i,FzRR)^2;
            
            rF(i,latMuFL) = car.tire.latMu + car.tire.muLoadSenLin * rF(i,FzFL) ...
                + car.tire.muLoadSenSq * rF(i,FzFL)^2;
            rF(i,latMuFR) = car.tire.latMu + car.tire.muLoadSenLin * rF(i,FzFR) ...
                + car.tire.muLoadSenSq * rF(i,FzFR)^2;
            rF(i,latMuRL) = car.tire.latMu + car.tire.muLoadSenLin * rF(i,FzRL) ...
                + car.tire.muLoadSenSq * rF(i,FzRL)^2;
            rF(i,latMuRR) = car.tire.latMu + car.tire.muLoadSenLin * rF(i,FzRR) ...
                + car.tire.muLoadSenSq * rF(i,FzRR)^2;
            
            %calculate maximum tire forces
            maxFxFL = rF(i,FzFL) * rF(i,longMuFL);
            maxFxFR = rF(i,FzFR) * rF(i,longMuFR);
            maxFxRL = rF(i,FzRL) * rF(i,longMuRL);
            maxFxRR = rF(i,FzRR) * rF(i,longMuRR);
            
            maxFyFL = rF(i,FzFL) * rF(i,latMuFL);
            maxFyFR = rF(i,FzFR) * rF(i,latMuFR);
            maxFyRL = rF(i,FzRL) * rF(i,latMuRL);
            maxFyRR = rF(i,FzRR) * rF(i,latMuRR);
            
            maxTotalFy = (maxFyFL + maxFyFR + maxFyRL + maxFyRR);
            
            pctFyFL = maxFyFL/maxTotalFy;
            pctFyFR = maxFyFR/maxTotalFy;
            pctFyRL = maxFyRL/maxTotalFy;
            pctFyRR = maxFyRR/maxTotalFy;
            
            totalFy = rF(i,latAccel) * car.ch.effMass;
            
            rF(i,FyFL) = pctFyFL*totalFy;
            rF(i,FyFR) = pctFyFR*totalFy;
            rF(i,FyRL) = pctFyRL*totalFy;
            rF(i,FyRR) = pctFyRR*totalFy;
            

            %determine available longitudinal traction based on
            %instantaneous traction ellipse 
            %(how much friction remains for long accel)

            sqrtArgFL = (1 - (rF(i,FyFL)/maxFyFL)^2)*(maxFxFL^2);
            sqrtArgFR = (1 - (rF(i,FyFR)/maxFyFR)^2)*(maxFxFR^2);
            sqrtArgRL = (1 - (rF(i,FyRL)/maxFyRL)^2)*(maxFxRL^2);
            sqrtArgRR = (1 - (rF(i,FyRR)/maxFyRR)^2)*(maxFxRR^2);

            if(sqrtArgFL > 0)
                rF(i,FtFL) = realsqrt(sqrtArgFL);
            else
                rF(i,FtFL) = 0;
            end
            if(sqrtArgFR > 0)
                rF(i,FtFR) = realsqrt(sqrtArgFR);
            else
                rF(i,FtFR) = 0;
            end
            if(sqrtArgRL > 0)
                rF(i,FtRL) = realsqrt(sqrtArgRL);
            else
                rF(i,FtRL) = 0;
            end
            if(sqrtArgRR > 0)
                rF(i,FtRR) = realsqrt(sqrtArgRR);
            else
                rF(i,FtRR) = 0;
            end
            
            if(~powerLimiter)
                if(iteration == 1)
                    %calculate available motor force at tires
                    rF(i,FmFL) = calcFm(car.ptF,car.tire,rF(i-1,accV),rF(i,backEMFF));
                    rF(i,FmFR) = calcFm(car.ptF,car.tire,rF(i-1,accV),rF(i,backEMFF));
                    rF(i,FmRL) = calcFm(car.ptR,car.tire,rF(i-1,accV),rF(i,backEMFR));
                    rF(i,FmRR) = calcFm(car.ptR,car.tire,rF(i-1,accV),rF(i,backEMFR));
                else
                    rF(i,FmFL) = calcFm(car.ptF,car.tire,rF(i,accV),rF(i,backEMFF));
                    rF(i,FmFR) = calcFm(car.ptF,car.tire,rF(i,accV),rF(i,backEMFF));
                    rF(i,FmRL) = calcFm(car.ptR,car.tire,rF(i,accV),rF(i,backEMFR));
                    rF(i,FmRR) = calcFm(car.ptR,car.tire,rF(i,accV),rF(i,backEMFR));
                end
            end
                
            %calculate actual Fx's
            rF(i,FxFL) = min(rF(i,FtFL),rF(i,FmFL));
            rF(i,FxFR) = min(rF(i,FtFR),rF(i,FmFR));
            rF(i,FxRL) = min(rF(i,FtRL),rF(i,FmRL));
            rF(i,FxRR) = min(rF(i,FtRR),rF(i,FmRR));
            
            
            %calculate new long accel
            newLongAccel = (rF(i,FxFL) + rF(i,FxFR) + rF(i,FxRL) + rF(i,FxRR) - rF(i,aeroDrag)) ...
                /car.ch.effMass;
            
            if(iteration == 1)
                %calculate error
                error = abs(newLongAccel - rF(i-1,longAccel));

                %reset rF(i,longAccel) for next iteration- * by factor to converge 
                rF(i,longAccel) = newLongAccel - (newLongAccel - rF(i-1,longAccel))* 0.9;
            else
                %calculate error
                error = abs(newLongAccel - rF(i,longAccel));

                %reset rF(i,longAccel) for next iteration- * by factor to converge 
                rF(i,longAccel) = newLongAccel - (newLongAccel - rF(i,longAccel))*0.9;
            end
            
            
            %calculate motor losses based on tire force/speed
            rF(i,lossFL) = calcLoss(rF(i,FxFL),car.ptF,car.tire,rF(i,backEMFF));
            rF(i,lossFR) = calcLoss(rF(i,FxFR),car.ptF,car.tire,rF(i,backEMFF));
            rF(i,lossRL) = calcLoss(rF(i,FxRL),car.ptR,car.tire,rF(i,backEMFR));
            rF(i,lossRR) = calcLoss(rF(i,FxRR),car.ptR,car.tire,rF(i,backEMFR));

            %calculate mechanical power at each tire
            rF(i,mPowerFL) = rF(i-1,v)*rF(i,FxFL);
            rF(i,mPowerFR) = rF(i-1,v)*rF(i,FxFR);
            rF(i,mPowerRL) = rF(i-1,v)*rF(i,FxRL);
            rF(i,mPowerRR) = rF(i-1,v)*rF(i,FxRR);
            rF(i,mPowerTotal) = rF(i,mPowerFL) + rF(i,mPowerFR) + ...
                rF(i,mPowerRL) + rF(i,mPowerRR);


            %calculate electrical power for each tire
            %ePower = mPower + Losses

            rF(i,ePowerFL) = rF(i,mPowerFL) + rF(i,lossFL);
            rF(i,ePowerFR) = rF(i,mPowerFR) + rF(i,lossFR);        
            rF(i,ePowerRL) = rF(i,mPowerRL) + rF(i,lossRL);        
            rF(i,ePowerRR) = rF(i,mPowerRR) + rF(i,lossRR);
            rF(i,ePowerTotal) = rF(i,ePowerFL) + rF(i,ePowerFR) + rF(i,ePowerRL) ...
                + rF(i,ePowerRR);
            
            %set the regen flag
            regen = false;
            
            %engage power limiter if necessary
            if(rF(i,ePowerTotal) > car.acc.powerLimiter*1000)
                powerLimiter = true;
                
                minFxF = min(rF(i,FxFL),rF(i,FxFR));
                minFxR = min(rF(i,FxRL),rF(i,FxRR));
                
                %recalculate motor forces
                [FmF,FmR] = calcPowerLimitFm(car,rF(i-1,v),minFxF,minFxR,...
                    car.acc.powerLimiter*1000, regen);
                rF(i,FmFL) = FmF;
                rF(i,FmFR) = FmF;
                rF(i,FmRL) = FmR;
                rF(i,FmRR) = FmR;
            end
            
            [newAccV,rF(i,accI)] = accVI(rF(i,ePowerTotal),...
                car.acc.maxVoltage, car.acc.packR, regen);

            if(iteration == 1)
                %calculate error
                error2 = abs(newAccV - rF(i-1,accV));

                %reset rF(i,accV) for next iteration- * by factor to converge 
                rF(i,accV) = newAccV - (newAccV - rF(i-1,accV))* 0.9;
            else
                %calculate error
                error2 = abs(newAccV - rF(i,accV));

                %reset rF(i,accV) for next iteration- * by factor to converge 
                rF(i,accV) = newAccV - (newAccV - rF(i,accV))* 0.9;
            end


            iteration = iteration + 1;
            
            
            %check if loop exit conditions have been satisfied
            if(error <= sim.wtErrorThresh && error2 <= sim.vErrorThresh)
                extraLoop = extraLoop + 1;
            end

        end

        
        rF(i,a) = rF(i,longAccel);
        rF(i,v) = rF(i-1,v) + rF(i-1,a)*sim.dt;
        rF(i,x) = rF(i-1,x) + rF(i,v)*sim.dt;
        
        rF(i,t) = rF(i-1,t) + sim.dt;
        
        if(~calibrateCurve)
            %update yaw rate
            rF(i,yawRate) = rF(i,v)/radius * hand;

            %update heading of vehicle
            rF(i, heading) = rF(i, yawRate)*sim.dt + rF(i-1,heading);

            %update x and y of vehicle in world coordinate system
            rF(i,xCS) = (rF(i,x) - rF(i-1,x))*cos(rF(i,heading)) + rF(i-1,xCS);
            rF(i,yCS) = (rF(i,x) - rF(i-1,x))*sin(rF(i,heading)) + rF(i-1,yCS);
        end
        
%         disp([' accel = ' num2str(rF(i,a))])
%         disp([' vel = ' num2str(rF(i,v))])
%         keyboard
        i = i + 1;

        %keyboard
    end
    
    if(~calibrateCurve)
        %clean up last row due to overshoot
        fixEnd = zeros(1,length(rF(1,:)));
        for j = 1:length(fixEnd)
            fixEnd(1,j) = interp1(rF(i-2:i-1,x),rF(i-2:i-1,j),runUntilx+rF(startIndex-1,x));
        end
        rF(i-1,:) = fixEnd;
    end
    

end