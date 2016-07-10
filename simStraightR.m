function [rR, i] = simStraightR(varargin)
    
    p = inputParser;
    p.addOptional('sim',@isstruct)
    p.addOptional('car', WR217e, @isstruct);
    p.addOptional('plots',false, @islogical);
    p.addOptional('debug',false, @islogical);
    p.addOptional('x',75,@isnumeric);
    p.addOptional('startIndex',2,@isnumeric);
    p.addOptional('r',zeros(10000,62),@ismatrix);
    p.parse(varargin{:});
    plots = p.Results.plots;
    debug = p.Results.debug;
    car = p.Results.car;
    runUntilx = p.Results.x;
    i = p.Results.startIndex;
    startIndex = p.Results.startIndex;
    rR = p.Results.r;
    sim = p.Results.sim;
    
    %TODO: fix case where one motor drags the other above accV
    %force motors to equal drag or normal torque, whichever is less
    
    %close all
    
    %% Assign Variable Columns
    
    assignVarColumns;
    %% Begin Simulation

    if( i <= 2)
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
    vMax = min(vMaxGuessAero,min(vMaxGuessF,vMaxGuessR));
    
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
        
        
        %reset iteration and error
        iteration = 1;
        error = 1;
        powerLimiter = false;
        extraLoop = 0;
        
        %solve for longitudinal acceleration at this time step
        while(error > sim.wtErrorThresh  ...
            || extraLoop <= 1)
        
            if(iteration == 1)
                %calculate longitudinal weight transfer
                rR(i,longWT) = -(car.ch.cgZ/car.ch.wb)*car.ch.effMass*rR(i-1,longAccel);
            else
                rR(i,longWT) = -(car.ch.cgZ/car.ch.wb)*car.ch.effMass*rR(i,longAccel);
            end
            
            %calculate normal loads on tires
            rR(i,FzFL) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rR(i,aeroDownF) - rR(i,longWT))/2;
            rR(i,FzFR) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rR(i,aeroDownF) - rR(i,longWT))/2;
            rR(i,FzRL) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rR(i,aeroDownR) + rR(i,longWT))/2;
            rR(i,FzRR) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rR(i,aeroDownR) + rR(i,longWT))/2;
            
            %saturate at 0 (no wheelies allowed)
            rR(i,FzFL) = max(0,rR(i,FzFL));
            rR(i,FzFR) = max(0,rR(i,FzFR));
            rR(i,FzRL) = max(0,rR(i,FzRL));
            rR(i,FzRR) = max(0,rR(i,FzRR));
            
            %calculate effective mu at tires (load sensitivity)
            rR(i,longMuFL) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzFL) ...
                + car.tire.muLoadSenSq * rR(i,FzFL)^2;
            rR(i,longMuFR) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzFR) ...
                + car.tire.muLoadSenSq * rR(i,FzFR)^2;
            rR(i,longMuRL) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzRL) ...
                + car.tire.muLoadSenSq * rR(i,FzRL)^2;
            rR(i,longMuRR) = car.tire.longMu + car.tire.muLoadSenLin * rR(i,FzRR) ...
                + car.tire.muLoadSenSq * rR(i,FzRR)^2;
        
            %calculate available traction at tires
            rR(i,FtFL) = rR(i,FzFL) * rR(i,longMuFL);
            rR(i,FtFR) = rR(i,FzFR) * rR(i,longMuFR);
            rR(i,FtRL) = rR(i,FzRL) * rR(i,longMuRL);
            rR(i,FtRR) = rR(i,FzRR) * rR(i,longMuRR);
            
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
            newLongAccel = (rR(i,FxFL) + rR(i,FxFR) + rR(i,FxRL) + rR(i,FxRR) - rR(i,aeroDrag)) ...
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
            %ePower = mPower - Losses

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
                
                %recalculate motor forces
                [FmF,FmR] = calcPowerLimitFm(car,rR(i-1,v),rR(i,FmFL),rR(i,FmRL)...
                    ,car.acc.maxRegen*1000,regen);
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
        rR(i,yawRate) = 0;

        %update heading of vehicle
        rR(i, heading) = rR(i-1,heading);
        
        %update x and y of vehicle in world coordinate system
        rR(i,xCS) = (rR(i,x) - rR(i-1,x))*cos(rR(i,heading)) + rR(i-1,xCS);
        rR(i,yCS) = (rR(i,x) - rR(i-1,x))*sin(rR(i,heading)) + rR(i-1,yCS);

        i = i + 1;
        %keyboard
    end
    
    %clean up last row due to overshoot
    fixEnd = zeros(1,length(rR(1,:)));
    
    for q = 1:length(fixEnd)
        fixEnd(1,q) = interp1(rR(i-2:i-1,x),rR(i-2:i-1,q),runUntilx+rR(startIndex-1,x));
    end
    rR(i-1,:) = fixEnd;
    
    if(plots)
        rR = rR(1:i-1,:);
        
        
        %accel in g's
        aG = rR(:,a)./sim.g;

        %velocity in mph
        vMPH = rR(:,v) * convert('m/s','mph');

        %net mechanical power at car level, W
        netMechPowerCar = rR(:,a).*rR(:,v).*car.ch.effMass;

        %acc losses
        accLosses = rR(:,accI).^2 .*car.acc.packR;

        %total power leaving pack

        accPower = rR(:,ePowerTotal) - accLosses;

        eff = rR(:,mPowerTotal)./ accPower;

        aeroDownTotal = rR(:,aeroDownF) + rR(:,aeroDownR);
    
            %create track plot colored by velocity
        figure
        scatter(rR(:, xCS), rR(:,yCS), 2, rR(:,v));
        colormap jet
        axis equal
        title('Velocity')

        figure
        scatter(rR(:, xCS), rR(:,yCS), 2, rR(:,longAccel));
        colormap jet
        axis equal
        title('Long Gs')

        figure
        scatter(rR(:, xCS), rR(:,yCS), 2, rR(:,latAccel));
        colormap jet
        axis equal
        title('Lat Gs')

        figure
        scatter(rR(:, xCS), rR(:,yCS), 2, aeroDownTotal);
        colormap jet
        axis equal
        title('Downforce')

        figure
        plot(rR(:,x),rR(:,v))
        title('Velocity vs. Distance')

        figure
        plot(rR(:,t),rR(:,a)/sim.g)
        title('Long Accel')

        figure
        plot(rR(:,t),rR(:,latAccel)/sim.g)
        title('Lat Accel')

        figure
        plot(rR(:,latAccel)/sim.g, rR(:,longAccel)/sim.g,'o-')
        title('Long Accel & Lat Accel')

        figure
        hold all
        plot(rR(:,t),(rR(:,FxFL)+rR(:,FxFR)) ... 
            ./(rR(:,FxFL)+rR(:,FxFR)+rR(:,FxRL)+rR(:,FxRR)))
        plot(rR(:,t),(rR(:,FzFL)+rR(:,FzFR)) ... 
            ./(rR(:,FzFL)+rR(:,FzFR)+rR(:,FzRL)+rR(:,FzRR)))
        plot(rR(:,t),(rR(:,FtFL)+rR(:,FtFR)) ... 
            ./(rR(:,FtFL)+rR(:,FtFR)+rR(:,FtRL)+rR(:,FtRR)))
        legend('Powertrain Force % Front','Normal Force % Front', 'Vehicle Grip % Front')
        hold off

        figure
        hold all
        plot(rR(:,t),rR(:,FxFL))
        plot(rR(:,t),rR(:,FxRL))
        hold off

        figure
        plot(rR(:,t),rR(:,ePowerTotal))
        hold all
        plot(rR(:,t),netMechPowerCar)
        plot(rR(:,t),rR(:,mPowerTotal))
        plot(rR(:,t),accPower)
        hold off
        legend('Electrical Power','Car-level Mechanical Power', ...
            'Sum of Motor Mechanical Powers','Total Acc Power (incl. losses)')

        figure
        plot(rR(:,t),rR(:,accV))
        hold all
        plot(rR(:,t),rR(:,backEMFF)); 
        plot(rR(:,t),rR(:,backEMFR)); 
        hold off

        figure
        plot(rR(:,t),eff)

        figure
        plot(rR(:,t),rR(:,ePowerFL)); 
        hold all
        plot(rR(:,t),rR(:,ePowerRL)); 
        hold off

        %sanity check
        FzCheck = rR(:,FzFL) + rR(:,FzFR) + rR(:,FzRL) + rR(:,FzRR) - aeroDownTotal;

        figure
        plot(rR(:,t),FzCheck)

        figure
        plot(rR(:,t),rR(:,FzFL))
        hold all
        plot(rR(:,t),rR(:,FzRL))
        plot(rR(:,t),rR(:,FzFR))
        plot(rR(:,t),rR(:,FzRR))
        legend('FzFL','FzRL','FzFR','FzRR')
%        figure
%         plot(rR(:,t), aG)
%         ylim([ 0 1.1*max(aG)])
%         xlabel('Time (sec)')
%         ylabel('Longitudinal Acceleration (g)')
%     %     
%        figure
%         plot(rR(:,t), vMPH)
%         ylim([ 0 1.1*max(vMPH)])
%         xlabel('Time (sec)')
%         ylabel('Velocity (mph)')
    %     
    %     
    %     figure
    %     hold all
    %     plot(t,FtFL)
    %     plot(t,FmFL)
    %     hold off
    %     legend('Available Traction','Motor Force')
    %     title('Front Motor')
    %     
    %     figure
    %     hold all
    %     plot(t,FtRL)
    %     plot(t,FmRL)
    %     hold off
    %     legend('Available Traction','Motor Force')
    %     title('Rear Motor')
    %     

        %pctFzF = 
        figure
        hold all
        plot(rR(:,t),(rR(:,FmFL)+rR(:,FmFR)) ... 
            ./(rR(:,FmFL)+rR(:,FmFR)+rR(:,FmRL)+rR(:,FmRR)))
        plot(rR(:,t),(rR(:,FzFL)+rR(:,FzFR)) ... 
            ./(rR(:,FzFL)+rR(:,FzFR)+rR(:,FzRL)+rR(:,FzRR)))
        plot(rR(:,t),(rR(:,FtFL)+rR(:,FtFR)) ... 
            ./(rR(:,FtFL)+rR(:,FtFR)+rR(:,FtRL)+rR(:,FtRR)))
        legend('Powertrain Force % Front','Normal Force % Front', 'Vehicle Grip % Front')
        hold off

        figure
        hold all
        plot(rR(:,t),rR(:,FxFL))
        plot(rR(:,t),rR(:,FxRL))
        hold off

        figure
        plot(rR(:,t),rR(:,ePowerTotal))
        hold all
        %plot(rR(:,t),netMechPowerCar)
        plot(rR(:,t),rR(:,mPowerTotal))
        %plot(rR(:,t),accPower)
        hold off
        %legend('Electrical Power','Car-level Mechanical Power', ...
        %    'Sum of Motor Mechanical Powers','Total Acc Power (incl. losses)')

        figure
        plot(rR(:,t),rR(:,accV))
        hold all
        plot(rR(:,t),rR(:,backEMFF)); 
        plot(rR(:,t),rR(:,backEMFR)); 
        hold off

%         figure
%         plot(rR(:,t),eff)

        figure
        plot(rR(:,t),rR(:,ePowerFL)); 
        hold all
        plot(rR(:,t),rR(:,ePowerRL)); 
        hold off

    %     %what torque would have been
    %     torqueWouldR = (car.acc.maxVoltage - backEMFR)/car.ptR.motor.Rm * car.ptR.motor.Kt * car.ptR.gr;
    %     
    %     %what torque should have been
    %     torqueShouldR = (accV - backEMFR)/car.ptR.motor.Rm * car.ptR.motor.Kt * car.ptR.gr;
    %     
    %     %what torque was
    %     torqueWasR = FxRR * car.tire.radius;
    %     
    %     %what torque would have been
    %     torqueWouldF = (car.acc.maxVoltage - backEMFF)/car.ptF.motor.Rm * car.ptF.motor.Kt * car.ptF.gr;
    %     
    %     %what torque should have been
    %     torqueShouldF = (accV - backEMFF)/car.ptF.motor.Rm * car.ptF.motor.Kt * car.ptF.gr;
    %     
    %     %what torque was
    %     torqueWasF = FxFR * car.tire.radius;
    %     
    %     figure
    %     plot(t,torqueWouldR)
    %     hold all
    %     plot(t,torqueShouldR)
    %     plot(t,torqueWasR)
    %     hold off
    %     xlim([ 2.4 time])
    %     ylim([-500 2000])
    %     legend('would','should','was')
    %     
    %     figure
    %     plot(t,torqueWouldF)
    %     hold all
    %     plot(t,torqueShouldF)
    %     plot(t,torqueWasF)
    %     hold off
    %     xlim([ 2.4 time])
    %     ylim([-500 2000])
    %     legend('would','should','was')
    
    end
    
    if(debug)
        keyboard
    end
end