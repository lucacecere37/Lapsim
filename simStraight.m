function [rF, i] = simStraight(varargin)
    
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
    rF = p.Results.r;
    sim = p.Results.sim;

    
    %set the regen flag
    regen = false;
    
    
    %TODO: fix case where one motor drags the other above accV
    %force motors to equal drag or normal torque, whichever is less
    
    %close all
    
    %% Assign Variable Columns
    
    assignVarColumns;
    
    %% Begin Simulation

    if( i <= 2)
        %initial guess for long accel based on mu
        rF(1:2,longAccel) = car.tire.longMu ...
            * sim.g;
        rF(1:2,accV) = car.acc.maxVoltage;
    end
    
    while((rF(i-1,x) - rF(startIndex-1,x)) < runUntilx)
        
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
        
        
%         aeroDownF(i) = 0;
%         aeroDownR(i) = 0;
%         rF(i,aeroDrag) = 0;
        
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
            
            %calculate normal loads on tires
            rF(i,FzFL) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rF(i,aeroDownF) - rF(i,longWT))/2;
            rF(i,FzFR) = (car.ch.effMass * (1 - car.ch.staticWD) * sim.g ...
                + rF(i,aeroDownF) - rF(i,longWT))/2;
            rF(i,FzRL) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rF(i,aeroDownR) + rF(i,longWT))/2;
            rF(i,FzRR) = (car.ch.effMass * car.ch.staticWD * sim.g ...
                + rF(i,aeroDownR) + rF(i,longWT))/2;
            
            %saturate at 0 (no wheelies allowed)
            rF(i,FzFL) = max(0,rF(i,FzFL));
            rF(i,FzFR) = max(0,rF(i,FzFR));
            rF(i,FzRL) = max(0,rF(i,FzRL));
            rF(i,FzRR) = max(0,rF(i,FzRR));
            
            %calculate effective mu at tires (load sensitivity)
            rF(i,longMuFL) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzFL) ...
                + car.tire.muLoadSenSq * rF(i,FzFL)^2;
            rF(i,longMuFR) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzFR) ...
                + car.tire.muLoadSenSq * rF(i,FzFR)^2;
            rF(i,longMuRL) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzRL) ...
                + car.tire.muLoadSenSq * rF(i,FzRL)^2;
            rF(i,longMuRR) = car.tire.longMu + car.tire.muLoadSenLin * rF(i,FzRR) ...
                + car.tire.muLoadSenSq * rF(i,FzRR)^2;
            
            %calculate available traction at tires
            rF(i,FtFL) = rF(i,FzFL) * rF(i,longMuFL);
            rF(i,FtFR) = rF(i,FzFR) * rF(i,longMuFR);
            rF(i,FtRL) = rF(i,FzRL) * rF(i,longMuRL);
            rF(i,FtRR) = rF(i,FzRR) * rF(i,longMuRR);
            
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
            
        
            
            %engage power limiter if necessary
            if(rF(i,ePowerTotal) > car.acc.powerLimiter*1000)
                powerLimiter = true;
                
                %recalculate motor forces
                [FmF,FmR] = calcPowerLimitFm(car,rF(i-1,v),rF(i,FxFL),rF(i,FxRL),car.acc.powerLimiter*1000,regen);
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

        %update yaw rate
        rF(i,yawRate) = 0;

        %update heading of vehicle
        rF(i, heading) = rF(i-1,heading);
        
        %update x and y of vehicle in world coordinate system
        rF(i,xCS) = (rF(i,x) - rF(i-1,x))*cos(rF(i,heading)) + rF(i-1,xCS);
        rF(i,yCS) = (rF(i,x) - rF(i-1,x))*sin(rF(i,heading)) + rF(i-1,yCS);

        i = i + 1;

        %keyboard
    end
    
    %clean up last row due to overshoot
    fixEnd = zeros(1,length(rF(1,:)));
    for j = 1:length(fixEnd)
        fixEnd(1,j) = interp1(rF(i-2:i-1,x),rF(i-2:i-1,j),runUntilx+rF(startIndex-1,x));
    end
    rF(i-1,:) = fixEnd;
    
    if(plots)
       figure
        plot(rF(:,t), aG)
        ylim([ 0 1.1*max(aG)])
        xlabel('Time (sec)')
        ylabel('Longitudinal Acceleration (g)')
    %     
       figure
        plot(rF(:,t), vMPH)
        ylim([ 0 1.1*max(vMPH)])
        xlabel('Time (sec)')
        ylabel('Velocity (mph)')
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
        plot(rF(:,t),(rF(:,FxFL)+rF(:,FxFR)) ... 
            ./(rF(:,FxFL)+rF(:,FxFR)+rF(:,FxRL)+rF(:,FxRR)))
        plot(rF(:,t),(rF(:,FzFL)+rF(:,FzFR)) ... 
            ./(rF(:,FzFL)+rF(:,FzFR)+rF(:,FzRL)+rF(:,FzRR)))
        plot(rF(:,t),(rF(:,FtFL)+rF(:,FtFR)) ... 
            ./(rF(:,FtFL)+rF(:,FtFR)+rF(:,FtRL)+rF(:,FtRR)))
        legend('Powertrain Force % Front','Normal Force % Front', 'Vehicle Grip % Front')
        hold off

        figure
        hold all
        plot(rF(:,t),rF(:,FxFL))
        plot(rF(:,t),rF(:,FxRL))
        hold off

        figure
        plot(rF(:,t),rF(:,ePowerTotal))
        hold all
        plot(rF(:,t),netMechPowerCar)
        plot(rF(:,t),rF(:,mPowerTotal))
        plot(rF(:,t),accPower)
        hold off
        legend('Electrical Power','Car-level Mechanical Power', ...
            'Sum of Motor Mechanical Powers','Total Acc Power (incl. losses)')

        figure
        plot(rF(:,t),rF(:,accV))
        hold all
        plot(rF(:,t),rF(:,backEMFF)); 
        plot(rF(:,t),rF(:,backEMFR)); 
        hold off

        figure
        plot(rF(:,t),eff)

        figure
        plot(rF(:,t),rF(:,ePowerFL)); 
        hold all
        plot(rF(:,t),rF(:,ePowerRL)); 
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