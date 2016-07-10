function [time, noRegenEndurokWh, endurokWh] = simLap(varargin)
    
    %% Parse Inputs
    
    p = inputParser;
    p.addOptional('car', WR217e, @isstruct);
    p.addOptional('plots',false, @islogical);
    p.parse(varargin{:})
    car = p.Results.car;
    plots = p.Results.plots;
    
    i = 2;
    rF = zeros(1000000,62);
    
    %% Calculate Additional Car Parameters
    
    car = initializeCar(car);

    
    %% Simulation Parameters
    
    sim.dt = 0.001;
    sim.g = 9.81;
    sim.wtErrorThresh = 0.1; %0.05
    sim.vErrorThresh = 0.01;
    sim.aErrorThresh = 0.05;
    
    %% Variable Column Assignments
    
    assignVarColumns;
    
    %% Load track
    
    track = LincolnEnduro2015;
    rF(1:2,heading) = -pi/2;
    
    %% Perform Simulation Forwards
    
    for k = 1:length(track(:,1))
        disp(['Forwards simulating track feature ' num2str(k)...
            ' of ' num2str(length(track(:,1))) '...']);
        
        tic
        
        if(track(k,1) == 1)
            [rF, i] = simStraight('sim', sim, 'car', car, 'x', track(k,2),...
                'startIndex', i, 'r', rF);
        elseif(track(k,1) == 0)
            [rF, i] = simCurve('sim', sim, 'car', car, 'radius', track(k,2),...
                'theta', track(k,3), 'hand', track(k,4), ...
                'startIndex', i, 'r', rF);
        end
        
        toc

    end 
    
    
    %remove first row
    rF = rF(2:end,:);
    %remove extra rows
    rF = rF(1:i-2,:);
    
    %% Perform Simulation in Reverse
    
    %reset i
    i = 2;
    
    %create reverse matrix
    rR = zeros(1000000,62);
    
    %flip track feature order
    trackB = flipud(track);
    
    %reverse turn handedness
    trackB(:,4) = -trackB(:,4);
    
    %Some initial conditions based on final conditions of the forwards lap
    %rR(1,:) = rF(end,:);
    %make adjustments to reorient and reset the vehicle 
    rR(1,x) = 0;
    rR(1,longAccel) = 0;
    rR(1,heading) = rF(end, heading) + pi;
    rR(1,yawRate) = -rF(end, yawRate);
    rR(1,xCS) = rF(end,xCS);
    rR(1,yCS) = rF(end,yCS);
    rR(1,v) = rF(end,v);
    rR(1,latAccel) = rF(end,latAccel);
    rR(1,aeroDownF) = rF(end,aeroDownF);
    rR(1,aeroDownR) = rF(end,aeroDownR);
    rR(1,aeroDrag) = rF(end,aeroDrag);
    rR(1,backEMFF) = rF(end,backEMFF);
    rR(1,backEMFR )= rF(end,backEMFR);
    rR(1,accV) = rF(end,accV);
    rR(1,accI) = rF(end,accI);
    
    for k = 1:length(trackB(:,1))
        disp(['Reverse simulating track feature ' num2str(k)...
            ' of ' num2str(length(trackB(:,1))) '...']);
        
        tic
        
        if(trackB(k,1) == 1)
            [rR, i] = simStraightR('sim', sim, 'car', car, 'x', trackB(k,2),...
                'startIndex', i, 'r', rR);
        elseif(trackB(k,1) == 0)
            [rR, i] = simCurveR('sim', sim, 'car', car, 'radius', trackB(k,2),...
                'theta', trackB(k,3), 'hand', trackB(k,4), ...
                'startIndex', i, 'r', rR);
        end
        
        toc

    end 
    
    %remove first row
    rR = rR(2:end,:);
    %remove extra rows
    rR = rR(1:i-2,:);
    
    
    %% Superimpose F and R Simulations
    
    disp('Simulations complete. Processing output...')
    
    %reset rR distances to match rF distances
    for p = 1:length(rR(:,x))
        rR(p,x) = rR(end,x) - rR(p,x);
    end
    
    %initialize the transformed reverse matrix
    rRT = zeros(size(rF));

    %Transform the reverse matrix to match indices with the forwards matrix
    for k = 1:length(rR(1,:))
        rRT(:,k) = interp1(rR(:,x),rR(:,k),rF(:,x));
    end
    

    %initialize distance interp r matrix
    rD = zeros(size(rF));
    
    %we want to keep the simulation results from the run with the lower
    %instantaneous velocity
    
    rfFlagColNum = length(rF(1,:)) + 1;
    
    rD(:,rfFlagColNum) = zeros(length(rF(:,x)),1);
    
    %set a flag to -1 if motion came from reverse run
    for m = 1:length(rF(:,x))
        if(rF(m,v) < rRT(m,v))
            rD(m,1:end-1) = rF(m,:);
            rD(m,rfFlagColNum) = 1;
        else
            rD(m,1:end-1) = rRT(m,:);
            rD(m,rfFlagColNum) = -1;
        end
    end
    
    %remove last two rows because of NaNs -why does this happen?
    rD = rD(1:end-2,:);
    
    %correct positive/negative signs
    rD(:,a) = rD(:,a).*rD(:,rfFlagColNum);
    rD(:,longAccel) = rD(:,longAccel).*rD(:,rfFlagColNum);
    rD(:,a) = rD(:,a).*rD(:,rfFlagColNum);
    rD(:,FmFL) = rD(:,FmFL).*rD(:,rfFlagColNum);
    rD(:,FmFR) = rD(:,FmFR).*rD(:,rfFlagColNum);
    rD(:,FmRL) = rD(:,FmRL).*rD(:,rfFlagColNum);
    rD(:,FmRR) = rD(:,FmRR).*rD(:,rfFlagColNum);
    rD(:,FxFL) = rD(:,FxFL).*rD(:,rfFlagColNum);
    rD(:,FxFR) = rD(:,FxFR).*rD(:,rfFlagColNum);
    rD(:,FxRL) = rD(:,FxRL).*rD(:,rfFlagColNum);
    rD(:,FxRR) = rD(:,FxRR).*rD(:,rfFlagColNum);
    rD(:,mPowerFL) = rD(:,mPowerFL).*rD(:,rfFlagColNum);
    rD(:,mPowerFR) = rD(:,mPowerFR).*rD(:,rfFlagColNum);
    rD(:,mPowerRL) = rD(:,mPowerRL).*rD(:,rfFlagColNum);    
    rD(:,mPowerRR) = rD(:,mPowerRR).*rD(:,rfFlagColNum); 
    rD(:,mPowerTotal) = rD(:,mPowerTotal).*rD(:,rfFlagColNum);
    rD(:,ePowerFL) = rD(:,ePowerFL).*rD(:,rfFlagColNum);
    rD(:,ePowerFR) = rD(:,ePowerFR).*rD(:,rfFlagColNum);
    rD(:,ePowerRL) = rD(:,ePowerRL).*rD(:,rfFlagColNum);    
    rD(:,ePowerRR) = rD(:,ePowerRR).*rD(:,rfFlagColNum); 
    rD(:,ePowerTotal) = rD(:,ePowerTotal).*rD(:,rfFlagColNum);
    rD(:,accI) = rD(:,accI) .* rD(:,rfFlagColNum);
    
            
    %fix the t column because it is totally wrong due to distance interp
    rD(1,t) = 0;
    %throw out first row if zero velocity
    if rD(1,v) == 0
        rD = rD(2:end,:);
    end
    
    for g = 2:length(rD(:,t))
        
        dx = rD(g,x) - rD(g-1,x);
        instaVelocity = rD(g-1,v);
        
        rD(g,t) = dx/instaVelocity;
       
    end
    
    %now the t column is full of dts, change to running t
    rD(:,t) = cumsum(rD(:,t));
    
    time = rD(end,t);
    
    %interp again to get back to constant dt increments
    %this will be the final output r matrix
    
    tQuery = linspace(rD(1,t),time,(time-rD(1,t))/sim.dt)';
    
    r = zeros(length(tQuery),length(rD(1,:))-1);
    for k = 1:length(rD(1,:))-1
        r(:,k) = interp1(rD(:,t),rD(:,k),tQuery);
    end
    
    %% Perform Post-Processing Calculations
    
    %net mechanical power at car level, W
    netMechPowerCar = r(:,a).*r(:,v).*car.ch.effMass;
    
    %acc losses
    accLosses = r(:,accI).^2 .*car.acc.packR;
    
    %total power leaving pack
    
    accPower = r(:,ePowerTotal) + accLosses;
    
    aeroDownTotal = r(:,aeroDownF) + r(:,aeroDownR);
    
    aggkWh = cumtrapz(r(:,t),accPower).*convert('J','kWh');
    
    kWh = aggkWh(end);
    
    endurokWh = 22000/r(end,x) * kWh;
    
    noRegenAccPower = max(0,r(:,ePowerTotal)) + accLosses;
    
    noRegenAggkWh = cumtrapz(r(:,t),noRegenAccPower).*convert('J','kWh');
    
    noRegenkWh = noRegenAggkWh(end);
    
    noRegenEndurokWh = 22000/r(end,x) * noRegenkWh;
    
    motorTorqueFR = abs(r(:,mPowerFR)./r(:,v) .* car.tire.radius ./ car.ptF.gr);
    motorTorqueFL = abs(r(:,mPowerFL)./r(:,v) .* car.tire.radius ./ car.ptF.gr);
    motorTorqueRL = abs(r(:,mPowerRL)./r(:,v) .* car.tire.radius ./ car.ptR.gr);
    motorTorqueRR = abs(r(:,mPowerRR)./r(:,v) .* car.tire.radius ./ car.ptR.gr);
    
    mcCurrentFL = abs(r(:,ePowerFL)./r(:,accV));
    mcCurrentFR = abs(r(:,ePowerFR)./r(:,accV));
    mcCurrentRL = abs(r(:,ePowerRL)./r(:,accV));
    mcCurrentRR = abs(r(:,ePowerRR)./r(:,accV));
    
    if plots
        %% Default Plots 

        % %generate color map
        % rColumn = linspace(1,0)';
        % gColumn = linspace(0,1)';
        % bColumn = zeros(100,1);
        % %map = [rColumn, gColumn, bColumn];

        %create track plot colored by velocity
        figure
        scatter(r(:, xCS), r(:,yCS), 2, r(:,v));
        colormap jet
        axis equal
        title('Velocity')

        figure
        scatter(r(:, xCS), r(:,yCS), 2, r(:,longAccel));
        colormap jet
        axis equal
        title('Long Gs')

        figure
        scatter(r(:, xCS), r(:,yCS), 2, r(:,latAccel));
        colormap jet
        axis equal
        title('Lat Gs')

        figure
        scatter(r(:, xCS), r(:,yCS), 2, aeroDownTotal);
        colormap jet
        axis equal
        title('Downforce')

        figure
        plot(r(:,x),r(:,v))
        title('Velocity vs. Distance')

        figure
        plot(r(:,t),r(:,a)/sim.g)
        title('Long Accel')

        figure
        plot(r(:,t),r(:,latAccel)/sim.g)
        title('Lat Accel')

        figure
        plot(r(:,latAccel)/sim.g, r(:,longAccel)/sim.g,'o')
        title('Long Accel & Lat Accel')

        figure
        hold all
        plot(r(:,t),(r(:,FxFL)+r(:,FxFR)) ... 
            ./(r(:,FxFL)+r(:,FxFR)+r(:,FxRL)+r(:,FxRR)))
        plot(r(:,t),(r(:,FzFL)+r(:,FzFR)) ... 
            ./(r(:,FzFL)+r(:,FzFR)+r(:,FzRL)+r(:,FzRR)))
        plot(r(:,t),(r(:,FtFL)+r(:,FtFR)) ... 
            ./(r(:,FtFL)+r(:,FtFR)+r(:,FtRL)+r(:,FtRR)))
        legend('Powertrain Force % Front','Normal Force % Front', 'Vehicle Grip % Front')
        hold off

        figure
        hold all
        plot(r(:,t),r(:,FxFL))
        plot(r(:,t),r(:,FxRL))
        hold off

        figure
        plot(r(:,t),r(:,ePowerTotal))
        hold all
        plot(r(:,t),netMechPowerCar)
        plot(r(:,t),r(:,mPowerTotal))
        plot(r(:,t),accPower)
        hold off
        legend('Electrical Power','Car-level Mechanical Power', ...
            'Sum of Motor Mechanical Powers','Total Acc Power (incl. losses)')

        figure
        plot(r(:,t),r(:,accV))
        hold all
        plot(r(:,t),r(:,backEMFF)); 
        plot(r(:,t),r(:,backEMFR)); 
        hold off

    %     figure
    %     plot(r(:,t),eff)

        figure
        plot(r(:,t),r(:,ePowerFL)); 
        hold all
        plot(r(:,t),r(:,ePowerRL)); 
        hold off

        %sanity check
        FzCheck = r(:,FzFL) + r(:,FzFR) + r(:,FzRL) + r(:,FzRR) - aeroDownTotal;

        figure
        plot(r(:,t),FzCheck)

        figure
        plot(r(:,t),r(:,FzFL))
        hold all
        plot(r(:,t),r(:,FzRL))
        plot(r(:,t),r(:,FzFR))
        plot(r(:,t),r(:,FzRR))
        legend('FzFL','FzRL','FzFR','FzRR')

        figure
        plot(r(:,t),aggkWh)
        title('Agg kWh')


        %% Histograms

        figure
        histogram(mcCurrentFL)
        title('Motor Controller Current FL')

        figure
        histogram(mcCurrentFR)
        title('Motor Controller Current FR')

        figure
        histogram(mcCurrentRL)
        title('Motor Controller Current RL')

        figure
        histogram(mcCurrentRR)
        title('Motor Controller Current RR')

        figure
        histogram(motorTorqueFL)
        title('Motor Torque FL')

        figure
        histogram(motorTorqueFR)
        title('Motor Torque FR')

        figure
        histogram(motorTorqueRL)
        title('Motor Torque RL')

        figure
        histogram(motorTorqueRR)
        title('Motor Torque RR')

        keyboard
    end
end