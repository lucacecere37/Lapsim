function [Results] = sweepAccel(varargin)

%sweepAccel('var1','ptR.gr','var1Vals', 4:0.5:8,'var2','tire.longMu','var2Vals',[0.8:.1:1.6])
%% Parse Inputs

p = inputParser;
p.addOptional('var1',[])
p.addOptional('var2',[])
p.addOptional('var1Vals',[])
p.addOptional('var2Vals',[])
p.parse(varargin{:});

var1 = p.Results.var1;
var2 = p.Results.var2;
var1Vals = p.Results.var1Vals;
var2Vals = p.Results.var2Vals;

if(isempty(var1))
    error('No argument found for swept variable.')
end

numVars = ~isempty(var1) + ~isempty(var2);

%keyboard

%% Sim Setup

sim.dt = 0.001;
sim.g = 9.81;
sim.wtErrorThresh = 0.1; %0.05
sim.vErrorThresh = 0.01;
sim.aErrorThresh = 0.05;

assignVarColumns;

%% Sweep through variables

if(numVars == 1)
    
    Results = zeros(length(var1Vals),1);

    for kk = 1:length(var1Vals)
        disp(['Simulation running with ' var1 ' = ' num2str(var1Vals(kk)) ' ...'])
        
        car = WR217e;
        eval(['car.' var1 '=var1Vals(kk);']);
        car = initializeCar(car);
        tic
        [rF, i] = simStraight('sim', sim, 'car', car);%, 'debug', true, 'plots', true);
        toc
        Results(kk) = rF(i-1,t);
    end

    disp('Sweep complete.')

    figure
    plot(var1Vals,Results,'o-');
    title('Swept Var Results')
    
elseif(numVars ==2)
    Results = zeros(length(var1Vals),length(var2Vals));

    for kk = 1:length(var1Vals)
        for jj = 1:length(var2Vals)
            disp(['Simulation running with ' var1 ' = ' num2str(var1Vals(kk)) ' and ' var2 ' = ' num2str(var2Vals(jj)) ' ...'])

            car = WR217e;
            eval(['car.' var1 '=var1Vals(kk);']);
            eval(['car.' var2 '=var2Vals(jj);']);
            car = initializeCar(car);
            tic
            [rF, i] = simStraight('sim', sim, 'car', car);%, 'debug', true, 'plots', true);
            toc
            Results(kk,jj) = rF(i-1,t);
        end
    end

    disp('Sweep complete.')
    
    var1Vals = var1Vals';
    for ll = 1:length(var2Vals)
        var1Vals(:,ll) = var1Vals(:,1);
    end
    for ll = 1:length(var1Vals(:,1))
        var2Vals(ll,:) = var2Vals(1,:);
    end

    figure
    %surf(var1Vals,var2Vals,Results);
    surf(var1Vals,var2Vals,Results);
    title('Swept Var Results')
    
    [~, minIndCol] = min(Results);
    for ll = 1:length(minIndCol)
        resX(ll) = var2Vals(1,ll);
        resY(ll) = var1Vals(minIndCol(ll),ll);
    end
    figure
    plot(resX,resY)
    
    figure
    for jj = 1:length(var2Vals(1,:))
        hold all 
        plot(var1Vals(:,jj),Results(:,jj));
    end
    
end

keyboard