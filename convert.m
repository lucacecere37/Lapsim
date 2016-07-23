function [factor] = convert(arg1, arg2)

factor = [];

if(strcmp(arg1, 'rpm'))
    if(strcmp(arg2, 'rad/s') || strcmp(arg2, 'rad/sec'))
        factor = (2*pi)/60;    
    end
    return
end

if(strcmp(arg1, 'rad/s') || strcmp(arg1, 'rad/sec'))
    if(strcmp(arg2, 'rpm'))
        factor = 60/(2*pi);  
    end
    return
end

if(strcmp(arg1, 'mph'))
    if(strcmp(arg2, 'm/s'))
        factor = 0.447;
    elseif(strcmp(arg2,'kph'))
        factor = 1.609344;
    end
    return
end

if(strcmp(arg1, 'm/s'))
    if(strcmp(arg2, 'mph'))
        factor = 2.2369;
    end
    return
end

if(strcmp(arg1, 'deg'))
    if(strcmp(arg2, 'rad'))
        factor = (2*pi)/360;
    end
    return
end

if(strcmp(arg1, 'J'))
    if(strcmp(arg2, 'kWh'))
        factor = 2.778e-7;
    end
    return
end
        
if(isempty(factor))
    factor = 1;
    return
end


end