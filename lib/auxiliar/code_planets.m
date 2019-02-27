
planet_list

ind.planets = zeros(length(planet_fb),1);

for i  = 1:length(planet_fb)
    
    [truefalse,index] = ismember(planet_fb(i),planets);
    
    if truefalse
        
        ind.planets(i)  = index; 
        
    else
        
        error('The selected planet is not available')
        
    end
    
    
end