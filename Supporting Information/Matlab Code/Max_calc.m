function [E]=Max_calc(data,species)
step = 0.0000001;
max1 = max(data.Time);
x = linspace(0+step,max1-step);
L1 = length(species);
y = [];
for j=1:length(x)
    conc1 = 0;
    conc2 = 0;
    for i=1:L1
        conc1 = conc1 + Conc_Time(x(j),data,species{i});
    end
    y = [y, conc1];
            
end
    
E = max(y);

end
