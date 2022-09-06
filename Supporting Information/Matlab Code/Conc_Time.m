function conc=Conc_Time(time,data,species)
L = length(data.Time);
time_stamp = 1;
for i=2:L
    if data.Time(i)>time & data.Time(i-1)<time
        time_stamp = i;
    end
end
m = (data.Conc.(species)(time_stamp) - data.Conc.(species)(time_stamp - 1))/(data.Time(time_stamp) - data.Time(time_stamp - 1));
conc = data.Conc.(species)(time_stamp-1) + (time-data.Time(time_stamp-1))*m;
end
