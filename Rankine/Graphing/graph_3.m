rawdata = readtable("data.csv");
data = rawdata{:,:};

n = length(data(:,1));
m = length(data(1,:));

x = linspace(1, ((m-1)/3), ((m-1)/3));

for i = 1:n
    plot(x, data(i, 2:((m+2)/3)), x, data(i, ((m+5)/3):((2*m+1)/3)), x, data(i, ((2*m+4)/3):m))
    
    
    set(gcf,'position',[200,10,1200,900])
    legend(strcat(string(data(i,1)), ' = time'))
    pause(0.001)
end

