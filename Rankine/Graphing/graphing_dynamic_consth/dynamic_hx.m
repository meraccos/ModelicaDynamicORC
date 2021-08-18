rawdata = readtable("data.csv");
data = rawdata{:,:};

n = length(data(:,1));
m = length(data(1,:));

x = linspace(1, ((m-1)/2), ((m-1)/2));

for i = 1:n
    plot(x, data(i, 2:((m+1)/2)), x, data(i, ((m+3)/2):m))
    
    
    set(gcf,'position',[200,10,1200,900])
    legend(strcat(string(data(i,1)), ' = time'))
    pause(0.0001)
end

