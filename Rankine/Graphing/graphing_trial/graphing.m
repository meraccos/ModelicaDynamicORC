data = readtable("data.csv");
A = data{2,:};
x = linspace(1, length(A), length(A));
plot(x(2:end), A(2:end))