xin = importdata('data/x.dat');
x = complex(xin(:,1),xin(:,2));

yin = importdata('data/y.dat');
Y = complex(yin(:,1),yin(:,2));

X = fft(x,length(Y));

figure(1); plot(abs(x));
figure(2); stem(abs(X));
figure(3); stem(abs(Y));