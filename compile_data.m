xin = importdata('data/x.dat');
x = complex(xin(:,1),xin(:,2));

yin = importdata('data/y.dat');
Y = complex(yin(:,1),yin(:,2));

X = fft(x,length(Y));
ny = (0:(length(Y)-1))/length(Y);
nx = (0:(length(X)-1))/length(X);

figure(1); plot(abs(x));
figure(2); plot(nx,20*log10(abs(X))); title('Matlab fft');
figure(3); plot(ny,20*log10(abs(Y))); title('Radix-6 fft');