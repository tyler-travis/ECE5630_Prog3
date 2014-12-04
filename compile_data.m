xin = importdata('data/x.dat');
x = complex(xin(:,1),xin(:,2));

yin = importdata('data/y.dat');
y = complex(yin(:,1),yin(:,2));

X = fft(x,216);

figure(1); plot(x);
figure(2); plot(abs(X));
figure(3); plot(abs(y));