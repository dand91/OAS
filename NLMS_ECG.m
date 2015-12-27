
s = 5000;

load pulsedata.mat;

ekg = ekg(1:s,1);
ekgEchotemp = downsample(ekg(1:s,1),2);
delay = 5;
ekgEcho = vertcat(ekgEchotemp,ekgEchotemp);
ekgEcho = vertcat([zeros(1,delay)]',ekgEcho(1:end-delay));
n1 = randn(1,s) + Laplace(1, 1, s) + 0.1*ekgEcho';

%%

warning('off','all')

load pulsedata.mat;

s = 5000;

ekg = ekg(1:s,1);

a = 0.2;

sV = filter([1-a a-1],[1 a-1], ekg' + n1);
nV = filter([1-a a-1],[1 a-1], n1 ); 

imp = [1; zeros(49,1)];

plot( filter([1-a a-1],[1 a-1], imp))

%zplane([1-a a-1],[1 a-1]);
%plot(ekg);
%plot(conv(sV,0.5*ones(1,20)))
%plot(sV);
%plot(conv(sV,0.5*ones(1,10)))
%xlabel('Real');
%ylabel('Imaginary');
%axis([0 5000 -12 8 ]);

%%

%initiate variables
%for i = 10:50

u= n1;
d= nV;
un = [];
dn = [];

%size of filter one and two

M1 = 20;
M2 = 4;

%preallocate memory

w1=zeros(M1,1);
w2=zeros(M2,1);

%input signal length

N=length(u);

%make sure that u and d are colon vectors

u=u(:);
d=d(:);

%step size

mu1 = 1;
mu2 = 1;

%lms theoretical bound

%mu1 =(2/(M1*var(u)));
%mu2 =(2/(M2*var(u)));

%guard value for nlms

a1 = 1000;
a2 = 1000;

for n=M1:N
    
    %System identification
    
    uvec=u(n:-1:n-M1+1);
    e1(n)=d(n)-w1'*uvec;
    w1=w1+mu1/(a1+uvec'*uvec)*uvec*conj(e1(n));
    
    %Filtering
    
    un = filter(1,w1,nV(1:n));
    dn = filter(1,w1,sV(1:n));
    
    %Noise cancelation
    
    uvec2=un(end:-1:end-M2+1)';
    y2(n) = w2'*uvec2;
    e2(n) = dn(end)'-y2(n);
    w2=w2+mu2/(a2+uvec2'*uvec2)*uvec2*conj(e2(n));
    
end

%MSE2(i) = sum((ekg(1:length(e2(1500:end)))'-e2(1500:end)).^2)/length(e2(1500:end));

%end

%plot(10:50 ,MSE2(10:50));

plot(e2);
ylabel('Amplitude');
xlabel('Time');
hold on
plot(ekg)
zoom on
grid on
%legend('20','25','30');
axis([0 5000 -8 8]);