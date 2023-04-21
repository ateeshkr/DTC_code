clear all;
%% Mass of graphene
mg=10.*.2695.*2.32e-16;
%% Frequencies of two graphene drum
OmegaG1 =2.*pi.*2.17e6;
OmegaG2 =2.*pi.*2.23e6;
%% Parametric drive frequency
OmegaF=OmegaG1+OmegaG2;
%% Damping
gammag=2.*pi.*22037;
%% Duffing coefficient and nonlinear damping coefficient
beta=5.8e15;
eta=1.8e10; 
%% Effective coupling
alpha=1.*4.*pi.^2.*2.358e-5;%1.8
%% Parametric drive strength
% epsilon=0.01;%0
epsilon=input('Parametric drive strength (provide 0 or 0.01) =  ')
%% Thermalnoise strength
Fth=5e-2;
%% Differential equations
f = @(t,y) [y(2); -gammag.*y(2)-(eta./mg).*y(1).*y(1).*y(2)-(beta./mg).*y(1)^3-(OmegaG1.*OmegaG1+(epsilon/mg).*cos(OmegaF.*t)).*y(1)+(alpha./mg).*y(3);
            y(4); -gammag.*y(4)-(eta./mg).*y(3).*y(3).*y(4)-(beta./mg).*y(3)^3-(OmegaG2.*OmegaG2+(epsilon/mg).*cos(OmegaF.*t+0*pi/4)).*y(3)+(alpha./mg).*y(1)]; % define function f(t,y)

%% Initial conditions
y0 = [1e-15;0;1e-15;0]; 

%% time steps
t0 = 0;T = 5e-4;  
N = 1.*10.^5;

%% Averaging loop
index=1:1:50;
%% Noise increment loop 
mul=0:0.1:2;
%% Data matrix initialization
all_Av_peak=zeros([length(index) length(mul)]);
Av_peak=zeros([1 length(mul)]);
av_tdarray=zeros([length(mul) (N)/2+1]);

for j=1:length(index)
    tdarray=zeros([length(mul) (N)/2+1]);
    mean_noise_array=zeros([length(mul) 1]);
    Av_peak=zeros([length(mul) 1]);
    
    for i=1:length(mul)
      epsilon_p=mul(i);
      [ts,ys,freq,P1] = rk4_sto(f,[t0,T],y0,N,epsilon_p,Fth);
      tdarray(i,:)=P1;
      if i==1
          [val, idx] = max(tdarray(i,:));    
      end 
    end
    
    for i=1:length(mul)
        Av_peak(i)=tdarray(i,idx);
    end
    
     all_Av_peak(j,:)=Av_peak;
     av_tdarray=av_tdarray+tdarray;

end
av_tdarray=av_tdarray/length(index);
peak=mean(all_Av_peak);

%% output plot 
figure(3)
plot(mul,peak)
set(gca,'YScale','log')

Mat1=av_tdarray(:,1001:1201);
frequency=freq(:,1001:1201);

disp('Collect mul (noise strength), frequency and Mat1 (frequency along x axis and noise strength along y axis) ')

%% RK4 SDE function with fourier transform of time series data
function [ts,ys,freq,P1] = rk4_sto(f,tv,y0,N,epsilon_p,Fth)

  t0 = tv(1); T = tv(2);
  dt = (T-t0)/N;                    % stepsize h
  ts = zeros(N+1,1); ys = zeros(N+1,length(y0));
  t = t0; y = y0;                  % initial point
  ts(1) = t; ys(1,:) = y';
  n_matrix=[0 1 0 1];
  Fth = Fth.*n_matrix;
  Fadd = epsilon_p.*n_matrix;
%   st=sqrt(dt).*Fadd.*randn(N+1,1)*1;
  st=sqrt(dt).*Fth.*randn(N+1,1)+sqrt(dt).*Fadd.*randn(N+1,1);

  for i=1:N
      k1 = dt*f(t,y);
      k2 = dt*f(t+0.5*dt,y+0.5*k1+0.5*st(i,:));
      k3 = dt*f(t+0.5*dt,y+0.5*k2+0.5*st(i,:));
      k4 = dt*f(t+dt,y+k3+st(i,:));
%       y = y + dt*f(t,y)+ st(i,:)';
      y = y + (k1+2*k2+2*k3+k4)/6+ st(i,:)';
      t = t + dt;
    
      ts(i+1) = t; ys(i+1,:) = y';   % store y(1),y(2) in row of array ys
%       var(i)=k1(2);%st(i,2);
  end
    % frequency
    Fs=1/dt;
    L=(length(ts)-1);
    freq = Fs*(0:(L/2))/L;
    
    % PSD
    ffty1= fft(ys(:,1));
    P2 = abs(ffty1/L);%/CG;
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
end