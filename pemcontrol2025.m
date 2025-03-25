clear
close all
% clc
%%


%%

wavset=[760.153,780.124,800.4016,821.467,840.439,859.362];
% wavset=[751,778,786,810,813,880,912,940]
slopevals=[];%Creating an array for the slopes generateed in the for loop
interceptvals=[];%Creating an array for the y intercepts generateed in the for loop
figure;
hold on;

ydataT=cell(1,6);
xdataT=cell(1,6);
for i=1:size(wavset, 2)
    elk="7_11_"+num2str(wavset(1,i))+" nm.xlsx";
   % elk="voltagefor"+num2str(wavset(1,i))+".xlsx";
    data=xlsread(elk);
    xdata = data(:,1); % the primes are to make these into a column array, as the fit function only wants column arrays
    ydata = data(:,2);
    
    % For the equation, the syntax can matter heavily. A linear fit would be "A*x + B", as we always use x as the data variable
    % THe fitting parameters are the Lower Bounds, the Upper Bounds, and the Starting Point. For an example, assume the slope in the linear function is around 3. Thus, since we think it's around 3, we set the starting point as 3, and then the upper bound at 5, and lower bound at 1, so that the fitting method can extend around that point to give us a precise value.
    fo=fitoptions('Method','NonlinearLeastSquares',...
                'Lower',([4, -5]),....
                'Upper',([5.5, 5]),....
                'StartPoint',([5, 0]));
    ft=fittype("A*x+B",...
                'option',fo,'independent','x');
    
    [c, gof]=fit(xdata, ydata,ft);
    % c is a structure variable that contains the information on what the fitting parameters of A and B are, and also holds statistic information on how it is plotted. We can either take the parameters out, or plot the fit
    slope = c.A;
    offset = c.B;
    slopevals=[slopevals,slope];
    interceptvals=[interceptvals,offset];
    % gof is parameters that show how good the fit is, such as R^2 etc. It's not often necessary.
    
    
    xspec=linspace(xdata(1,1), xdata(end,1), 400);
    yfit=slope.*xspec+offset;
    
    % plot(c)
    hold on;
    % subplot(5, 4, i); %gets a subplot of all the plots
    figure %gets individual plots
    hold on;
    plot(xdata, ydata);
    plot(xspec, yfit);
    legend(wavset(1,i)+"nm m= "+slope);

    plot(xdata, ydata)
    
    ydataT{1,i}=ydata;
    xdataT{1,i}=xdata;

end
%%
%Finding the dependence of wavelength on slope

hold on
fo=fitoptions('Method','NonlinearLeastSquares',...
            'Lower',([0,-5]),....
            'Upper',([5.5,5]),....
            'StartPoint',([0.05,-4.1]));
ft=fittype("A*(x)+B",...
            'option',fo,'independent','x');

[d, gof]=fit(wavset', slopevals',ft);

figure;
plot(d);
hold on;
plot(wavset, slopevals);
Aa=d.A;
Ba=d.B;
%Finding the y intercept depending on the wavelength / This is not included
%due to it being minute to the overall voltage
% fo=fitoptions('Method','NonlinearLeastSquares',...
%             'Lower',([0, -5, 700]),....
%             'Upper',([5.5, 5, 760]),....
%             'StartPoint',([0.5, 4.2, 740]));
% ft=fittype("A*(x+C)+B",...
%             'option',fo,'independent','x');
% 
% [f, gof]=fit(wavset', interceptvals',ft)
% 
% figure
% plot(f)
% hold on;
% plot(wavset, interceptvals)
% Ab=f.A;
% Bb=f.B;
% Cb=f.C;

%%
xspec=linspace(xdata(1,1), xdata(end,1), 400);
yfit=slopevals(1,1).*xspec+offset;
tfit=xspec*(Aa*760+Ba)+offset;
figure;
hold on;
plot(xspec,yfit);
plot(xspec,tfit);
 
%%
ziDAQ('connect',  '127.0.0.1', 8004, 5);
example_poll('dev3205');
ziDAQ('subscribe','/dev3205/demods/2')
ziDAQ('subscribe','/dev3205/demods/3')
result = ziDAQ('listNodes', '/dev3205/demods/', int64(0))
d = daq("mcc");
%Add channels
addoutput(d,"Board0","Ao0","Voltage");
d.Channels;
%%
% wavelengt=810;
% retarded=0.223;
% voltage=retarded*(Aa*wavelengt+Ba)+offset;
% besselval=besselj(2,2*pi*retarded)/besselj(4,2*pi*retarded)
% dcOutput = voltage;
%         write(d,dcOutput);
% pause(10)
%         data=ziDAQ('poll', 10, 0.01,[0]);
%         x2data = data.dev3205.demods(3).sample.x;
%         y2data = data.dev3205.demods(3).sample.y;
%         r2=sqrt(mean(x2data)^2+mean(y2data)^2);
%         x3data = data.dev3205.demods(4).sample.x;
%         y3data = data.dev3205.demods(4).sample.y;
%         r3=sqrt(mean(x3data)^2+mean(y3data)^2);
%         therval=r2/r3
% 
%  pause(30)
%  while (therval > besselval+0.03) || (therval < besselval-0.03)
%         dcOutput = voltage;
%         pause(3)
%         write(d,dcOutput);
% 
%         data=ziDAQ('poll', 1, 0.03,[0]);
%         x2data = data.dev3205.demods(3).sample.x;
%         y2data = data.dev3205.demods(3).sample.y;
%         r2=sqrt(mean(x2data)^2+mean(y2data)^2);
%         x3data = data.dev3205.demods(4).sample.x;
%         y3data = data.dev3205.demods(4).sample.y;
%         r3=sqrt(mean(x3data)^2+mean(y3data)^2);
%         therval=r2/r3;
%         if therval > besselval
%             voltage = voltage + 0.005;
%         elseif therval < besselval
%             voltage = voltage - 0.005;
%         end
%         pause(1)
%         % result = ziDAQ('listNodes','/dev3205/demods/' , int64(16));
%         % ziDAQ('poll', duration, 10 ms, [flags]);
%     end
%% While loop
therval=0;
wavelengt=820.508;
offset=-0.0553;
Aa=0.0056;
Ba=0.0567;
retardance=[0.368,0.375,0.383,0.390,0.4,0.45,0.55];
% retardance=[0.368]
voltagefor814=zeros(size(retardance))
avratio=zeros(size(retardance))
newvoltage=[1.57728758500000,1.60563823900000,1.62311225800000,1.67233954700000,1.68281745400000,1.88498316600000,2.30668836300000];
% for k=1:linspace(1,2,1)
for i=1:size(retardance, 2)
    voltage=newvoltage(1,i);
    besselval=besselj(2,2*pi*retardance(1,i))/besselj(4,2*pi*retardance(1,i))


    while (therval > besselval+0.03) || (therval < besselval-0.03)
        dcOutput = voltage;
   
        write(d,dcOutput);
     pause(3)
        data=ziDAQ('poll', 5, 0.03,[0]);
        x2data = data.dev3205.demods(3).sample.x;
        y2data = data.dev3205.demods(3).sample.y;
        r2=sqrt(mean(x2data)^2+mean(y2data)^2);
        x3data = data.dev3205.demods(4).sample.x;
        y3data = data.dev3205.demods(4).sample.y;
        r3=sqrt(mean(x3data)^2+mean(y3data)^2);
        therval=r2/r3
        if therval > besselval
            voltage = voltage + 0.008;
        elseif therval < besselval
            voltage = voltage - 0.008;
        end
        pause(1)
        % result = ziDAQ('listNodes','/dev3205/demods/' , int64(16));
        % ziDAQ('poll', duration, 10 ms, [flags]);
    end
    voltagefor814(1, i)=voltage
    avratio(1,i)=therval
end
% end
%% Test

retardance=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55]
avratio=zeros(size(retardance));
therval=0;
bval=0;
sevals=0;
wavelengt=800;
for i=1:size(retardance, 2)
 voltage=retardance(1,i)*(Aa*878.9+Ba)+offset;%slopevals(1,3)*retardance(1,i)+offset;
besselval=besselj(2,2*pi*retardance(1,i))/besselj(4,2*pi*retardance(1,i))
   dcOutput = voltage;
   write(d,dcOutput);
   pause(15)
   data=ziDAQ('poll', 5, 0.03,[0]);
        x2data = data.dev3205.demods(3).sample.x;
        y2data = data.dev3205.demods(3).sample.y;
        r2=sqrt(mean(x2data)^2+mean(y2data)^2);
        x3data = data.dev3205.demods(4).sample.x;
        y3data = data.dev3205.demods(4).sample.y;
        r3=sqrt(mean(x3data)^2+mean(y3data)^2);
        therval=r2/r3;
                x2se = std(x2data)/sqrt(length(x2data));
                y2se=std(y2data)/sqrt(length(y2data));
                x3se = std(x3data)/sqrt(length(x3data));
                y3se=std(y3data)/sqrt(length(y3data));
                r2se=sqrt(((mean(x2data)*x2se)^2/((mean(x2data))^2+mean(y2data)^2))+((mean(y2data)*y2se)^2/(mean(x2data)^2+(mean(y2data)^2))));
                r3se=sqrt(((mean(x3data)*x3se)^2/((mean(x3data))^2+mean(y3data)^2))+((mean(y3data)*y3se)^2/(mean(x3data)^2+(mean(y3data)^2))));
                thervalse=sqrt((r2se/r3)^2+((r2*r3se)/r3^2));
        avratio(1,i)=therval
        bval(1,i)=besselval
        sevals(1,i)=thervalse
end
%% Standard Error
 data=ziDAQ('poll', 1, 0.03,[0]);
        x2data = data.dev3205.demods(3).sample.x;
        y2data = data.dev3205.demods(3).sample.y;
        r2=sqrt(mean(x2data)^2+mean(y2data)^2);
        x3data = data.dev3205.demods(4).sample.x;
        y3data = data.dev3205.demods(4).sample.y;
        r3=sqrt(mean(x3data)^2+mean(y3data)^2);
        therval=r2/r3;
x2se = std(x2data)/sqrt(length(x2data));
y2se=std(y2data)/sqrt(length(y2data));
x3se = std(x3data)/sqrt(length(x3data));
y3se=std(y3data)/sqrt(length(y3data));
r2se=sqrt(((mean(x2data)*x2se)^2/((mean(x2data))^2+mean(y2data)^2))+((mean(y2data)*y2se)^2/(mean(x2data)^2+(mean(y2data)^2))));
r3se=sqrt(((mean(x3data)*x3se)^2/((mean(x3data))^2+mean(y3data)^2))+((mean(y3data)*y3se)^2/(mean(x3data)^2+(mean(y3data)^2))));
thervalse=sqrt((r2se/r3)^2+((r2*r3se)/r3^2));
%% Plot
wavelengt=878.8
xspec=linspace(xdata(1,1), xdata(end,1), 400);
bval=0;
for i=1:size(xspec,2)
besselval=besselj(2,2*pi*xspec(1,i))/besselj(4,2*pi*xspec(1,i));
bval(1,i)=besselval;
end
elk="7_11_"+wavelengt+" error nm.xlsx";
   % elk="voltagefor"+num2str(wavset(1,i))+".xlsx";
    data=xlsread(elk);
    xdata = data(:,1); % the primes are to make these into a column array, as the fit function only wants column arrays
    ydata = data(:,2);
    errordata=data(:,3);
figure;
hold on;
plot(xspec,bval,':')
errorbar(xdata,ydata,errordata,'.','Vertical')
title('Bessel Check at 878 nm')
xlabel('Retardance')
ylabel('I(2nd Harmonic)/I(4th Harmonic)')
legend('theoretical','PEM')
%%
stop(d)
ziDAQ('disconnectDevice', 'dev3205');