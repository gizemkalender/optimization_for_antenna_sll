using PyPlot
 N=24;
d=0.5;
alpha=-90:0.1:90;
Ed=abs(sin(pi*d*N*(sind(alpha)))./sin(pi*d*sind(alpha)))/N
Edl=20*log10(Ed)
subplot(4,1,1);
PyPlot.plot(alpha,Edl)
PyPlot.title("Radiation pattern dB");
subplot(4,1,2)
Ed_rad=abs(sin(pi*d*N*(sind(alpha)))./sin(pi*d*sind(alpha)))/(N)
PyPlot.plot(alpha,Ed_rad)
PyPlot.title("Normalized Radiation pattern dB");
subplot(4,1,3)
alpha_rad=-pi:0.01:pi
GAIN=abs(sind(pi*d*N*(sin(alpha_rad)))./sind(pi*d*sin(alpha_rad))).^2/N
GAIN_rad_db=10*log10(GAIN)
PyPlot.plot(alpha_rad,GAIN_rad_db)
PyPlot.title("Gain dB");
subplot(4,1,4)
dis=0.001:0.001:1
# alpha_rad=pi/2
GAIN=abs(sin(pi*dis*N) ./ sin(pi*dis)).^2 /N
# GAIN_rad_db=10*log10(GAIN)
PyPlot.plot(dis,GAIN)


#
# d=0.5;
# theta=-90:0.1:90;
# an=1
# AF=abs(an*sin(pi*d*sind(theta))+an*sin(3*pi*d*sind(theta))+an*sin(5*pi*d*sind(theta))+an*sin(7*pi*d*sind(theta))
# +an*sin(9*pi*d*sind(theta))+an*sin(11*pi*d*sind(theta)))
# AF_db=20*log10(AF)
# subplot(2,1,2);
# PyPlot.plot(theta,AF_db)
# # y=theta*0-20
# # PyPlot.plot(theta,y)
# # z=theta*0
# # PyPlot.plot(theta,z)
