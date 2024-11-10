function [sigma] = GrapheneConductivity(Temp,gs,mc,frequency)
kb=1.3806503e-23;
e_charge=1.60217646e-19;
h=6.626068e-34;
hb=1.0546E-34;
heVs=4.135667516e-15;
Gs=2*pi*(gs)/(heVs);
Mc=mc*e_charge;
sigma = -1i*e_charge^2*kb*Temp*(Mc/(kb*Temp)+2*log(exp(-Mc/(kb*Temp))+1))/(pi*hb^2*(2*pi*frequency-1i*2*Gs));
end

