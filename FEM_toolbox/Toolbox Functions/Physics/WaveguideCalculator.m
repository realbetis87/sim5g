function [Beta,Fc] = WaveguideCalculator(a,b,frequency,mode)
    m=0;n=0;ElectromagneticConstants;k0=2*pi*frequency/c0;
    switch mode
        case 'TE10',m=1;n=0;
        case 'TE20',m=2;n=0;
        case 'TE30',m=3;n=0;
        case 'TE11',m=1;n=1;
        case 'TE21',m=2;n=1;
        case 'TE32',m=3;n=2;
    end
    kc=sqrt((m*pi/a)^2+(n*pi/b)^2);Fc=kc/(2*pi*c0);Beta=sqrt(k0^2-kc^2);
end