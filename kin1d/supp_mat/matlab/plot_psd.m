clear all


% plot PSD at three levels of filiq


load ./outQcldb2.dat
load ./outNcldb2.dat


i = 1;
for x = 5951:50:6000;
    
    z(1:50,i) = outQcldb2(x:x+49,3)*10.^-3;
    qi(1:50,i) = outQcldb2(x:x+49,10)*10.^-5;
    qil(1:50,i) = outQcldb2(x:x+49,9)*10.^-5;
    time(1:50,i) = outQcldb2(x:x+49,2);
    temp(1:50,i) = outQcldb2(x:x+49,6);
    qr(1:50,i) = outQcldb2(x:x+49,7)*10.^-5;
    qg(1:50,i) = outQcldb2(x:x+49,8)*10.^-5;
    bg(1:50,i) = outNcldb2(x:x+49,6)*10.^-8;    
    ns(1:50,i) = outNcldb2(x:x+49,5);
    nr(1:50,i) = outNcldb2(x:x+49,4);
    lam(1:50,i) = outNcldb2(x:x+49,11);
    mu(1:50,i) = outNcldb2(x:x+49,10);
    
%     i = i+1;
end
frim = qg./(qi-qil);
fliq = qil./qi;
rhorim = qg./bg;


