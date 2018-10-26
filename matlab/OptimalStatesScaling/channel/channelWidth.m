function a = channelWidth(Le,Ra,T0,S0,dSdz,dTdz)

gamma = Le*dSdz;
beta = sqrt(Ra*gamma);

C = beta^2*(S0-T0)/(dTdz - Le*dSdz);
a = (1/beta)*acos(1/(C+1));


C = beta^2*(S0+T0)/(dTdz + Le*dSdz);
a =(1/beta)*acos(1/(1-C));

end