function [mid1, mid2]=sbb_branch(z_best,a,b)
if z_best>a &&z_best<b
    mid1=a+8*(z_best-a)/9;
    mid2= z_best+(b-z_best)/9;
elseif z_best==a
    mid1 = a+(b-a)/9;
    mid2= a+4*(b-a)/9;
else
    mid1 = a+4*(b-a)/9;
    mid2= a+8*(b-a)/9;
end
end