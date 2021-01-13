 function [Aclosed, Bclosed, Cclosed, Dclosed] = closedloop_from_plantplusctrler(Ap, Bp, Cp, Dp, Ac, Bc1,Bc2,Cc, Dc1, Dc2, yp2yclosed)
%[Ap, Bp, Cp, Dp] is the state-space model of the plant
%[Ac, Bc1, Bc2, Cc, Dc1, Dc2] is the state-space model of a generic controller
% yp2yclosed is a matrix defining the relation between the plant output(yp) and the closed-loop output(yclosed)
% yclosed = yp2yclosed*yp;
%yp2yclosed can be set to identity if the plant output and the closed-loop output are the same

Z = eye(size(Dc1,1))- Dc1*Dp;
Aclosed = [Ap+Bp*(Z\Dc1*Cp) , Bp*(Z\Cc) ; Bc1*(eye(size(Dp,1))+Dp*(Z\Dc1))*Cp , Ac+Bc1*Dp*(Z\Cc)];
Bclosed = [Bp*(Z\Dc2) ; Bc2+Bc1*Dp*(Z\Dc2)];

%%If the output of the closed-loop system is to be taken as the same as the output of the plant
Cclosed_large = [(eye(size(Dp,1))+ Dp*(Z\Dc1))*Cp , Dp*(Z\Cc)] ;
Dclosed_large = Dp*(Z\Dc2) ;

%%If the output of the closed-loop system is to be taken as a part of the plant output (or say some linear combination of the plant's output variables)
Cclosed = yp2yclosed*Cclosed_large;
Dclosed = yp2yclosed*Dclosed_large;