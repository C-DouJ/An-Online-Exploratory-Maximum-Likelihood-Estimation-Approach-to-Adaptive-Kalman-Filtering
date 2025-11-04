function [xkk,Pkk,Pk1k]=kf_std(xkk,Pkk,F,H,B,u,G,z,Q,R)

%%%%%Time update
xk1k=F*xkk+B*u;

Pk1k=F*Pkk*F'+G*Q*G';

%%%%%Measurement update
Pzzk1k=H*Pk1k*H'+R;

Pxzk1k=Pk1k*H';

Kk=Pxzk1k*inv(Pzzk1k);

xkk=xk1k+Kk*(z-H*xk1k);

Pkk=Pk1k-Kk*H*Pk1k;