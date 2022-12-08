function dYdt = ode_sys(t,Y,l,K,y1,y2,m0)

C=Y(3)-Y(1)*Y(1);
Cl=Y(4)-2*Y(2)*Y(1);
Cdl=Y(4)-Y(2)*Y(2); 
G=G_scalar(l,Y(1),K);
G1=G1_scalar(Y(1),K);
G2=G2_scalar(Y(1),K);
G_con=G_scalar(l,m0,K);
G1_con=(G1_scalar(m0,K));
G2_con=(G2_scalar(m0,K));

Gl=G_scalar(l,Y(2),K);
G1l=G1_scalar(Y(2),K);
G2l=G2_scalar(Y(2),K);

GE=G_scalar(l,Y(3),K);
G1E=G1_scalar(Y(3),K);
G2E=G2_scalar(Y(3),K);

GEl=G_scalar(l,Y(4),K);
G1El=G1_scalar(Y(4),K);
G2El=G2_scalar(Y(4),K);
ytot=y(y1,y2,l);


% the equations are m,m_\lambda, E, E_lambda

dYdt =   [ -C*G_con(1)*(ytot(1)-G(1));         %C*G^T*\gammainv*(y-G(m))
      
          -(-C*G_con(1)*Gl(1)+Cl*G_con(1)*(ytot(1)-G(1)))+...                                     %c*at*(b*integ-A*m)+cl*(b-A*m)
          -(C*(G1_con(1)-G2_con(1))*(ytot(1)-G(1))+C*G_con(1)*((y1(1)-y2(1))-(G1(1)-G2(1))));     %+c*(b-A*m)+c*(bl-Al*m)  integ=0

          2*C*G_con(1)*(ytot(1)*Y(1)-GE(1));

          ( 2*C*G_con(1)*(ytot(1)-GEl(1))+2*Cdl*(G_con(1))*(ytot(1)*Y(1)-GE(1)))+... % 2*c*at*(b*m-A*El)+2*cl*at*(b*m-A*E)
          ( 2*C*(G1_con(1)-G2_con(1))*(ytot(1)*Y(1)-GE(1)) +  2*C*G_con(1)*((y1(1)-y2(1))*Y(1)-(G1E(1)-G2E(1)))) ]; 
end