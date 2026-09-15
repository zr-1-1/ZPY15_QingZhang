function dy=atk_calibrated_dynamics(t,y,c)
% Empirical inertial J2 adapter fitted to saved ATK states, SI, t since epoch.
% Keep model_rv and its existing libraries unchanged.
tau=t/86400; p=c.PolePolynomial*[1;tau;tau^2];
axis=[p;sqrt(1-sum(p.^2))]; r=y(1:3); n=norm(r); z=dot(r,axis);
acc=-c.Mu/n^3*r+1.5*c.Mu*c.J2*c.GravityRadius^2/n^5*((5*z^2/n^2-1)*r-2*z*axis);
dy=[y(4:6);acc];
end
