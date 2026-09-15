function orbital_elements = rv2coe(r,v)
% ECI 位置速度 -> [a,e,i,Omega,omega,M]，m、m/s、rad；支持椭圆轨道。
% 有向 atan2 保留径向运动象限；退化轨道返回可重建同一状态的根数。
global mu
validateattributes(r,{'numeric'},{'real','finite','numel',3});
validateattributes(v,{'numeric'},{'real','finite','numel',3});
r=r(:); v=v(:); radius=norm(r); h=cross(r,v); hn=norm(h);
assert(radius>0 && hn>0,'位置与角动量必须非零。');
energy=dot(v,v)/2-mu/radius;
assert(energy<0,'当前函数只支持有界椭圆轨道。');
a=-mu/(2*energy); ev=cross(v,h)/mu-r/radius; e=norm(ev);
assert(e<1,'当前函数只支持偏心率小于 1 的轨道。');
i=atan2(norm(h(1:2)),h(3));
n=cross([0;0;1],h); nn=norm(n); tol=1e-12;
if nn/hn<tol
    Omega=0;
    if e<tol
        omega=0; nu=atan2(sign(h(3))*r(2),r(1));
    else
        omega=atan2(sign(h(3))*ev(2),ev(1));
        nu=oriented_angle(ev,r,h);
    end
else
    Omega=atan2(n(2),n(1));
    if e<tol
        omega=0; nu=oriented_angle(n,r,h);
    else
        omega=oriented_angle(n,ev,h); nu=oriented_angle(ev,r,h);
    end
end
eccentric=atan2(sqrt(1-e^2)*sin(nu),e+cos(nu));
M=eccentric-e*sin(eccentric);
orbital_elements=[a,e,i,mod(Omega,2*pi),mod(omega,2*pi),mod(M,2*pi)];
end

function angle = oriented_angle(a,b,normal)
angle=atan2(dot(cross(a,b),normal)/norm(normal),dot(a,b));
end
