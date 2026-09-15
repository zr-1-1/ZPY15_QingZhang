function test_rv2coe_roundtrip()
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
global mu
mu=3.986004418e14;
% 径向速度正负、极轨道、正反向圆/椭圆赤道轨道。
for inc=[0,pi/2,98*pi/180,pi]
    for ecc=[0,.02]
        for anomaly=[.4,2,4,5.5]
            state=orb_elements2rv([7200e3,ecc,inc,1.3,.7,anomaly]);
            rebuilt=orb_elements2rv(rv2coe(state(1:3),state(4:6)));
            assert(norm(rebuilt(1:3)-state(1:3))<1e-5);
            assert(norm(rebuilt(4:6)-state(4:6))<1e-8);
        end
    end
end
t=readtable(fullfile(root,'data','debris_orbits.csv'));
for id=1:height(t)
    state=t{id,7:12}; rebuilt=orb_elements2rv(rv2coe(state(1:3),state(4:6)));
    assert(norm(rebuilt(1:3)-state(1:3))<1e-5);
    assert(norm(rebuilt(4:6)-state(4:6))<1e-8);
end
fprintf('test_rv2coe_roundtrip passed (32 synthetic and %d CSV states).\n',height(t));
end
