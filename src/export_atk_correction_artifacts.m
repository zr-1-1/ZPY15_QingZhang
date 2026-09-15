function export_atk_correction_artifacts()
load('data/atk_correction/corrected_plan.mat','plan','calibration');
load('data/atk_correction/adapter_verification.mat','verification','auditTable');
load('data/atk_correction/export_crosscheck.mat','crosscheck','roundingTable');
rows=zeros(3,7);
for ship=1:3,rows(ship,:)=[ship,plan.Ships(ship).InitialState'];end
% Full precision for reconstructing the XML's actual initial Cartesian state.
fid=fopen('data/atk_correction/initial_states_J2000.csv','w');
fprintf(fid,'Mother,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps\n');
for ship=1:3,fprintf(fid,'%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n',rows(ship,:));end
fclose(fid);
pulses=readtable('data/atk_correction/maneuvers.csv');
pulses.UTC=plan.EpochUTC+seconds(pulses.TimeSeconds);pulses.UTC.Format='yyyy-MM-dd HH:mm:ss';
pulses.Norm_mps=vecnorm(pulses{:,{'DeltaVx_mps','DeltaVy_mps','DeltaVz_mps'}},2,2);
writetable(pulses,'data/atk_correction/maneuvers_with_UTC.csv');
summary=struct('NativeATKVerified',false,'Score',verification.Score, ...
    'EligiblePerShip',verification.EligiblePerShip,'MinimumAltitudes_m',crosscheck.MinimumAltitudes, ...
    'MaxIndependentAltitudeDifference_m',crosscheck.MaxAltitudeDifference_m, ...
    'MaxOriginalATKTrajectoryPositionDifference_m',auditTable.MaxPositionDifference_m', ...
    'CalibrationPolePolynomial',calibration.PolePolynomial, ...
    'OriginalATKClearanceSetsReproduced',crosscheck.OriginalATKClearanceSetsReproduced, ...
    'RoundingSensitivity',table2struct(roundingTable));
fid=fopen('data/atk_correction/summary.json','w');fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));fclose(fid);
f=figure('Visible','off','Color','w','Position',[50,50,1100,850]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;plot((0:86400)/3600,verification.Altitudes/1000,'LineWidth',1.1);hold on;
yline(200,'r--','200 km');xlabel('Hours since 2030-11-14 08:00 UTC');ylabel('Geocentric altitude (km)');
legend('Mother1','Mother2','Mother3','Location','best');grid on;
title('Corrected plan: empirical ATK adapter; native ATK rerun pending');
nexttile;
for ship=1:3
    counts=arrayfun(@(t)sum(verification.Events.Mother==ship & verification.Events.TimeSeconds<=t),0:60:86400);
    stairs((0:60:86400)/3600,counts,'LineWidth',1.3);hold on;
end
xlabel('Hours since epoch');ylabel('Distinct cleared debris');grid on;
legend(arrayfun(@(s)sprintf('Mother%d: %d',s,verification.EligiblePerShip(s)),1:3,'UniformOutput',false),'Location','best');
exportgraphics(f,'docs/ATK修正方案复核.png','Resolution',160);close(f);
end
