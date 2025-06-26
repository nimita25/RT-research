function [Dij,totalNumBixel] = aggregatedij(ptid_1,dij_cell,gantryAngles,scale,bw,lss)
%% Aggregated into a large Dij
Dij = [];
% if ~isdeployed
    totalNumBixel = 0;
    sampleAngle = downsample(gantryAngles,scale);
    Id_sampleAgle = sampleAngle/10 + 1;
    for j = Id_sampleAgle
        %             Dij = cat(2,Dij,dij_cell{j});
        totalNumBixel = totalNumBixel + dij_cell{1,j}.totalNumOfBixels;
        Dij = cat(2,Dij,dij_cell{1,j}.physicalDose{1});
    end
    save([ptid_1 '\','51426287_Dij_proton_',num2str(bw),num2str(lss),'_Spec_energy_NumAngle_',num2str(numel(gantryAngles)/scale),'.mat'],'Dij','-v7.3')
% end
end