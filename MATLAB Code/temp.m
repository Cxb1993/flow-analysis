% for i = 1:23
%     figure;
%     pressure = (i * 100) + 100;
%     if (pressure < 1000)
%         pressureStr = strcat('0', num2str(pressure));
%     else
%         pressureStr = num2str(pressure);
%     end
    
%     master(strcat('tracking\data\', pressureStr, ' mbar 0.0_MMStack_Pos0.ome.mat'), strcat('tracking\images\', pressureStr, ' mbar 0.0_MMStack_Pos0.ome.tif'), pressureStr);
    master(strcat('tracking\data\0000 mbar 0003 +0_MMStack_Pos0.ome.mat'), strcat('tracking\images\0000 mbar 0003 +0_MMStack_Pos0.ome.tif'), 1);  
    drawnow;
% end