
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;

%==============================
% Forward Signals
%==============================
openfig('Example51_Sensor1_signal.fig'); 
set(gca,'FontSize',13)
saveas(gcf, 'Example51_Sensor1_signal.fig');
saveas(gcf, 'Example51_Sensor1_signal', 'epsc');

openfig('Example51_Sensor2_signal.fig'); 
set(gca,'FontSize',13)
saveas(gcf, 'Example51_Sensor2_signal.fig');
saveas(gcf, 'Example51_Sensor2_signal', 'epsc');

openfig('Example51_Sensor3_signal.fig'); 
set(gca,'FontSize',13)
saveas(gcf, 'Example51_Sensor3_signal.fig');
saveas(gcf, 'Example51_Sensor3_signal', 'epsc');
