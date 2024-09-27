clear all; clc; close all;

folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
filename = "DMD-AUGDMD-Seq1.eps";
fullFilePath = fullfile(folderPath, filename);

F_DMD_AUG_Seq1 = load("f_dmd_aug_win500_r5_seq1_predictionWin250.mat");


F_DMD_STAN_Seq1 = load('f_dmd_standard_win500_r5_seq1_predictionWin250.mat');

myVideo = VideoWriter('myVideoFile','MPEG-4'); %open video file
myVideo.FrameRate = 100;  %can adjust this, 5 - 10 works well for me
open(myVideo)

f = figure;
f.Color = 'w';
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 16])
lin_width = 1.2;
lin_width_tick = 1;
font_size = 15;

t = tiledlayout("vertical");
t.TileSpacing = 'compact';
t.Padding = 'compact';

numSteps_approx = F_DMD_AUG_Seq1.t_extra_0;
numSteps_prediction = length(F_DMD_AUG_Seq1.t)-numSteps_approx;
axis([0 .15 -20 20 ])

h1 = animatedline('Color',"#4daf4a",'LineWidth',lin_width);
h2 = animatedline('Color',"#984ea3",'LineWidth',lin_width);
h3 = animatedline('Color',"#ff7f00",'LineWidth',lin_width);
ylabel(t, '$\dot{z}$(t) [mm/s]', 'Interpreter', 'latex','FontSize',font_size);
xlabel(t, 'Time [s]', 'Interpreter', 'latex','FontSize',font_size);
t = F_DMD_STAN_Seq1.t;


grid on

textHandle = text(0.05, -18, 'Measurement', 'Interpreter', 'latex','FontSize',font_size);

for step = 1:numSteps_approx
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
    if step == 1
        xline(F_DMD_STAN_Seq1.t(end)-F_DMD_STAN_Seq1.extraTime,'LineWidth',lin_width,...
 'Interpreter','latex','HandleVisibility','off');
        leg = legend('Measured','DMD','AugDMD','Orientation', 'Horizontal');
       % leg.Layout.Tile = 'north';
        leg.Interpreter = 'latex';
        leg.FontSize = font_size;
        leg.Direction = "reverse";
    end
    

    y1 = F_DMD_STAN_Seq1.F_AUG_1(step)*10;

    
  
    addpoints(h1,t(step),y1);


    pause(.03)
    drawnow

    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);

end
hold on
close(myVideo);

delete(textHandle); % Remove the text

textHandle = text(0.05, -18, 'Approximation', 'Interpreter', 'latex','FontSize',font_size);
plot(F_DMD_AUG_Seq1.t(1:numSteps_approx),F_DMD_STAN_Seq1.F_DMD_1(1:numSteps_approx)*10,'Color',"#984ea3",'LineWidth',lin_width,'HandleVisibility','off');
% frame = getframe(gcf); 
% writePauseFrames(myVideo, frame, 2, myVideo.FrameRate);
plot(F_DMD_STAN_Seq1.t(1:numSteps_approx),F_DMD_AUG_Seq1.F_DMD_1(1:numSteps_approx)*10,'Color',"#ff7f00",'LineWidth',lin_width,'HandleVisibility','off');

f = gcf;
exportgraphics(f,'Fig1.png','Resolution',300)
delete(textHandle); % Remove the text
% frame = getframe(gcf); %get frame
% writePauseFrames(myVideo, frame, 2, myVideo.FrameRate);



textHandle = text(0.14, -18, 'Prediction', 'Interpreter', 'latex','FontSize',font_size);
plot(F_DMD_AUG_Seq1.t(numSteps_approx:end),F_DMD_STAN_Seq1.F_DMD_1(numSteps_approx:end)*10,'Color',"#984ea3",'LineWidth',lin_width,'HandleVisibility','off');

% frame = getframe(gcf); %get frame
% writePauseFrames(myVideo, frame, 2, myVideo.FrameRate);


plot(F_DMD_STAN_Seq1.t(numSteps_approx:end),F_DMD_AUG_Seq1.F_DMD_1(numSteps_approx:end)*10,'Color',"#ff7f00",'LineWidth',lin_width,'HandleVisibility','off');

f = gcf;
exportgraphics(f,'Fig2.png','Resolution',300)

delete(textHandle); % Remove the text
textHandle = text(0.14, -18, 'Measurement', 'Interpreter', 'latex','FontSize',font_size);

myVideo2 = VideoWriter('myVideoFile2','MPEG-4'); %open video file
myVideo2.FrameRate = 100;  %can adjust this, 5 - 10 works well for me
open(myVideo2)


for step = numSteps_approx:numSteps_approx+numSteps_prediction
    y1 = F_DMD_STAN_Seq1.F_AUG_1(step)*10;

    
  
    addpoints(h1,t(step),y1);


    pause(.03)
    drawnow
    frame = getframe(gcf); %get frame
    writeVideo(myVideo2, frame);

end

%writePauseFrames(myVideo, frame, 2, myVideo.FrameRate);


close(myVideo2);
% % Function to simulate pause by writing duplicate frames
% function writePauseFrames(video, frame, pauseDuration, frameRate)
%     numFrames = round(pauseDuration * frameRate);
%     for i = 1:numFrames
%         writeVideo(video, frame);
%     end
% end

