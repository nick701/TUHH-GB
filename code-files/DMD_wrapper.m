function DMD_wrapper(filePath, n, r, window, nstep, s, plot_geophones, plot_NRMSE, standard_DMD, aug_DMD, prediction_error, print_NRMSE, print_prediction_error, online_DMD, horizon)
    % Load data from the specified file using the full path
    data = load(filePath);  % Load the specified .mat file (e.g., Seq1, Seq2, or Seq3)

    % Extract the sequence name without ".mat" extension
    [~, name, ~] = fileparts(filePath); 

    % Set variables from loaded data
    G1 = downsample(data.G1V3, n);
    G2 = downsample(data.G2V3, n);
    G3 = downsample(data.G3V3, n);
    G4 = downsample(data.G4V3, n);
    G5 = downsample(data.G5V3, n);
    time = downsample(data.Time, n);
    dt = time(2) - time(1);
    f_all = [G1'; G2'; G3'; G4'; G5']';  % Matrix for all data

    % Assign remaining parameters as specified in the GUI
    num_geo = 5;
    l = ceil(window / 2);  % Number of rows of the augmented matrix
    initstep = 1;
    step = window;
    s_0 = 1;
    s_f = 101;
    s_step = 10;
    t_extra_0 = length(0:dt:(l - 1) * dt);
    folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures';

    % Plotting and control flags from GUI
    plot_geophones = plot_geophones;
    plot_NRMSE = plot_NRMSE;
    standard_DMD = standard_DMD;
    aug_DMD = aug_DMD;
    prediction_error = prediction_error;
    print_NRMSE = print_NRMSE;
    print_prediction_error = print_prediction_error;

    % DMD algorithm settings
    s_analysis = 1;  % Default analysis setting
    extraTime = horizon;  % Prediction horizon

    % Run DMD.m script
    run('DMD.m');
end
