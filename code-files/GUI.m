function DMD_GUI
    % Create a figure for the GUI
    fig = uifigure('Name', 'DMD Parameter Testing', 'Position', [100, 100, 400, 600]);
    
    % Define the data directory path
    dataDir = '/Users/NickT/Documents/TUHH-MATLAB-GitHub-Repo/TUHH-GB-Repo/mat-files/';
    
    % Dropdown for selecting data sequence
    uilabel(fig, 'Position', [10, 560, 200, 22], 'Text', 'Data Sequence Name');
    seqDropdown = uidropdown(fig, 'Position', [220, 560, 100, 22], 'Items', {'Seq1', 'Seq2', 'Seq3'});
    
    % Text fields for numerical inputs
    uilabel(fig, 'Position', [10, 520, 200, 22], 'Text', 'Downsampling Factor (n)');
    nEdit = uieditfield(fig, 'numeric', 'Position', [220, 520, 100, 22], 'Value', 100);

    uilabel(fig, 'Position', [10, 480, 200, 22], 'Text', 'Rank Truncation (r)');
    rEdit = uieditfield(fig, 'numeric', 'Position', [220, 480, 100, 22], 'Value', 5);

    uilabel(fig, 'Position', [10, 440, 200, 22], 'Text', 'DMD Window Size');
    windowEdit = uieditfield(fig, 'numeric', 'Position', [220, 440, 100, 22], 'Value', 500);

    uilabel(fig, 'Position', [10, 400, 200, 22], 'Text', 'Window Step Size (nstep)');
    nstepEdit = uieditfield(fig, 'numeric', 'Position', [220, 400, 100, 22], 'Value', 1);

    uilabel(fig, 'Position', [10, 360, 200, 22], 'Text', 'Delay Coordinates (s)');
    sEdit = uieditfield(fig, 'numeric', 'Position', [220, 360, 100, 22], 'Value', 50);

    % Checkboxes for toggling flags
    plotGeophonesCheck = uicheckbox(fig, 'Position', [10, 320, 150, 22], 'Text', 'Plot Geophones');
    plotNRMSECheck = uicheckbox(fig, 'Position', [10, 290, 150, 22], 'Text', 'Plot NRMSE');
    standardDMDCheck = uicheckbox(fig, 'Position', [10, 260, 150, 22], 'Text', 'Standard DMD');
    augDMDCheck = uicheckbox(fig, 'Position', [10, 230, 150, 22], 'Text', 'Augmented DMD');
    predictionErrorCheck = uicheckbox(fig, 'Position', [10, 200, 150, 22], 'Text', 'Prediction Error');
    printNRMSECheck = uicheckbox(fig, 'Position', [10, 170, 150, 22], 'Text', 'Print NRMSE');
    printPredictionErrorCheck = uicheckbox(fig, 'Position', [10, 140, 150, 22], 'Text', 'Print Prediction Error');
    onlineDMDCheck = uicheckbox(fig, 'Position', [10, 110, 150, 22], 'Text', 'Online DMD');

    % Text field for prediction horizon
    uilabel(fig, 'Position', [10, 80, 200, 22], 'Text', 'Prediction Horizon');
    horizonEdit = uieditfield(fig, 'numeric', 'Position', [220, 80, 100, 22], 'Value', 0);

    % Run button
    runButton = uibutton(fig, 'push', 'Position', [150, 40, 100, 30], 'Text', 'Run Analysis', 'ButtonPushedFcn', @(runButton, event) runAnalysis);

    % Callback function for Run button
    function runAnalysis
        % Collect input values from the GUI
        name = seqDropdown.Value;
        filePath = fullfile(dataDir, [name, '.mat']); % Build full path to the .mat file
        n = nEdit.Value;
        r = rEdit.Value;
        window = windowEdit.Value;
        nstep = nstepEdit.Value;
        s = sEdit.Value;
        plot_geophones = plotGeophonesCheck.Value;
        plot_NRMSE = plotNRMSECheck.Value;
        standard_DMD = standardDMDCheck.Value;
        aug_DMD = augDMDCheck.Value;
        prediction_error = predictionErrorCheck.Value;
        print_NRMSE = printNRMSECheck.Value;
        print_prediction_error = printPredictionErrorCheck.Value;
        online_DMD = onlineDMDCheck.Value;
        horizon = horizonEdit.Value;

        % Run the DMD script with these parameters, including the file path
        % Assuming your DMD.m uses function structure with these inputs:
        DMD(filePath, n, r, window, nstep, s, plot_geophones, plot_NRMSE, standard_DMD, aug_DMD, prediction_error, print_NRMSE, print_prediction_error, online_DMD, horizon);
    end
end
