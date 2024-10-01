% Define the function for the backward_forward method for distribution network analysis.
function [v, iteration] = zbus_method(load_data, line_data, slack_bus_voltage, tolerance, max_iter, runs)
    % Arguments include the constant slack bus voltage, convergence value, and max. iterations,
    % all of which can be changed when the function is called.
    
    prompt = {'Ensure the argument data i.e, load and line data are in this format:'
                   'load data: column1 - bus index, column2 - real power (P), and column3 - reactive power (Q)'
                   'line data: column1 - sending bus index, column2 - receiving bus index, column3 - resistance (R), and reactance (X)'
                   '\n'
                   'The default value are slack_bus_voltage = 1.5, tolerance = 1e-6, and max_iter = 100'
                   'If you decide to change any of the constant argument value, such as slack_bus_voltage, tolerance, or max_iter, just pass the new value in their respective position'
                   'If you will jump any of the order, kindly put the constant value of the one you are jumping in its position and input the new value of the argument you are changing'
                   'For example: [v, iteration] = bfs_method(load_data, line_data, 1.5, 1e-4), this is because, I want to keep the slack bus the same but want to change the tolerance'
                   'I did not put the value of max_iter because, it is the last one and I want to keep it the same.'
                   };
    disp(prompt);
    user_input = input('Press Enter to continue or type ''quit'' to exit; ', 's');

     % Check if the user wants to quit
    if strcmp(user_input, 'quit')% strcmp is the MATLAB function for comparing two string  first to last
        disp('Exiting the function as per user request probably because the data does not conform with the requirement');
        iteration = []; % Return an empty array
        v = [];
        return;
    end 
    
    if nargin < 3
        slack_bus_voltage = 1.05; % Default slack bus voltage
    end
    if nargin < 4
        tolerance = 1e-6; % Default tolerance
    end
    
    if nargin < 5
        max_iter = 100; % Maximum iterations
    end
    
    if nargin < 6
        runs = 10; % Maximum iterations
    end
    
    computation_times = [];
    for run = 1:runs
    % System base values
    V_base = 12.66; % Nominal voltage in kV
    S_base = 1; % Base power in MVA
    Z_base = (V_base^2) / S_base; % Base impedance in ohms
    
    % Determine the number of buses from the load data
    num_buses = size(load_data, 1);
    
    % Initialize bus voltages
    v = ones(num_buses, 1) + 1j * zeros(num_buses, 1); % Let all the initial voltage for all the buses be 1.0 p.u.
    v(1) = slack_bus_voltage; % Set the slack bus voltage
    
    % Create load vectors
    P_load = load_data(:, 2) / 1000; % Convert kW to MW and to per-unit
    Q_load = load_data(:, 3) / 1000; % Convert kVar to MVar and to per-unit
    S_load = P_load + 1j * Q_load; % Complex power (S = P + jQ)
    
    % Construct the Y_bus matrix using the helper function script
    Y_bus = build_Y_bus(line_data, num_buses, Z_base);
    
    % Separate Y_bus into components of slack and non-slack
    Y_ss = Y_bus(1, 1);
    Y_sm = Y_bus(1, 2:end);
    Y_ms = Y_bus(2:end, 1);
    Y_mm = Y_bus(2:end, 2:end);
    
    % Initialize arrays to store maximum errors and iteration times
    max_voltage_errors = [];
    cumulative_iter_times = [];
    
    % Measure the computational time
    start_time = tic;
    
    % Z_bus (Fixed-Point iteration) method
    
    for iteration = 1:max_iter % Iterate up to max_iter times to perform the BFS algorithm.
        % Start timing for the iteration
        iter_start_time = tic;
        
        % Calculate current injection for all buses using I_i = sum(Y_ij * V_j)
        I = Y_bus * v;
        
        % Separate currents into slack and non-slack parts
        I_s = I(1); % Slack current
        I_m = I(2:end); % Non-slack currents
        v_m_old = v(2:end); % Old voltages for non-slack buses
       
        % Update non-slack bus voltages
       for i = 1:length(v_m_old)
        % Calculate the summation part of the formula: sum(Y_ij * V_j) where j ≠ i
        sum_term = 0;
        for j = 1:length(v)
            if j ~= (i + 1) % Adjust to skip the diagonal term
                sum_term = sum_term + Y_bus(i + 1, j) * v(j);
            end
        end
        
        % Update the voltage for the i-th non-slack bus
        v_m_new(i) = (1 / Y_bus(i + 1, i + 1)) * (I_m(i) - sum_term);
        end
        
        % Update the voltages
        v(2:end) = v_m_new.';

        % End timing the iteration
        iter_time = toc(iter_start_time);
        
        % Check for convergence
        max_diff = max(abs(v_m_new - v_m_old));
        
        % Store max voltage difference and cumulative iteration time
        max_voltage_errors = [max_voltage_errors, max_diff];
        if iteration == 1
            cumulative_iter_times = [cumulative_iter_times, iter_time];
        else
            cumulative_iter_times = [cumulative_iter_times, cumulative_iter_times(end) + iter_time];
        end
             
        % Print max voltage difference at an iteration and its respective
        % computation time till that iteration.
        fprintf('Iteration %d: max voltage difference = %.10f\n', iteration, max_diff);
        fprintf('Time taken till %d iteration: %.4f seconds\n', iteration, iter_time);
        
        if max_diff <= tolerance % Check if the maximum voltage difference between iterations is less than the tolerance.
            break; % Break the loop if convergence is achieved.
        end
    end
    
    % End computational time measurement
    total_computation_time = toc(start_time);
    computation_times = [computation_times, total_computation_time];
    end    

    % Prepare results and print them
    fprintf('Converged in %d iterations.\n', iteration); % Print the number of iterations it took to converge.
    fprintf('Bus Voltages:\n'); % Print a header for the bus voltages.
    for i = 1:num_buses % Iterate over each bus voltage.
        magnitude = abs(v(i)); % Calculate the magnitude of the voltage.
        phase_angle = rad2deg(angle(v(i))); % Calculate the angle of the voltage in degrees.
        rectangular_form = sprintf('%.4f + %.4fj', real(v(i)), imag(v(i)));
        polar_form = sprintf('%.4f ∠ %.2f°', magnitude, phase_angle);
        fprintf('Bus %d: Rectangular form: %s p.u., Polar form: %s\n', i, rectangular_form, polar_form); % Print the voltage in both forms.
    end
    
    fprintf('Total computational time: %.4f seconds.\n', total_computation_time);
   
    % Compute system loss, substation power, and other results
    
    % Calculate system loss
    total_active_loss = 0;
    total_reactive_loss = 0;
    for i = 1:size(line_data, 1)
        sending_bus = line_data(i,1);
        receiving_bus = line_data(i,2);
        resistance = line_data(i, 3); % Get the resistance (R)
        reactance = line_data(i, 4) * 1j; % Get the reactance (X) and multiply by 1j
        impedance = (resistance + reactance)/Z_base;
        
        % Calculate the current flowing in the line
        I_line = (v(sending_bus) - v(receiving_bus)) / impedance;
        
        % Calculate power losses
        P_loss = abs(I_line)^2 * resistance; % Active power loss
        Q_loss = abs(I_line)^2 * reactance; % Reactive power loss
        
        % Accumulate total power losses
        total_active_loss = total_active_loss + P_loss;
        total_reactive_loss = total_reactive_loss + Q_loss;
        fprintf('Loss in branch %d-%d: %.4f p.u. + %.4fj p.u.\n', sending_bus, receiving_bus, P_loss, Q_loss);
    end

    fprintf('Total active power loss: %.4f p.u. = %.4f kW\n', total_active_loss, total_active_loss * 1000);
    fprintf('Total reactive power loss: %.4f p.u. = %.4f kVAR\n', total_reactive_loss, total_reactive_loss * 1000);
    
    % Calculate substation power
    substation_active_power = sum(P_load) + total_active_loss;
    substation_reactive_power = sum(Q_load) + total_reactive_loss;
    fprintf('Substation active power: %.4f p.u. = %.4f MW\n', substation_active_power, substation_active_power);
    fprintf('Substation reactive power: %.4f p.u. = %.4f MVAR\n', substation_reactive_power, substation_reactive_power);
    
    % Assume Y_bus and v are defined

    % Debug output
    disp('Y_bus Matrix:');
    disp(Y_bus); % Print Y_bus matrix

    disp('Initial Voltages:');
    disp(v); % Print bus voltages

    disp('Current Injections I:');
    disp(I); % Print current injections

    
    % Find minimum and maximum voltages and their corresponding bus indices
    [min_voltage, min_index] = min(abs(v));
    [max_voltage, max_index] = max(abs(v));

    % Print the minimum and maximum voltages along with the bus indices
    fprintf('Minimum voltage = %.4f pu at bus %d\n', min_voltage, min_index);
    fprintf('Maximum voltage = %.4f pu at bus %d\n', max_voltage, max_index);

    % Store the formatted computation times as strings so as to return it
    % as a formated list of string in 4dp.
    formatted_times = arrayfun(@(x) sprintf('%.4f', x), computation_times, 'UniformOutput', false);
    
    % Print the formatted list of computation times
    fprintf('Computation times for each run: [%s]\n', strjoin(formatted_times, ', '));
    
    % Calculate and print the average computation time
    average_time = mean(computation_times);
    fprintf('Average computation time for the model after %d runs: %.4f seconds\n', runs, average_time);
    
    % Extract the magnitudes of the voltages
    voltage_magnitudes = abs(v);
    voltage_angles = angle(v);
      
    % Create the graph using the provided line data to visualize the system
    % network
    G = digraph(line_data(:, 1), line_data(:, 2));
    % Plot the graph
    figure;
    h = plot(G);
    % Customize the appearance
    h.NodeColor = 'r';  % Node color
    h.EdgeColor = 'b';  % Edge color
    h.ArrowSize = 10;   % Arrow size
    
    % Create a plot of bus voltages
    figure;
    subplot(1, 2, 1); % 1 row, 2 columns, first subplot
    plot(1:length(voltage_magnitudes), voltage_magnitudes, '-');
    xlabel('Bus Number');
    ylabel('Voltage Magnitude (p.u.)');
    title('Bus Voltage Profile');
    % grid on;
    
    % Create a plot of bus voltages angle
    % figure;
    subplot(1, 2, 2); % 1 row, 2 columns, second subplot
    plot(1:length(voltage_angles), voltage_angles, '-');
    xlabel('Bus Number');
    ylabel('Voltage Angle (radian)');
    title('Bus Voltage Angle');
    % grid on;
    
    % Plot the maximum voltage error vs. computation time till iteration
    figure;
    plot(cumulative_iter_times, max_voltage_errors, '-o');
    xlabel('Computation Time (seconds)');
    ylabel('Maximum Error');
    title('Maximum Error vs. Computation Time');
    % grid on;
    
    % Write results to a text file with UTF-8 encoding
    fid_txt = fopen(sprintf('bus_voltage_%d_bus.txt', num_buses), 'w');
    fprintf(fid_txt, 'Converged in %d iterations.\n', iteration);
    fprintf(fid_txt, 'Bus Voltages:\n');
    fprintf(fid_txt, 'Bus Number\tRectangular Form\tPolar Form\n');
    for i = 1:num_buses
        fprintf(fid_txt, '%d\t%.4f + %.4fj\t%.4f ∠ %.2f°\n', i, real(v(i)), imag(v(i)), magnitude, phase_angle);
    end
    fclose(fid_txt);
    
    % Write results to a CSV file with UTF-8 encoding
    fid_csv = fopen(sprintf('bus_voltage_%d_bus.csv', num_buses), 'w');
    fprintf(fid_csv, 'Bus Number,Rectangular Form,Polar Form\n');
    for i = 1:num_buses
        fprintf(fid_csv, '%d,%.4f + %.4fj,%.4f ∠ %.2f°\n', i, real(v(i)), imag(v(i)), magnitude, phase_angle);
    end
    fclose(fid_csv);

end

    