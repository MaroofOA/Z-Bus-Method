function [v, iteration] = zbus_debug(load_data, line_data, slack_bus_voltage, tolerance, max_iter)

    % Set default values if arguments are not provided
    if nargin < 3, slack_bus_voltage = 1.05; end
    if nargin < 4, tolerance = 1e-6; end
    if nargin < 5, max_iter = 100; end
    
    % System base values
    V_base = 12.66; % Nominal voltage in kV
    S_base = 1; % Base power in MVA
    Z_base = (V_base^2) / S_base; % Base impedance in ohms

    num_buses = size(load_data, 1); % Number of buses
    P_load = load_data(:, 2) / 1000; % Real power in p.u.
    Q_load = load_data(:, 3) / 1000; % Reactive power in p.u.
      
      
    % Initialize bus voltages
    v = ones(num_buses, 1) + 1j * zeros(num_buses, 1); % All initial voltages = 1.0 p.u.
    v(1) = slack_bus_voltage; % Set the slack bus voltage

    % Construct the Y_bus matrix
    Y_bus = build_Y_bus(line_data, num_buses, Z_base);

    % Split Y_bus into components for slack and non-slack buses
    Y_ms = Y_bus(2:end, 1); % Interaction between non-slack buses and slack bus
    Y_mm = Y_bus(2:end, 2:end); % Interaction between non-slack buses

    % Z-bus method loop
    for iteration = 1:max_iter

        % Correct current injection formula considering negative loads
        I = -(P_load - 1j * Q_load) ./ conj(v); % Power/current injection

        % Separate currents for non-slack buses
        I_m = I(2:end); % Non-slack currents
        V_old = v(2:end); % Old voltages for non-slack buses
            
        % Update non-slack bus voltages
        v(2:end) = Y_mm \ (I_m - Y_ms * v(1));

        % Check for convergence
        max_diff = max(abs(v(2:end) - V_old));

        if max_diff <= tolerance % Convergence check
            break;
        end
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
  
    % Compute system loss, substation power, and other results
    % Calculate system loss
    [total_active_loss, total_reactive_loss] = calculate_system_loss(num_buses, line_data, v, Z_base);
    
    fprintf('Converged in %d iterations.\n', iteration); % Print the number of iterations it took to converge.
    
    fprintf('Total active power loss: %.4f p.u. = %.4f kW\n', total_active_loss, total_active_loss * 1000);
    fprintf('Total reactive power loss: %.4f p.u. = %.4f kVAR\n', total_reactive_loss, total_reactive_loss * 1000);
    
    % Calculate substation power
    substation_active_power = sum(P_load) + total_active_loss;
    substation_reactive_power = sum(Q_load) + total_reactive_loss;
    fprintf('Substation active power: %.4f p.u. = %.4f MW\n', substation_active_power, substation_active_power);
    fprintf('Substation reactive power: %.4f p.u. = %.4f MVAR\n', substation_reactive_power, substation_reactive_power);
    
    % Find minimum and maximum voltages and their corresponding bus indices
    [min_voltage, min_index] = min(abs(v));
    [max_voltage, max_index] = max(abs(v));

    % Print the minimum and maximum voltages along with the bus indices
    fprintf('Minimum voltage = %.4f pu at bus %d\n', min_voltage, min_index);
    fprintf('Maximum voltage = %.4f pu at bus %d\n', max_voltage, max_index);

    % Extract the magnitudes of the voltages
    voltage_magnitudes = abs(v);
    voltage_angles = angle(v);
          
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
    
end

    