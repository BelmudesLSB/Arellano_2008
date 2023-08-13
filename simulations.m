function [Consumption, Bonds, Income] = simulations(V_star, Vr_star, Bp_star, Dp_star, Q_star, T_sim_star, P_star, Y_star, B_star, ny_star, nb_star, Burn_star)
    %%% Perform calculations here using the equilibrium objects:
    % V_star: Value function at equilibra.
    % Vr_star: Value function of repayment at equilibra.
    % Bp_star: Policy function for bonds function at equilibra.
    % Dp_star: Policy function for default at equilibra.
    % Q_star: Equilibrium price of debt.
    % T_sim_star: Number of periods to simulate.
    % P: Transition probabilities constructed following Tauchen.
    % Y_star: Income Grid.
    % B_star: Bond grid.
    % Burn_star: Observations at the start of the simulation deleted

    %Set seed for random sampling:
    rng(123);

    % Simulations are vectors of size (1xT_sim).

    % Initialize the vectors:
    Consumption = zeros(1,T_sim_star);
    Bonds = zeros(1,T_sim_star);
    Income_idx =zeros(1,T_sim_star);
    Bonds_idx = zeros(1,T_sim_star);
    Default_p = zeros(1,T_sim_star);
    Y_effective = zeros(1,T_sim_star);

    % Get the simulation for income: Let the simulation start in the worst
    % possible outcome.

    % Starting values for the simulation.
    Income_idx(1) = 1;
    Bonds_idx(1) = 1;

    for t = 2:T_sim_star
        % Check the default decision today, given output and bonds:
        Default_p(t-1) = Dp_star(Income_idx(t-1), Bonds_idx(t-1));
        % Draw a random number to determine the next state:
        next_state = randsample(ny_star, 1, true, P_star(Income_idx(t-1), :));
        % Update the current state to the next state
        Income_idx(t) = next_state;
        Bonds_idx(t) = Bp_star(Income_idx(t-1), Bonds_idx(t-1));
    end

    % Burn Some observations:

    Income = Y_star(Income_idx(Burn_star+1:end));
    %display(Default_p);
    
   
end