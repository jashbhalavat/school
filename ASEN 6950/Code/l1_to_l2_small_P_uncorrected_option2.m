%% Plot uncorrected trajectory - option 2, single cross

% Found these nubmers by manually looking through plot:
% find(abs(saved_final_state_init_single_cross(:,1)-(-0.0811833)) < 1e-5)
% find(abs(saved_final_state_final_single_cross(:,1)-(-0.0869406)) < 1e-5)
% find(abs(saved_final_state_final_single_cross(:,1)-(-0.0810933)) < 1e-5)
uncorrected_init_idx = 2876;
uncorrected_final_idx = 1904;
uncorrected_init_i = saved_final_state_init_single_cross(uncorrected_init_idx,4);
uncorrected_init_j = saved_final_state_init_single_cross(uncorrected_init_idx,5);
% uncorrected_init_seconds = saved_final_state_init_mult_cross(uncorrected_init_idx,6);
uncorrected_final_i = saved_final_state_final_single_cross(uncorrected_final_idx,4);
uncorrected_final_j = saved_final_state_final_single_cross(uncorrected_final_idx,5);
% uncorrected_final_seconds = saved_final_state_final_mult_cross(uncorrected_final_idx,6);
uncorrected_init_state0 = xout_lyapunov_init(uncorrected_init_j,:);
uncorrected_final_state0 = xout_lyapunov_final(uncorrected_final_j,:);
uncorrected_init_thrust = thrust_direction(:,uncorrected_init_i);
uncorrected_final_thrust = thrust_direction(:,uncorrected_final_i);

figure(8)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
plot(xout_lyapunov_init(1:uncorrected_init_j,1), xout_lyapunov_init(1:uncorrected_init_j,2), 'blue', 'LineWidth', 2)
plot(xout_lyapunov_final(uncorrected_final_j:end,1), xout_lyapunov_final(uncorrected_final_j:end,2), 'red', 'LineWidth', 2)
grid on
title("Uncorrected Trajectory Option 2")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)


fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
[uncorrected_init_tout, uncorrected_init_xout] = ode113(fun, [0, 45], [uncorrected_init_state0, 1], options_single_cross);
x_1_f = xout(end,:)';
plot(uncorrected_init_xout(:,1), uncorrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
[uncorrected_final_tout, uncorrected_final_xout] = ode113(fun, [0, -25], [uncorrected_final_state0, 0.9], options_single_cross);
plot(uncorrected_final_xout(:,1), uncorrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
hold off
grid on
legend("L2", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")

