diff_s=s-test_x;
cumulative_diff = cumsum(abs(diff_s)) * dt;
plot(cumulative_diff)
plot(diff_s)
total_energy = sum(abs(s).^2) * dt;
error_energy = sum(abs(diff_s).^2) * dt;
error_energy_ratio = error_energy / total_energy * 100;
