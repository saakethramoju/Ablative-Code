
L = 1.0;  % length (m)
T_total = 1.0;  % time (s)
N = 10;  % Number of nodes
M = 100;  % Number of time steps

dx = L / (N - 1);
dt = T_total / (M - 1);
alpha = 0.01;  

T = zeros(M, N);
rho = ones(M, N);  

T(1, :) = sin(pi * linspace(0, L, N));

T(:, 1) = 0; %replace with surface energy balance**
T(:, end) = 273; % K

q_diff = 1.0;
q_rad = 1.0;
m_dot_c = 1.0;
m_dot_g = 1.0;
q_star = 1.0;
rho_w = 1.0;
h_w = 1.0;
h_c = 1.0;
h_g = 1.0;
q_cond = 1.0;
E_i = 1.0;
R = 1.0;
k_i = 1.0;
delta_theta = dt;

r = alpha * dt / (2 * dx^2);

for m = 2:M
    for i = 1:N
        if M ~= 1
            rho(m, i) = rho(m-1, i) + density_change_rate_mi_not_1(rho(m-1, i), rho_ri, T(m-1, i), M);
        else
            rho(m, i) = rho(m-1, i) + density_change_rate_mi_1(rho(m-1, i), rho_ri, T(m-1, i));
        end
    end
    
    T(m, :) = crank_nicolson_step(T(m-1, :), rho(m, :), N, r, q_diff, q_rad, m_dot_c, m_dot_g, q_star, q_cond, rho_w, h_w, h_c, h_g, dx, dt);
end

%plotting stuff eventually

function rate = density_change_rate_mi_not_1(rho_i, rho_ri, T_i, m_i) %wrong way of doing it, need to calculate explicit, ask Dec how to incorporate, also ask why that one paper has it different based on if we are at first node and no other paper has that
    term1 = (rho_i - rho_ri) / delta_theta;
    term2 = (rho_ri - rho_i) / delta_theta + (rho_i - rho_ri) * (1 - m_i) / m_i * k_i * exp(-E_i / (R * T_i));
    rate = term1 + term2;
end

function rate = density_change_rate_mi_1(rho_i, rho_ri, T_i)
    rate = (rho_i - rho_ri) / delta_theta * (exp(-E_i / (R * T_i) * delta_theta) - 1);
end


function T_new = crank_nicolson_step(T_prev, rho_prev, N, r, q_diff, q_rad, m_dot_c, m_dot_g, q_star, q_cond, rho_w, h_w, h_c, h_g, dx, dt)
    A = zeros(N, N);
    B = zeros(N, 1);

    for n = 2:N-1
        k_n_1 = k_i;  % n-1
        k_n = k_i;    % n
        k_n1 = k_i;   % n+1
        
        A_n_1 = 1.0;  %node n-1
        A_n = 1.0;    %node n
        A_n1 = 1.0;   %node n+1

        delta_x_1 = dx; %n-1 and n
        delta_x = dx;   %n and n+1
        
        h_n = h_g;   %node n
        h_n1 = h_c;  %n+1
        
        term1 = (T_prev(n-1) - T_prev(n)) / delta_x_1 / (k_n * A_n_1);
        term2 = (T_prev(n) - T_prev(n+1)) / delta_x / (k_n * A_n);
        term3 = h_n * (T_prev(n) - T_prev(n));
        term4 = h_n1 * (T_prev(n+1) - T_prev(n+1));
        term5 = (rho_prev(n+1) - rho_prev(n)) / delta_theta;
        term6 = (rho_prev(n) - rho_prev(n-1)) / delta_theta;
        
        A(n, n-1) = -r;
        A(n, n) = 1 + 2 * r;
        A(n, n+1) = -r;

        B(n) = r * T_prev(n-1) + (1 - 2 * r) * T_prev(n) + r * T_prev(n+1) + term3 + term4 + term5 + term6; %wrong fdm formula fix this
    end
    
    A(1, 1) = 1 + r;
    A(1, 2) = -r;
    B(1) = r * T_prev(2) + (1 - r) * T_prev(1) ...
        + (q_diff + q_rad + m_dot_c * h_c + m_dot_g * h_g - q_star - q_cond) * dt / (rho_w * h_w);
    A(N, N-1) = -r;
    A(N, N) = 1 + r;
    B(N) = r * T_prev(N-1) + (1 - r) * T_prev(N);
    

    T_new = thomas_algorithm(A, B, N);
end

function x = thomas_algorithm(A, B, N) %test with simple heat conduction to see if works
    c_prime = zeros(N-1, 1);
    d_prime = zeros(N, 1);

    c_prime(1) = A(1, 2) / A(1, 1);
    d_prime(1) = B(1) / A(1, 1);

    for i = 2:N-1
        denom = A(i, i) - A(i, i-1) * c_prime(i-1);
        c_prime(i) = A(i, i+1) / denom;
        d_prime(i) = (B(i) - A(i, i-1) * d_prime(i-1)) / denom;
    end

    d_prime(N) = (B(N) - A(N, N-1) * d_prime(N-1)) / (A(N, N) - A(N, N-1) * c_prime(N-1));

    x = zeros(N, 1);
    x(N) = d_prime(N);

    for i = N-1:-1:1
        x(i) = d_prime(i) - c_prime(i) * x(i+1); 
    end
end




