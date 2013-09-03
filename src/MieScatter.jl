
module MieScatter

export compute_mie, size_parameter

size_parameter(radius::Real, wavelength::Real) = 2pi * radius / wavelength

function compute_mie(size_param::Real, ref_idx::Number, N_angles::Integer)
    y = size_param * ref_idx
    x_stop = size_param + 4 * size_param^0.3333 + 2.0

    y_modulus = abs(y)
    n_max = int(max(x_stop, y_modulus) + 15)
    delta_angle = pi / 2 / (N_angles-1)
    
    mu_table = [cos((j-1)*delta_angle) for j = 1:N_angles]
    
    d = zeros(Complex128, n_max+1)
    
    # TODO: verify this loop
    for n = n_max-1:-1:1
        rn = n+1
        d[n] = (rn/y) - (1 / (d[n+1] + rn/y))
    end
    
    pi0 = zeros(Float64, N_angles)
    pi1 = ones(Float64, N_angles)
    piX = zeros(Float64, N_angles)
    
    s1 = zeros(Complex128, 2*N_angles-1)
    s2 = zeros(Complex128, 2*N_angles-1)
    tau = zeros(Float64, N_angles)
    
    psi0 = cos(size_param)
    psi1 = sin(size_param)
    chi0 = -sin(size_param)
    chi1 = cos(size_param)
    
    S = zeros(Float64, 2*N_angles-1, 4)
    
    xi1 = Complex128(psi1, -chi1)
    Qsca = 0.0
    
    for n = 1:int(x_stop)
        fn = (2n+1) / (n*(n+1))
        psi = (2n-1) * psi1/size_param - psi0
        chi = (2n-1) * chi1/size_param - chi0
        xi = Complex128(psi, -chi)
        t_a = d[n] / ref_idx + n/size_param
        t_b = d[n] * ref_idx + n/size_param
        an = (t_a * psi - psi1) / (t_a * xi - xi1)
        bn = (t_b * psi - psi1) / (t_b * xi - xi1)
        Qsca += (2n+1) * (abs(an)^2 + abs(bn)^2)
        
        for j = 1:N_angles
            piX[j] = pi1[j]
            jj = 2*N_angles - j
            tau[j] = n * mu_table[j] * piX[j] - (n+1)*pi0[j]
            t = (-1)^n
            p = (-1)^(n-1)
            s1[j] = s1[j] + fn * (an*piX[j] + bn*tau[j])
            s2[j] = s2[j] + fn * (an*tau[j] + bn*piX[j])
            if j != jj
                s1[jj] = s1[jj] + fn * (an*piX[j]*p + bn*tau[j]*t)
                s2[jj] = s2[jj] + fn * (an*tau[j]*t + bn*piX[j]*p)
            end
        end
        
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = Complex128(psi1, -chi1)
        
        for i = 1:N_angles
            pi1[i] = (2n+1)/n * mu_table[i] * pi1[i] - (n+1)*pi0[i] / n
            pi0[i] = piX[i]
        end
    end
    
    Qsca *= 2/size_param^2
    Qext = (4 / size_param^2) * real(s1[1])
    Qback = (4 / size_param^2) * abs(s1[end])^2
    
    for i = 1:2*N_angles-1
        S[i,1] = 0.5 * (abs(s1[i])^2 + abs(s2[i])^2)
        S[i,2] = -0.5 * (abs(s1[i])^2 + abs(s2[i])^2)
        S[i,3] = real(s2[i] * conj(s1[i]))
        S[i,4] = imag(s2[i] * conj(s1[i]))
    end
    
    return S, Qsca, Qext, Qback
    
end


end # module
