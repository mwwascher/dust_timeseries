using CSV
using Distributions
using DataFrames
using Plots
using Random
using StatsPlots


function f_lam(t, onset, p, h, r)
    if t < onset
        return 0
    elseif onset < t <= p
        return exp(log(h + 1) / (p - onset) * (t - onset)) - 1
    elseif p < t <= r
        return exp(log(h + 1) / (r - p) * (r - t)) - 1
    else
        return 0
    end
end

function lam_int(a, b, lam_t)
    net = a:0.01:b  # Create a range with step size 0.01
    out = 0.0       # Initialize output
    for i in net
        out += 0.01 * f_lam(i, lam_t[:onset], lam_t[:peaktime], lam_t[:peakh], lam_t[:recovery])
    end
    return out
end

function lam_arr(a, b, lam_t)
    net = a:0.01:b  # Create a range with step size 0.01
    f = Vector{Union{Float64, Missing}}(undef, length(net))  # Initialize a vector
    
    # Populate f using a for loop
    for (counter, i) in enumerate(net)
        f[counter] = f_lam(i, lam_t[:onset], lam_t[:peaktime], lam_t[:peakh], lam_t[:recovery])
    end
    
    # Normalize f to create probabilities for multinomial sampling
    probabilities = f ./ sum(f)
    
    # Multinomial sampling
    idx = sample(1:length(net), ProbabilityVector(probabilities))  # `idx` is the sampled index
    
    return net[idx]  # Return the corresponding value from `net`
end

function p_int(a, b, lam_t, mu, xi)
    net = a:0.01:b  # Create a range with step size 0.01
    out = 0.0       # Initialize the output
    
    for i in net
        out += 0.01 * exp(xi * (i - b)) * f_lam(i, lam_t[:onset], lam_t[:peaktime], lam_t[:peakh], lam_t[:recovery])
    end
    
    return out / mu  # Return the normalized result
end

function MCMC_timecal()

    cd("C:\\Users\\matth\\OneDrive\\Documents\\JuliaWD\\")

    # Parameters
    xi = log(10) / 7
    mu_p = 4.2
    sigma_p = 0.5
    mu_r = 7.3
    sigma_r = 0.6

    nreps = 15000
    #Library
    #T = [7,13,14,29]
    #n = [33,29,8,37]
    #dust_cal = [505084,348497,752,226785]
    #amt = [fill(156,n[1]);fill(201,n[2]);fill(56,n[3]);fill(202,n[4])]

    #Union
    #T = [6,21,28,43]
    #n = [28,19,15,31]
    #dust_cal = [15920,10350,2415,625]
    #amt = [fill(51,n[1]);fill(52,n[2]);fill(53,n[3]);fill(48,n[4])]

    #RPAC
    #T = [7,13,20,28]
    #n = [35,23,16,18]
    #dust_call = [3224,623,1318,484]
    #amt = [50,58,55,51]

    #RPAC
    T = [7]
    n = [35]
    dust_cal = [3224]
    amt = [fill(50,n[1])]
    dust_old = 0

    period = length(T)
    dust_cal_adj = deepcopy(dust_cal)
    for j in 1:period
        if j == 1
            dust_cal_adj[j] = dust_cal[j] - min(convert(Int64, round(.9*dust_old)),rand(Binomial(dust_old, exp(-xi*(T[j])))))
        else
            dust_cal_adj[j] = dust_cal[j] - min(convert(Int64, round(.9*dust_cal[j])),rand(Binomial(dust_cal[j-1], exp(-xi*(T[j] - T[j-1])))))
        end
    end


    # Initialize data structures
    sum_n = sum(n)
    lam_t = DataFrame(onset = Vector{Union{Float64, Missing}}(undef, sum_n),
                    peaktime = Vector{Union{Float64, Missing}}(undef, sum_n),
                    peakh = Vector{Union{Float64, Missing}}(undef, sum_n),
                    recovery = Vector{Union{Float64, Missing}}(undef, sum_n))

    lam_t_new = deepcopy(lam_t)
    mu_all = Matrix{Union{Float64, Missing}}(undef, sum_n, period)
    p_all = Matrix{Union{Float64, Missing}}(undef, sum_n, period)
    mu_all_new = deepcopy(mu_all)
    p_all_new = deepcopy(p_all)

    beta = 0.14
    beta_samp = Vector{Float64}(undef, nreps)
    loglik = -1e15  # Initial log-likelihood

    # Simulation loop
    for l in 1:nreps
        beta_new = 0.0
        while beta_new <= 0
            beta_new = beta + rand(Normal(0, 0.05))
        end

        for i in 1:period
            counter = 1
            for j in 1:sum_n
                p = rand(Normal(mu_p, sigma_p))
                r = rand(Normal(mu_r, sigma_r))
                h = rand(Exponential(1/beta))
                u = rand(Uniform(0, p + r))
                onset = T[i] - u
                
                lam_t.onset[counter] = onset
                lam_t.peaktime[counter] = onset + p
                lam_t.peakh[counter] = h
                lam_t.recovery[counter] = onset + p + r
                if i == 1
                    mu_all[counter, i] = lam_int(0,T[i], lam_t[counter, :])
                    p_all[counter, i] = p_int(0, T[i], lam_t[counter, :], mu_all[counter, i], xi)
                else
                    mu_all[counter, i] = lam_int(T[i-1],T[i], lam_t[counter, :])
                    p_all[counter, i] = p_int(T[i-1], T[i], lam_t[counter, :], mu_all[counter, i], xi)
                end

                p_new = rand(Normal(mu_p, sigma_p))
                r_new = rand(Normal(mu_r, sigma_r))
                h_new = rand(Exponential(1/beta_new))
                u_new = rand(Uniform(0, p_new + r_new))
                onset_new = T[i] - u_new

                lam_t_new.onset[counter] = onset_new
                lam_t_new.peaktime[counter] = onset_new + p_new
                lam_t_new.peakh[counter] = h_new
                lam_t_new.recovery[counter] = onset_new + p_new + r_new
                if i == 1
                    mu_all_new[counter, i] = lam_int(0, T[i], lam_t_new[counter, :])
                    p_all_new[counter, i] = p_int(0, T[i], lam_t_new[counter, :], mu_all_new[counter, i,], xi)
                else
                    mu_all_new[counter, i] = lam_int(T[i-1], T[i], lam_t_new[counter, :])
                    p_all_new[counter, i] = p_int(T[i-1], T[i], lam_t_new[counter, :], mu_all_new[counter, i,], xi)
                end

                counter += 1
            end
        end

        grand_mean = sum(skipmissing(amt .* mu_all .* p_all))
        grand_mean_new = sum(skipmissing(amt .* mu_all_new .* p_all_new))
        
        dust_cal_adj_new = deepcopy(dust_cal)
        for j in 1:period
            if j == 1
                dust_cal_adj_new[j] = dust_cal[j]
            else
                dust_cal_adj_new[j] = dust_cal[j] - min(convert(Int64, round(.9*dust_cal[j])),rand(Binomial(dust_cal_adj_new[j-1], exp(-xi*(T[j] - T[j-1])))))
            end
        end

        loglik = logpdf(Poisson(grand_mean[1]), sum(dust_cal_adj))
        loglik_new = logpdf(Poisson(grand_mean_new[1]), sum(dust_cal_adj_new))

        ar = min(1, exp(loglik_new - loglik))
        flip = rand(Bernoulli(ar))

        if flip == 1
            beta = beta_new
        end

        beta_samp[l] = beta
    end

    # Visualization
    histogram(beta_samp, bins=30, xlabel="Beta", ylabel="Frequency", title="Beta Samples")
    plot(beta_samp, xlabel="Iteration", ylabel="Beta", title="Trace Plot", lw=2)
    df = DataFrame(Value = beta_samp)
    CSV.write("beta_samp1.csv", df)
end
