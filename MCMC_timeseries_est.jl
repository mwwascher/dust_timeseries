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

function MCMC_timeest()

    cd("C:\\Users\\matth\\OneDrive\\Documents\\JuliaWD\\")

    # Parameters
    xi = log(10) / 7
    mu_p = 4.2
    sigma_p = 0.5
    mu_r = 7.3
    sigma_r = 0.6

    

    df = CSV.read("beta_post.csv", DataFrame)
    beta_in = df[:, :x]

    nreps = 15000

    A = DataFrame(DataFrame(zeros(nreps, 1), :auto))

    #Library
    #T = [7,14,28,39]
    #n = [10,10,10,10]
    #dust_cal = [15982,12476,65556,30213]
    #amt = [fill(200,n[1]);fill(154,n[2]);fill(197,n[3]);fill(199,n[4])]

    #Union
    #T = 
    #n = [10,10,10,10]
    #dust_cal = [11734,2532,930,1511]
    #amt = [fill(49,n[1]);fill(51,n[2]);fill(47,n[3]);fill(50,n[4])]

    #RPAC
    #T = [6,21,33]
    #n = [10,10,10]
    #dust_cal = [4736,1340,757]
    #amt = [fill(48,n[1]);fill(53,n[2]);fill(52,n[3])]



    T = [7]
    n = [10]
    dust_cal = [15982]
    dust_old = [17598]
    amt = [fill(200,n[1])]

    period = length(T)
    dust_cal_adj = deepcopy(dust_cal)
    for j in 1:period
        if j == 1
            dust_cal_adj[j] = dust_cal[j] - min(convert(Int64, round(.9*dust_cal[j])),rand(Binomial(dust_old[1], exp(-xi*(T[j])))))
        else
            dust_cal_adj[j] = dust_cal[j] - min(convert(Int64, round(.9*dust_cal[j])),rand(Binomial(dust_cal[j-1], exp(-xi*(T[j] - T[j-1])))))
        end
    end


    # Initialize data structures
    #sum_n = sum(n)
    #lam_t = DataFrame(onset = Vector{Union{Float64, Missing}}(undef, sum_n),
    #                peaktime = Vector{Union{Float64, Missing}}(undef, sum_n),
    #                peakh = Vector{Union{Float64, Missing}}(undef, sum_n),
    #                recovery = Vector{Union{Float64, Missing}}(undef, sum_n))

    #lam_t_new = deepcopy(lam_t)
    #mu_all = Matrix{Union{Float64, Missing}}(undef, sum_n, period)
    #p_all = Matrix{Union{Float64, Missing}}(undef, sum_n, period)
    #mu_all_new = deepcopy(mu_all)
    #p_all_new = deepcopy(p_all)

    #beta = mean(beta_in)
    #beta_samp = Vector{Float64}(undef, nreps)
    #n_samp = DataFrame([fill(0, nreps) for _ in 1:4], [:n1, :n2, :n3, :n4])
    n_samp = DataFrame([fill(0, nreps) for _ in 1:1], [:n1])
    loglik = -1e15  # Initial log-likelihood
    counter_n = 1

    # Simulation loop
    for l in 1:nreps
        
        beta = rand(beta_in)
        n_new = n + rand([-5,-4,-3,-2,-1,0,1,2,3,4,5],1)#rand([-2,-1,0,1,2],4)
        while(maximum(n_new) > 50 || minimum(n_new) < 1)
            n_new = n + rand([-5,-4,-3,-2,-1,0,1,2,3,4,5],1)#rand([-2,-1,0,1,2],4)
        end
        # if counter_n > 4
        #     counter_n = 1
        # end
        # n_new = n
        # n_new[counter_n] = rand(1:50)
        # counter_n = counter_n + 1

        # Initialize data structures for current and proposed n
        sum_n = sum(n)
        #amt = [fill(49,n[1]);fill(51,n[2]);fill(47,n[3]);fill(50,n[4])]
        amt = [fill(49,n[1])]
        lam_t = DataFrame(onset = Vector{Union{Float64, Missing}}(undef, sum_n),
                    peaktime = Vector{Union{Float64, Missing}}(undef, sum_n),
                    peakh = Vector{Union{Float64, Missing}}(undef, sum_n),
                    recovery = Vector{Union{Float64, Missing}}(undef, sum_n))

        #lam_t_new = deepcopy(lam_t)
        mu_all = Matrix{Union{Float64, Missing}}(undef, sum_n, period)
        p_all = Matrix{Union{Float64, Missing}}(undef, sum_n, period)

        sum_n_new = sum(n_new)
        #amt_new = [fill(49,n_new[1]);fill(51,n_new[2]);fill(47,n_new[3]);fill(50,n_new[4])]
        amt_new = [fill(49,n[1])]
        lam_t_new = DataFrame(onset = Vector{Union{Float64, Missing}}(undef, sum_n_new),
                    peaktime = Vector{Union{Float64, Missing}}(undef, sum_n_new),
                    peakh = Vector{Union{Float64, Missing}}(undef, sum_n_new),
                    recovery = Vector{Union{Float64, Missing}}(undef, sum_n_new))
        mu_all_new = Matrix{Union{Float64, Missing}}(undef, sum_n_new, period)
        p_all_new = Matrix{Union{Float64, Missing}}(undef, sum_n_new, period)

        for i in 1:period
            counter = 1
            counter_new = 1
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

                # p_new = rand(Normal(mu_p, sigma_p))
                # r_new = rand(Normal(mu_r, sigma_r))
                # h_new = rand(Exponential(beta_new))
                # u_new = rand(Uniform(0, p_new + r_new))
                # onset_new = T[i] - u_new

                # lam_t_new.onset[counter] = onset_new
                # lam_t_new.peaktime[counter] = onset_new + p_new
                # lam_t_new.peakh[counter] = h_new
                # lam_t_new.recovery[counter] = onset_new + p_new + r_new
                # if i == 1
                #     mu_all_new[counter, i] = lam_int(0, T[i], lam_t_new[counter, :])
                #     p_all_new[counter, i] = p_int(0, T[i], lam_t_new[counter, :], mu_all_new[counter, i,], xi)
                # else
                #     mu_all_new[counter, i] = lam_int(T[i-1], T[i], lam_t_new[counter, :])
                #     p_all_new[counter, i] = p_int(T[i-1], T[i], lam_t_new[counter, :], mu_all_new[counter, i,], xi)
                # end

                counter += 1
            end



            for j in 1:sum_n_new
                p = rand(Normal(mu_p, sigma_p))
                r = rand(Normal(mu_r, sigma_r))
                h = rand(Exponential(1/beta))
                u = rand(Uniform(0, p + r))
                onset = T[i] - u
                
                lam_t_new.onset[counter_new] = onset
                lam_t_new.peaktime[counter_new] = onset + p
                lam_t_new.peakh[counter_new] = h
                lam_t_new.recovery[counter_new] = onset + p + r
                if i == 1
                    mu_all_new[counter_new, i] = lam_int(0,T[i], lam_t_new[counter_new, :])
                    p_all_new[counter_new, i] = p_int(0, T[i], lam_t_new[counter_new, :], mu_all_new[counter_new, i], xi)
                else
                    mu_all_new[counter_new, i] = lam_int(T[i-1],T[i], lam_t_new[counter_new, :])
                    p_all_new[counter_new, i] = p_int(T[i-1], T[i], lam_t_new[counter_new, :], mu_all_new[counter_new, i], xi)
                end

                counter_new += 1
            end
        end

        #grand_mean = sum(skipmissing(amt .* mu_all .* p_all))
        #grand_mean_new = sum(skipmissing(amt_new .* mu_all_new .* p_all_new))

        p_means = vec(sum(mu_all, dims = 1) .* sum(p_all, dims = 1))
        p_means_new = vec(sum(mu_all_new, dims = 1) .* sum(p_all_new, dims = 1))

        A[l,:] = p_means_new
        
        dust_cal_adj_new = deepcopy(dust_cal)
        for j in 1:period
            if j == 1
                dust_cal_adj_new[j] = dust_cal[j]
            else
                dust_cal_adj_new[j] = dust_cal[j] - min(convert(Int64, round(.9*dust_cal[j])),rand(Binomial(dust_cal_adj_new[j-1], exp(-xi*(T[j] - T[j-1])))))
            end
        end

        loglik = 0
        for i in 1:period
            loglik = loglik + logpdf(Poisson(p_means[i]), dust_cal_adj[i])
        end

        loglik_new = 0
        for i in 1:period
            loglik_new = loglik_new + logpdf(Poisson(p_means_new[i]), dust_cal_adj[i])
        end

        #loglik = logpdf(Poisson(grand_mean), sum(dust_cal_adj))
        #loglik_new = logpdf(Poisson(grand_mean_new), sum(dust_cal_adj_new))

        ar = min(1, exp(loglik_new - loglik))
        flip = rand(Bernoulli(ar))

        if flip == 1
            n = n_new
        end

        n_samp[l,:] = n
    end

    # Visualization
    #histogram(beta_samp, bins=30, xlabel="Beta", ylabel="Frequency", title="Beta Samples")
    #plot(beta_samp, xlabel="Iteration", ylabel="Beta", title="Trace Plot", lw=2)
    df = DataFrame(Value = n_samp)
    CSV.write("n_samp_3.csv", n_samp)
    CSV.write("A_3.csv", A)
end
