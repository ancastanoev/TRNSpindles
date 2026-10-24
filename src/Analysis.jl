"""
    Analysis.jl

Signal-analysis utilities for TRN spindle research.

Functions
─────────
- `dominant_frequency`       — FFT-based dominant frequency of a time series
- `spindle_proximity_score`  — scalar score: how close is a cycle to a spindle?
- `classify_regime`          — label a solution as N2/transition/N3
- `extract_peaks`            — simple peak-detection on a univariate signal
- `mean_cycle_period`        — average inter-peak period
- `autocorrelation_freq`     — frequency via zero-crossing of autocorrelation
"""

# Regime enum

"""
    SpindleRegime

Classification of a dynamical regime based on dominant frequency.

| Regime        | Frequency     | Sleep stage analogy |
|---------------|---------------|---------------------|
| N2_SPINDLE    | 7–15 Hz       | N2 sleep            |
| TRANSITION    | 1–7 Hz        | transition          |
| N3_SLOW_WAVE  | 0.1–1 Hz      | N3 / slow-wave      |
| INFRASLOW     | < 0.1 Hz      | infraslow (< 0.1 Hz)|
| UNKNOWN       | unclassifiable|                     |
"""
@enum SpindleRegime begin
    N2_SPINDLE
    TRANSITION
    N3_SLOW_WAVE
    INFRASLOW
    UNKNOWN
end

# Dominant frequency via FFT

"""
    dominant_frequency(signal, dt; fmin=0.01, fmax=20.0) → Float

Return the dominant frequency (Hz) of `signal` sampled at interval `dt` (s).
Only frequencies in [fmin, fmax] Hz are considered.

Uses a Hann window to reduce spectral leakage.
"""
function dominant_frequency(
    signal::AbstractVector{F}, dt::F; fmin::F=F(0.01), fmax::F=F(20.0)
)::F where {F<:AbstractFloat}
    n = length(signal)
    win = DSP.hanning(n)
    y = win .* (signal .- mean(signal))
    Y = abs.(rfft(y))
    freqs = rfftfreq(n, 1/dt)

    mask = (freqs .>= fmin) .& (freqs .<= fmax)
    if !any(mask)
        return zero(F)
    end
    idx = argmax(Y[mask])
    return freqs[mask][idx]
end

# Peak detection

"""
    extract_peaks(signal; min_prominence=0.1) → Vector{Int}

Return indices of local maxima in `signal`.
A peak at index i requires signal[i] > signal[i-1] and signal[i] > signal[i+1].
"""
function extract_peaks(
    signal::AbstractVector{F}; min_prominence::F=F(0.1)
) where {F<:AbstractFloat}
    peaks = Int[]
    for i in 2:(length(signal) - 1)
        if signal[i] > signal[i - 1] && signal[i] > signal[i + 1]
            push!(peaks, i)
        end
    end
    # filter by prominence (height above global minimum)
    global_min = minimum(signal)
    filter!(
        i -> (signal[i] - global_min) >= min_prominence * (maximum(signal) - global_min),
        peaks,
    )
    return peaks
end

"""
    mean_cycle_period(t, peaks) → Float

Mean inter-peak period (in model time units) from a list of peak indices.
"""
function mean_cycle_period(
    t::AbstractVector{F}, peaks::Vector{Int}
)::F where {F<:AbstractFloat}
    (isempty(peaks) || length(peaks) < 2) && return zero(F)
    diffs = diff(t[peaks])
    return mean(diffs)
end

# Autocorrelation-based frequency

"""
    autocorrelation_freq(signal, dt; max_lag=nothing) → Float

Estimate dominant frequency via the first zero-crossing lag of the normalised
autocorrelation function.  Frequency = 1 / (2 · lag).
"""
function autocorrelation_freq(
    signal::AbstractVector{F}, dt::F; max_lag::Union{Nothing,Int}=nothing
)::F where {F<:AbstractFloat}
    n = length(signal)
    lim = isnothing(max_lag) ? n ÷ 2 : min(max_lag, n-1)
    y = signal .- mean(signal)
    σ2 = var(y)
    σ2 < eps(F) && return zero(F)

    acf = [dot(y[1:(n - k)], y[(1 + k):n]) / ((n-k) * σ2) for k in 0:lim]

    # first zero crossing
    for k in 1:lim
        if acf[k] >= 0 && acf[k + 1] < 0
            # linear interpolation
            lag = (k + acf[k] / (acf[k] - acf[k + 1])) * dt
            return 1 / (2lag)
        end
    end
    return zero(F)
end

# Spindle proximity score

"""
    spindle_proximity_score(Tmax, Rmax;
                            T_target=(4.0,6.0), R_target=(10.0,15.0)) → Float ∈ [0,1]

Score how close a limit-cycle peak is to the spindle-compatible target windows
for the thalamic relay (T) and TRN (R) populations.

Returns 1.0 for a perfect biological spindle, 0.0 for an escape state.

Algorithm
─────────
d_T = max(T_lo − Tmax, 0, Tmax − T_hi) / T_hi
d_R = max(R_lo − Rmax, 0, Rmax − R_hi) / R_hi
score = exp(−(d_T² + d_R²) / 2σ²)   with σ = 0.3
"""
function spindle_proximity_score(
    Tmax::F,
    Rmax::F;
    T_target::Tuple{F,F}=(F(4), F(6)),
    R_target::Tuple{F,F}=(F(10), F(15)),
    σ::F=F(0.3),
)::F where {F<:AbstractFloat}
    T_lo, T_hi = T_target
    R_lo, R_hi = R_target

    dT = max(T_lo - Tmax, zero(F), Tmax - T_hi) / T_hi
    dR = max(R_lo - Rmax, zero(F), Rmax - R_hi) / R_hi

    return exp(-(dT^2 + dR^2) / (2σ^2))
end

"""
    joint_spindle_score(ρ, Rmin, Rmax;
                        ρ_target=(0.3,0.6),
                        R_target=(10.0,15.0),
                        Rmin_target=(0.0,4.0)) → Float ∈ [0,1]

Joint spindle score incorporating the delay parameter ρ, the TRN oscillation
amplitude Rmax, and the TRN minimum Rmin.
"""
function joint_spindle_score(
    ρ::F,
    Rmin::F,
    Rmax::F;
    ρ_target::Tuple{F,F}=(F(0.3), F(0.6)),
    R_target::Tuple{F,F}=(F(10), F(15)),
    Rmin_target::Tuple{F,F}=(F(0), F(4)),
) where {F<:AbstractFloat}
    in_ρ = ρ_target[1] <= ρ <= ρ_target[2]
    in_Rmax = R_target[1] <= Rmax <= R_target[2]
    in_Rmin = Rmin_target[1] <= Rmin <= Rmin_target[2]

    n_satisfied = in_ρ + in_Rmax + in_Rmin
    return F(n_satisfied) / 3
end

# Regime classifier

"""
    classify_regime(freq; spindle_band=(7.0,15.0), slow_band=(0.1,1.0),
                    infraslow_max=0.1) → SpindleRegime

Classify a dominant frequency into a sleep-stage regime.

Regime bands:
- small delays → spindle-compatible / N2-like   (7–15 Hz)
- intermediate  → transition regime              (1–7 Hz)
- large delays  → slow-wave / N3-like            (0.1–1 Hz)
- very slow     → infraslow                      (< 0.1 Hz)
"""
function classify_regime(
    freq::F;
    spindle_band::Tuple{F,F}=(F(7), F(15)),
    transition_lo::F=F(1.0),
    slow_band::Tuple{F,F}=(F(0.1), F(1.0)),
)::SpindleRegime where {F<:AbstractFloat}
    freq <= 0 && return UNKNOWN
    freq >= spindle_band[1] && freq <= spindle_band[2] && return N2_SPINDLE
    freq >= transition_lo && freq < spindle_band[1] && return TRANSITION
    freq >= slow_band[1] && freq < transition_lo && return N3_SLOW_WAVE
    freq < slow_band[1] && return INFRASLOW
    return UNKNOWN
end

# High-level analysis of a SimResult

"""
    analyze_timeseries(res, pop_idx, dt;
                       transient_frac=0.5) → NamedTuple

Analyse the time series for population index `pop_idx` from a `SimResult`.

Returns
───────
A NamedTuple with fields:
- `freq_fft`     — FFT dominant frequency
- `freq_acf`     — autocorrelation frequency
- `freq_peaks`   — peak-based frequency
- `regime`       — `SpindleRegime`
- `amplitude`    — peak-to-trough amplitude
- `mean_val`     — temporal mean
"""
function analyze_timeseries(
    res::SimResult{F}, pop_idx::Int, dt::F; transient_frac::F=F(0.5)
)::NamedTuple where {F<:AbstractFloat}
    n_skip = round(Int, transient_frac * size(res.u, 2))
    sig = res.u[pop_idx, (n_skip + 1):end]
    t_sig = res.t[(n_skip + 1):end]

    f_fft = dominant_frequency(sig, dt)
    f_acf = autocorrelation_freq(sig, dt)
    peaks = extract_peaks(sig)
    f_peak = if length(peaks) >= 2
        1 / mean_cycle_period(t_sig, peaks)
    else
        zero(F)
    end
    reg = classify_regime(f_fft)
    amp = maximum(sig) - minimum(sig)
    μ = mean(sig)

    return (
        freq_fft=f_fft,
        freq_acf=f_acf,
        freq_peaks=f_peak,
        regime=reg,
        amplitude=amp,
        mean_val=μ,
    )
end
