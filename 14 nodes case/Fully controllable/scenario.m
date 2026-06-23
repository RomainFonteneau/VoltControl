function mpc = scenario(mpc, k)
%SCENARIO Apply the setpoints corresponding to a summer-day cycle at step k.
%   mpc = SCENARIO(mpc, k) returns the matpower case mpc with the load and
%   generation setpoints corresponding to time step k (0-indexed, 5-minute
%   time step). Step k = 0 is 00:00 and step k = 287 is 23:55, so 288
%   steps cover a full summer day. For k outside [0, 287], mpc is returned
%   unchanged.
%
%   The distribution generators at buses 12 and 13 are modeled as
%   aggregations of wind and solar production, each with its own
%   short-term fluctuation to reflect the fact that the two sites are not
%   exposed to the same weather at the same time. The generator at bus 10
%   is a hydraulic plant dispatched to follow the residual demand (reduced
%   output around solar noon, ramped up for the morning and evening
%   peaks).
%
%   Each distribution load (buses 11 to 14) is assigned a load type
%   (residential, industrial or mixed). The active power profile and the
%   power factor of each load are evaluated for the current time of day
%   from its load type, reflecting how the mix of appliances and
%   equipment behind each bus changes over the day.

define_constants;
 
 if k==100
     mpc.branch([17;18],BR_STATUS)=0;
 end

n_steps = 288;
if k < 0 || k > n_steps - 1
    return;
end

% Reference (nominal) load and generation data, used as the basis for the
% time-of-day setpoints below.
mpc0 = Extented_case();
n_bus = size(mpc0.bus, 1);

% Time of day, in hours (step k = 0 is midnight)
t = k / 12;

% --- Renewable availability profiles, in per unit of installed capacity ---

% Solar: clear-sky irradiance shape between sunrise and sunset, typical of
% a summer day in France (sunrise 05:30, sunset 21:30, solar noon at
% 13:30 local summer time). The exponent gives a broad plateau around
% solar noon rather than a single sharp peak. A derating factor accounts
% for inverter losses and the reduced efficiency of PV panels at high
% summer temperatures.
t_sunrise = 5.5;
t_sunset  = 21.5;
solar_shape_exponent = 1.5;
pv_derating = 0.85;

if t > t_sunrise && t < t_sunset
    solar_pu = pv_derating * sin(pi * (t - t_sunrise) / (t_sunset - t_sunrise)) ^ solar_shape_exponent;
else
    solar_pu = 0;
end

% Wind: mild daily variation, slightly stronger overnight and weaker in
% the afternoon, on top of which a slow smooth fluctuation represents the
% short-term variability of the wind resource. The two wind farms (buses
% 12 and 13) use independent fluctuations, since they are not exposed to
% the same gusts and weather fronts at the same time.
wind_daily = 0.40 + 0.10 * cos(2*pi * (t - 3) / 24);
wind_pu_12 = min(1, max(0, wind_daily + 0.3 * smooth_fluctuation(t, 0.0)));
wind_pu_13 = min(1, max(0, wind_daily + 0.3 * smooth_fluctuation(t, 4.7)));

% Solar production is also affected by passing clouds, independently at
% each of the two sites.
solar_pu_12 = solar_pu * (1 - 0.25 * max(0, smooth_fluctuation(t, 2.3)));
solar_pu_13 = solar_pu * (1 - 0.25 * max(0, smooth_fluctuation(t, 6.1)));

% --- Distribution generators (wind + solar aggregation) -----------------

% Bus 12: solar-dominant mix (55 MW PV + 20 MW wind installed capacity)
mpc = set_gen_power(mpc, 12, 55*solar_pu_12 + 20*wind_pu_12);

% Bus 13: more balanced mix (30 MW PV + 40 MW wind installed capacity)
mpc = set_gen_power(mpc, 13, 30*solar_pu_13 + 40*wind_pu_13);

% --- Hydraulic plant (bus 10), dispatched to follow residual demand -----
% Output is reduced around solar noon, when the PV plants cover most of
% the demand, and ramped up for the morning and evening peaks. The dip
% and the two ramps are based on the clear-sky solar profile rather than
% on the cloud-affected output of the individual PV sites, reflecting a
% dispatch decision made ahead of time on the basis of the expected solar
% resource for the whole area.
g10 = find(mpc0.gen(:, GEN_BUS) == 10);
hydro_capacity = mpc0.gen(g10, PMAX);

morning_ramp = exp(-((t - 8)/2)^2);
evening_ramp = exp(-((t - 20)/2)^2);
hydro_pu = 0.5 - 0.35*solar_pu + 0.20*morning_ramp + 0.35*evening_ramp;

mpc = set_gen_power(mpc, 10, hydro_capacity * hydro_pu);

% --- Distribution loads (buses 11, 12, 13, 14) ---------------------------
% Each load bus is associated with a load type: residential, industrial,
% or mixed. The active power and power factor profiles for each type are
% evaluated at the current time of day and applied on top of the nominal
% PD of mpc0; QD is derived from the resulting Pd and power factor.

load_type = zeros(n_bus, 1);
load_type(11) = 1; % residential
load_type(12) = 3; % mixed residential/commercial/industrial
load_type(13) = 3; % mixed residential/commercial/industrial
load_type(14) = 2; % industrial

for bus = 1:n_bus
    if load_type(bus) == 0
        continue;
    end

    [load_pu, cos_phi] = load_profile(load_type(bus), t);

    Pd = mpc0.bus(bus, PD) * load_pu;
    mpc.bus(bus, PD) = Pd;
    mpc.bus(bus, QD) = Pd * tan(acos(cos_phi));
end

end

function mpc = set_gen_power(mpc, bus, p_value)
%SET_GEN_POWER Set the active power setpoint of the generator connected to
%bus, clipped to its [PMIN, PMAX] range.

define_constants;

g = find(mpc.gen(:, GEN_BUS) == bus);
mpc.gen(g, PG) = min(mpc.gen(g, PMAX), max(mpc.gen(g, PMIN), p_value));

end

function y = smooth_fluctuation(t, phase_offset)
%SMOOTH_FLUCTUATION Smooth, quasi-random fluctuation around zero, with
%values roughly within [-1, 1], evaluated at time of day t (in hours).
%phase_offset shifts the fluctuation in time, so that calling this
%function with different offsets gives fluctuations that look
%independent of each other.
%
%   The fluctuation is built as a sum of cosines with periods that are
%   not simple multiples of each other or of 24 hours, so that the
%   resulting signal looks irregular over the course of a day while
%   remaining smooth and deterministic.

periods    = [5.3, 11.7, 2.1, 17.9]; % hours
amplitudes = [0.45, 0.30, 0.15, 0.10];
phases     = [0.7, 2.3, 4.1, 1.0];

y = 0;
for i = 1:length(periods)
    y = y + amplitudes(i) * cos(2*pi * t / periods(i) + phases(i) + phase_offset);
end

end

function [load_pu, cos_phi] = load_profile(type, t)
%LOAD_PROFILE Active power multiplier and power factor for a given load
%type at time of day t (in hours, 0 = midnight).
%
%   type = 1: residential load. Active power follows a night-time dip, a
%   morning peak, a small midday peak separated from the morning peak by
%   a mid-morning dip, an afternoon dip, and the main peak in the evening.
%   The power factor is close to unity for most of the day and degrades
%   in the afternoon and evening, when air-conditioning and other
%   motor-driven appliances are heavily used.
%
%   type = 2: industrial load. Active power rises in the morning to a
%   plateau, dips around lunchtime, rises again to a second plateau in
%   the afternoon, then falls back to its night-time level in the
%   evening. The power factor is lower overall and improves during the
%   two production plateaus, when machinery operates closer to its rated
%   load.
%
%   type = 3: mixed load (residential, commercial and light industrial).
%   Active power has a moderate night-time dip, rises in the morning to a
%   broad daytime level with a small midday bump, and peaks in the late
%   afternoon and evening. The power factor stays close to its nominal
%   value, with a small improvement during business hours and a small
%   degradation in the evening.

baseline = 0.40;

switch type
    case 1 % residential
        load_pu = baseline ...
            - 0.15 * gaussian_bump(t,  3.0, 2.0) ...
            + 0.30 * gaussian_bump(t,  7.5, 1.0) ...
            - 0.12 * gaussian_bump(t, 10.5, 1.3) ...
            + 0.18 * gaussian_bump(t, 12.5, 1.0) ...
            - 0.12 * gaussian_bump(t, 16.0, 1.5) ...
            + 0.50 * gaussian_bump(t, 20.0, 1.2);

        cos_phi = 0.97 ...
            - 0.04 * gaussian_bump(t, 16.0, 2.0) ...
            - 0.05 * gaussian_bump(t, 20.0, 1.5);

    case 2 % industrial
        load_pu = 0.5*baseline ...
            + 0.60 * gaussian_bump(t,  10.5, 2.5) ...
            - 0.10 * gaussian_bump(t, 12.5, 1.0) ...
            + 0.60 * gaussian_bump(t, 15.0, 3.0);

        cos_phi = 0.80 ...
            + 0.08 * gaussian_bump(t,  9.5, 3.0) ...
            + 0.08 * gaussian_bump(t, 15.5, 3.0);

    case 3 % mixed
        load_pu = baseline ...
            - 0.10 * gaussian_bump(t,  3.5, 2.0) ...
            + 0.35 * gaussian_bump(t,  9.0, 2.5) ...
            + 0.10 * gaussian_bump(t, 12.5, 1.2) ...
            + 0.40 * gaussian_bump(t, 18.5, 2.5);

        cos_phi = 0.95 ...
            + 0.02 * gaussian_bump(t,  9.0, 3.0) ...
            - 0.02 * gaussian_bump(t, 18.5, 3.0);

    otherwise
        load_pu = 1;
        cos_phi = 1;
end

end

function y = gaussian_bump(t, t0, sigma)
%GAUSSIAN_BUMP Unit-height Gaussian bump centered on t0 with standard
%deviation sigma, evaluated at t.

y = exp(-((t - t0) / sigma)^2);

end