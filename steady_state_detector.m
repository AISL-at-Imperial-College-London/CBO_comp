function [sys, x0, str, ts] = steady_state_detector(t, x, u, flag)

persistent buf_y buf_len ss_state ss_count

% ==========================
% User parameters
% ==========================
Ts        = 1;      % sample time (s)
Nwin      = 20;     % window length (samples)
Mconfirm  = 5;      % require conditions to hold this many consecutive checks


th_range_enter = 0.015;
th_slope_enter = 5e-4;
th_mad_enter   = 0.004;


th_range_exit = 1.5 * th_range_enter;
th_slope_exit = 1.5 * th_slope_enter;
th_mad_exit   = 1.5 * th_mad_enter;

switch flag

    % --------------------------------------
    % Initialization
    % --------------------------------------
    case 0
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = 0;
        sizes.NumOutputs     = 1;
        sizes.NumInputs      = 1;   % single measurement y
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        sys = simsizes(sizes);

        x0  = [];
        str = [];
        ts  = [Ts 0];

        buf_y    = NaN(Nwin,1);
        buf_len  = 0;
        ss_state = 0;
        ss_count = 0;

    % --------------------------------------
    % Update at each sample
    % --------------------------------------
    case 2
        y = u(1);

        % Guard against NaNs/Infs
        if ~isfinite(y)
            ss_state = 0;
            ss_count = 0;
            sys = [];
            return;
        end

        % Update buffer
        if buf_len < Nwin
            buf_len = buf_len + 1;
            buf_y(buf_len) = y;
        else
            buf_y(1:end-1) = buf_y(2:end);
            buf_y(end)     = y;
        end

        % Only evaluate when window is full
        if buf_len == Nwin && all(isfinite(buf_y))
            ywin = buf_y;

            y_range = max(ywin) - min(ywin);

            y_med = median(ywin);
            y_mad = 1.4826 * median(abs(ywin - y_med));  

            tt = (0:Nwin-1)' * Ts;
            p  = polyfit(tt, ywin, 1);   
            y_slope = abs(p(1));

            if ss_state == 0
                ok = (y_range < th_range_enter) && (y_mad < th_mad_enter) && (y_slope < th_slope_enter);
            else
                ok = (y_range < th_range_exit)  && (y_mad < th_mad_exit)  && (y_slope < th_slope_exit);
            end

            if ok
                ss_count = ss_count + 1;
            else
                ss_count = 0;
                ss_state = 0;   
            end

            if ss_state == 0 && ss_count >= Mconfirm
                ss_state = 1;
            end

        else
            ss_state = 0;
            ss_count = 0;
        end

        sys = [];

    % --------------------------------------
    % Output
    % --------------------------------------
    case 3
        sys = double(ss_state);

    otherwise
        sys = [];
end
end
