function [ yyD, yyR ] = normalizedValuesAsFuncRad(hexmat, yout0, R, k)
    % Rad = 6; k = 127;
    Rad = R;
    final = yout0(end, :);
    D = final(1:k); % levels of Delta in cells 1 to k
    R = final(k+1:2*k); % levels of Repressor in cells 1 to k
    N = final(2*k+1:3*k); % levels of Notch in cells 1 to k

    new_hexmat = hexmat(:, :, 1);
    new_hexmat(:, :, 2) = 0;

    % find radius
    new_hexmat(Rad+1, Rad+1, 2) = 0; % 0th
    for i=1:(2*Rad+1)
        for j=1:(2*Rad+1)
            if hexmat(i, j, 2) ~= 0
                new_hexmat(i, j, 2) = Rad; % Last radius
            end
        end
    end
    for r=Rad-1:-1:1
        for i=1:(2*Rad+1)
            for j=1:(2*Rad+1)
                if (i > 1 && i < 2*Rad+1 && j > 1 && j < 2*Rad+1) && ...
                   hexmat(i, j, 1) ~= 0 && hexmat(i, j, 2) == 0 && new_hexmat(i, j, 2) == 0 && ...
                  ((new_hexmat(i+1, j, 2) == (r+1) && (new_hexmat(i-1, j, 2) == 0 || new_hexmat(i-1, j, 2) == r)) || ... 
                   (new_hexmat(i-1, j, 2) == (r+1) && (new_hexmat(i+1, j, 2) == 0 || new_hexmat(i+1, j, 2) == r)) || ...
                   (new_hexmat(i, j+1, 2) == (r+1) && (new_hexmat(i, j-1, 2) == 0 || new_hexmat(i, j-1, 2) == r)) || ...
                   (new_hexmat(i, j-1, 2) == (r+1) && (new_hexmat(i, j+1, 2) == 0 || new_hexmat(i, j+1, 2) == r)))
                    new_hexmat(i, j, 2) = r;
                end
            end
        end
    end
    % Fix middle column
    for i=1:Rad
        if mod(i, 2) == 0
            new_hexmat(Rad+1-i, Rad+1, 2) = i;
            new_hexmat(Rad+1+i, Rad+1, 2) = i;
        end
    end

    radiusSums_D = zeros(1, Rad+1);
    radiusSums_R = zeros(1, Rad+1);
    numInRad = zeros(1, Rad+1);
    for i=1:(2*Rad+1)
        for j=1:(2*Rad+1)
            if new_hexmat(i, j, 1) ~= 0 % has an index
                radius = new_hexmat(i, j, 2);
                numInRad(1, radius+1) = numInRad(1, radius+1)+1;
                radiusSums_D(1, radius+1) = radiusSums_D(1, radius+1) + D(1, new_hexmat(i, j, 1));
                radiusSums_R(1, radius+1) = radiusSums_R(1, radius+1) + R(1, new_hexmat(i, j, 1));
            end
        end
    end
    avg_radiusSums_D = radiusSums_D./numInRad;
    avg_radiusSums_R = radiusSums_R./numInRad;
 
    tminD = 0.3; tmaxD = 1; rminD = min(avg_radiusSums_D); rmaxD = max(avg_radiusSums_D);
    norm_radiusSums_D = (((avg_radiusSums_D - rminD)./(rmaxD - rminD)) * (tmaxD-tminD)) + tminD;
    norm_radiusSums_D = [ flip(norm_radiusSums_D) norm_radiusSums_D(1, 2:end) ];
    
    tminR = 0.3; tmaxR = 1; rminR = min(avg_radiusSums_R); rmaxR = max(avg_radiusSums_R);
    norm_radiusSums_R = (((avg_radiusSums_R - rminR)./(rmaxR - rminR)) * (tmaxR-tminR)) + tminR;
    norm_radiusSums_R = [ flip(norm_radiusSums_R) norm_radiusSums_R(1, 2:end) ];
    
    norm_radiusSums_D = avg_radiusSums_D./max(avg_radiusSums_D);
    norm_radiusSums_D = [ flip(norm_radiusSums_D) norm_radiusSums_D(1, 2:end) ];
    norm_radiusSums_R = avg_radiusSums_R./max(avg_radiusSums_R);
    norm_radiusSums_R = [ flip(norm_radiusSums_R) norm_radiusSums_R(1, 2:end) ];
    x = [-6:1:0 1:6]./6; xx = -1:0.1:1;
    pD = polyfit(x, norm_radiusSums_D, 5);
    yyD = polyval(pD, xx);
    pR = polyfit(x, norm_radiusSums_R, 5);
    yyR = polyval(pR, xx);
end

function drawSecondGraph(norm_receptor_high, norm_receptor_low, norm_rep_high, norm_rep_low)
    figure;
    subplot(1,2,1)
    plot(1:13, norm_receptor_high, x, norm_receptor_low, 'LineWidth', 2.0)
    xlim([-1 1])
    ylim([0 1])
    xlabel('Dimensionless Radius')
    ylabel('Concentration')
    title('Simulated Receptor Repressor Concentration')
    legend('High (b=5)','Low (b=0)')
    grid on

    subplot(1,2,2)
    plot(x, norm_rep_high, x, norm_rep_low, 'LineWidth', 2.0)
    xlim([-1 1])
    ylim([0 1])
    xlabel('Dimensionless Radius')
    ylabel('Concentration')
    title('Simulated Signaling Repressor Concentration')
    legend('High (b=5)','Low (b=0)')
    grid on
end