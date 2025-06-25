function plotdvh(d0, var)
    % PLOTDVH Plots the DVH (Dose-Volume Histogram) for different structures.
    %
    % INPUT:
    %   d0   - Cell array of dose distributions for different structures.
    %   var  - Structure containing relevant information.
    %
    % OUTPUT:
    %   None

    n_c = var.n_c;
    dmax = var.dmax;

    n = 100;
    t = linspace(0, dmax, n);
    ptv1 = zeros(n, 1);
    body = zeros(n, 1);
    % bladder = zeros(n, 1);
    % rectum = zeros(n, 1);
    % femhead_lt = zeros(n, 1);
    % femhead_rt = zeros(n, 1);
    % penilebulb = zeros(n, 1);

    for i = 1:n
        ptv1(i) = numel(find(d0{1} >= t(i))) / n_c(1);
        body(i) = numel(find(d0{2} >= t(i))) / n_c(2);
        % bladder(i) = numel(find(d0{3} >= t(i))) / n_c(3);
        % rectum(i) = numel(find(d0{4} >= t(i))) / n_c(4);
        % femhead_lt(i) = numel(find(d0{5} >= t(i))) / n_c(5);
        % femhead_rt(i) = numel(find(d0{6} >= t(i))) / n_c(6);
        % penilebulb(i) = numel(find(d0{7} >= t(i))) / n_c(7);
    end

    figure;
    hold on;
    plot(t, ptv1, 'r');
    plot(t, body, 'k');
    % plot(t, bladder, 'b');
    % plot(t, rectum, 'b');
    % plot(t, femhead_lt, 'g');
    % plot(t, femhead_rt, 'g');
    % plot(t, penilebulb, 'y');
    hold off;
    drawnow;
end