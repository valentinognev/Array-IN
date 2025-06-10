function [f,ax] = plot_acc_pos(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end
    
    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end

    plotOpt.variable_name_mean = "r";
    plotOpt.variable_name_var = "r_std";
    plotOpt.transformation = @(x) x*1e3;
    plotOpt.title_label = "Accelerometer pos. 1";
    plotOpt.y_unit = "[mm]";
    plotOpt.pred_type = "filt";

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end


end
function [f,ax] = plot_angular_acceleration(data, t, cases_plot, caseOpt, plotType)

    if nargin < 5
        plotType = "normal";
    end

    plotOpt = struct;
    plotOpt.variable_name_mean = "omega_dot";
    plotOpt.variable_name_var = "omega_dot_std";
    plotOpt.transformation = @rad2deg;
    plotOpt.title_label = "Angular Acceleration";
    plotOpt.y_unit = "[deg/s^2]";
    plotOpt.pred_type = "pred";

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end

end
function [f,ax] = plot_angular_velocity(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end
    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end


    plotOpt.variable_name_mean = "w";
    plotOpt.variable_name_var = "w_std";
    plotOpt.transformation = @rad2deg;
    plotOpt.title_label = "Angular Velocity";
    plotOpt.y_unit = "[deg/s]";
    plotOpt.pred_type = "filt";

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "filt_and_pred"
        [f,ax] = plot_three_components_filt_pred(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end

end
function [f, ax] = plot_angular_velocity_filt_and_pred(data, t, cases_plot, caseOpt)


    plotOpt = struct;
    plotOpt.variable_name_mean = "w";
    plotOpt.variable_name_var = "w_std";
    plotOpt.transformation = @rad2deg;
    plotOpt.title_label = "Angular Velocity";
    plotOpt.y_unit = "[deg/s]";
    plotOpt.pred_type = "filt";

    f = figure;
    ax = zeros(3,1);
    for i = 1:3
        ax(i) = subplot(3,1,i);
        hold on;
        grid on;
    end

    plot_three_components_update(ax, data, t, cases_plot, caseOpt, plotOpt);


    plotOpt = struct;
    plotOpt.variable_name_mean = "w";
    plotOpt.variable_name_var = "w_std";
    plotOpt.transformation = @rad2deg;
    plotOpt.title_label = "Angular Velocity";
    plotOpt.y_unit = "[deg/s]";
    plotOpt.pred_type = "pred";

    plot_three_components_update(ax, data, t, cases_plot, caseOpt, plotOpt);


end
function [f,ax] = plot_bias_accelerometer(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end
    
    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end

    plotOpt.variable_name_mean = "b_a";
    plotOpt.variable_name_var = "b_a_std";
    plotOpt.transformation = @(x) x;
    plotOpt.title_label = "Bias accelerometer";
    plotOpt.y_unit = "[m/s^2]";
    plotOpt.pred_type = "filt";

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end


end
function [f,ax] = plot_bias_gyroscopes(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end

    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end


    plotOpt.variable_name = "b_g";
    plotOpt.transformation = @rad2deg;
    if isfield(plotOpt, "title")
        plotOpt.title_label = sprintf("%s, Bias gyroscopes",plotOpt.title);
    else
        plotOpt.title_label = "Bias gyroscopes";
    end
    plotOpt.y_unit = "[deg/s]";
    if ~isfield(plotOpt, "path")       
        plotOpt.path = ["filt"];
    end

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end


end
function [f,ax] = plot_bias_omega_dot(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end

    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end


    plotOpt.variable_name = "b_omega_dot";
    plotOpt.transformation = @rad2deg;
    if isfield(plotOpt, "title")
        plotOpt.title_label = sprintf("%s, Bias omega dot",plotOpt.title);
    else
        plotOpt.title_label = "Bias omega dot";
    end
    plotOpt.y_unit = "[deg/s^2]";    
    if ~isfield(plotOpt, "path")       
        plotOpt.path = ["filt"];
    end


    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end


end
function [f,ax] = plot_bias_s(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end
    
    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end

    plotOpt.variable_name = "b_s";
    plotOpt.transformation = @(x) x;
    if isfield(plotOpt, "title")
        plotOpt.title_label = sprintf("%s, Bias specific force",plotOpt.title);
    else
        plotOpt.title_label = "Bias specific force";
    end
    plotOpt.y_unit = "[m/s^2]";
    if ~isfield(plotOpt, "path")       
        plotOpt.path = ["filt"];
    end

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end


end
function [f,ax] = plot_navigation_acceleration(data, t, cases_plot, caseOpt, plotType)

    if nargin < 5
        plotType = "normal";
    end

    plotOpt = struct;
    plotOpt.variable_name_mean = "v_dot";
    plotOpt.variable_name_var = "v_dot_std";
    plotOpt.transformation = @(x) x;
    plotOpt.title_label = "Navigation Acceleration";
    plotOpt.y_unit = "[m/s^2]";
    plotOpt.pred_type = "pred";

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end

end
function [f,ax] = plot_navigation_position(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end

    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end


    plotOpt.variable_name = "p";
    plotOpt.transformation = @(x) x;
    if isfield(plotOpt, "title")
        plotOpt.title_label = sprintf("%s ",plotOpt.title);
    else
        plotOpt.title_label = "";
    end

    plotOpt.title_label = plotOpt.title_label + "Navigation position";
    plotOpt.y_unit = "[m]";
    if ~isfield(plotOpt, "path")       
        plotOpt.path = ["filt"];
    end


    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end


end
function [f,ax] = plot_navigation_position_v2(data, t,  caseOpt, plotType, extras)

    if nargin < 4
        plotType = "normal";
    end

    if nargin == 5
        plotOpt = extras;
        fprintf("Using Extras\n")
    else
        plotOpt = struct;
    end


    plotOpt.variable_name_mean = "p";
    plotOpt.variable_name_var = "p_std";
    plotOpt.transformation = @(x) x;
    if isfield(plotOpt, "title")
        plotOpt.title_label = sprintf("%s ",plotOpt.title);
    else
        plotOpt.title_label = "";
    end

    plotOpt.title_label = plotOpt.title_label + "Navigation position";
    plotOpt.y_unit = "[m]";
    plotOpt.pred_type = "filt";


    if plotType == "error"
        [f,ax] = plot_three_components_error_v2(data, t,  caseOpt, plotOpt);
    else
        error("wrong")
    end


end
function [f,ax] = plot_navigation_velocity(data, t, cases_plot, caseOpt, plotType , extras)

    if nargin < 5
        plotType = "normal";
    end

    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end


    plotOpt.variable_name_mean = "v";
    plotOpt.variable_name_var = "v_std";
    plotOpt.transformation = @(x) x;
    plotOpt.title_label = "Navigation Velocity";
    plotOpt.y_unit = "[m/s]";
    plotOpt.pred_type = "filt";

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end


end
function [f,ax] = plot_nine_components(dataTot, t, cases_plot, caseOpt, plotOpt)

pred_type = plotOpt.pred_type;
v_name = plotOpt.variable_name_mean;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;

if isfield(plotOpt,"components")
    directions = plotOpt.components;
else
    directions = reshape(["(1,1)","(1,2)","(1,3)";
        "(2,1)","(2,2)","(2,3)";
        "(3,1)","(3,2)","(3,3)"],[],1);
end

if isfield(plotOpt,"show_cov")
    show_cov = plotOpt.show_cov; 
else
    show_cov = true;
end

if isfield(plotOpt,"figure")
    f = plotOpt.figure; 
else
    f = figure();
end

if isfield(plotOpt,"axes")
    ax = plotOpt.axes; 
else
    ax = zeros(9,1);
    for i = 1:9
        ax(i) = subplot(3,3,i); hold on;
    end
end

for i = 1:9

    leg = zeros(2,length(cases_plot));
    for i_c = 1:length(cases_plot)
        c = cases_plot(i_c);
        res = dataTot.(c);
        res_pred = dataTot.(c).(pred_type);

        if isfield(caseOpt,c) && isfield(caseOpt.(c),"mean")
            resLineSpecs = get_ls_specs(caseOpt.(c).mean);
        else
            resLineSpecs = get_ls_specs(struct);
        end        
        if ~isfield(res_pred.mean,v_name)
            error("%s is not in %s",v_name, c)
        end
        data = f_trans(res_pred.mean.(v_name));
        leg(1,i_c) = plot(ax(i), t, data(i,:), ...
            "DisplayName", res.label, ...
            "Color", colorOrder(i_c,:), ...
            "Linestyle", resLineSpecs.ls_c, ...
            "Marker", resLineSpecs.m_c, ...
            "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

        if show_cov && isfield(res_pred.std, v_name)
            data_sig = f_trans(res_pred.std.(v_name));
            if isfield(caseOpt,c) && isfield(caseOpt.(c),"cov")
                resLineSpecs = get_ls_specs(caseOpt.(c).cov);
            else
                resLineSpecs = get_ls_specs(struct);
                resLineSpecs.ls_c = "--";
            end

            leg(2,i_c) = plot(ax(i), t, data(i,:) + 3*data_sig(i,:), ...
                "DisplayName", sprintf("%s $3\\sigma$", res.label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

            plot(ax(i), t, data(i,:) - 3*data_sig(i,:), ...
                "DisplayName", res.label, ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        else
            leg(2,i_c) = nan;
        end
    end
    grid(ax(i),"on")

    title(ax(i),title_label)
    if isfield(plotOpt,"legend")
        leg_opt = plotOpt.legend;
    else
        leg_opt = {"Location","best","NumColumns",1, 'interpreter',"latex"};
    end
    legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
    xlabel(ax(i),"Time [s]")
    ylabel(ax(i),sprintf("%s %s", directions(i), y_unit))
end

end
function [f,ax] = plot_rotation(data, t, cases_plot, caseOpt, plotType, extras)
% Only works for error

if nargin == 6
    plotOpt = extras;
else
    plotOpt = struct;
end


plotOpt.variable_name = "R";
plotOpt.transformation = @rad2deg;
plotOpt.y_unit = "[deg]";

if ~isfield(plotOpt, "path")       
    plotOpt.path = ["filt"];
end

if isfield(plotOpt, "title")
    plotOpt.title_label = sprintf("%s ",plotOpt.title);
else
    plotOpt.title_label = " ";
end

if nargin < 5
    plotType = "normal";
end

if plotType == "error"
    plotOpt.title_label = plotOpt.title_label + "Rotation Error";
    [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
elseif plotType == "error_log"
    plotOpt.title_label = plotOpt.title_label + "Rotation Error";
    [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
elseif plotType == "norm_error_log"
    plotOpt.title_label = plotOpt.title_label + "Rotation Error";
    [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
else
    plotOpt.title_label = plotOpt.title_label + "Rotation";
    plotOpt.components = ["roll","pitch","yaw"];
    plotOpt.transformation = @(x) rad2deg(my_rotm2eul(x));
    plotOpt.show_cov = false;
    [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
end

end
function [f,ax] = plot_rotation_matrix(dataTot, t, cases_plot, caseOpt)

    N_end = length(t);
    colorOrder = colororder;


    f = figure();
    ax = zeros(3,3);
    for i_c = 1:length(cases_plot)
        c = cases_plot(i_c);
        res = dataTot.(c);
        R = dataTot.(c).filt.R;
        R_vec = reshape(R,9,[]);

        if isfield(caseOpt,c) && isfield(caseOpt.(c),"mean")
            resLineSpecs = get_ls_specs(caseOpt.(c).mean);
        else
            resLineSpecs = get_ls_specs(struct);
        end

        for i = 1:9
            ax(i) = subplot(3,3,i); hold on;
            plot(t, R_vec(i,:), ...
                 "DisplayName", res.label, ...
                 "Color", colorOrder(i_c,:), ...
                 "Linestyle", resLineSpecs.ls_c, ...
                 "Marker", resLineSpecs.m_c, ...
                 "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
            grid on
        end

    end

end
function [f,ax] = plot_rotation_v2(data, t, caseOpt, plotType, extras)
% Only works for error
    if nargin < 4
        plotType = "normal";
    end

    if nargin == 5
        plotOpt = extras;
        fprintf("Using Extras\n")
    else
        plotOpt = struct;
    end


    plotOpt.variable_name_mean = "R";
    plotOpt.variable_name_var = "R_std";
    plotOpt.transformation = @rad2deg;
    if isfield(plotOpt, "title")
        plotOpt.title_label = sprintf("%s ",plotOpt.title);
    else
        plotOpt.title_label = "";
    end

    plotOpt.title_label = plotOpt.title_label + "Rotation";
    plotOpt.y_unit = "[deg]";
    plotOpt.pred_type = "filt";


    if plotType == "error"
        [f,ax] = plot_three_components_error_v2(data, t,  caseOpt, plotOpt);
    else
        error("wrong")
    end

end
function [f,ax] = plot_specific_force(data, t, cases_plot, caseOpt, plotType)

    if nargin < 5
        plotType = "normal";
    end

    plotOpt = struct;
    plotOpt.variable_name_mean = "s";
    plotOpt.variable_name_var = "s_std";
    plotOpt.transformation = @(x) x;
    plotOpt.title_label = "Specific force";
    plotOpt.y_unit = "[m/s^2]";
    plotOpt.pred_type = "pred";

    if plotType == "error"
        [f,ax] = plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "error_log"
        [f,ax] = plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt);
    elseif plotType == "norm_error_log"
        [f,ax] = plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt);
    else
        [f,ax] = plot_three_components(data, t, cases_plot, caseOpt, plotOpt);
    end

end
function [f,ax] = plot_T_a(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end
    
    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end

    plotOpt.variable_name_mean = "T_a";
    plotOpt.transformation = @(x) x;
    plotOpt.title_label = "T_a";
    plotOpt.y_unit = "[-]";
    plotOpt.pred_type = "filt";

    if plotType == "normal"
        [f,ax] = plot_nine_components(data, t, cases_plot, caseOpt, plotOpt);
    else
        error("Invalid choice")
    end


end
function [f,ax] = plot_T_g(data, t, cases_plot, caseOpt, plotType, extras)

    if nargin < 5
        plotType = "normal";
    end
    
    if nargin == 6
        plotOpt = extras;
    else
        plotOpt = struct;
    end

    plotOpt.variable_name_mean = "T_g";
    plotOpt.transformation = @(x) x;
    plotOpt.title_label = "T_g";
    plotOpt.y_unit = "[-]";
    plotOpt.pred_type = "filt";

    if plotType == "normal"
        [f,ax] = plot_nine_components(data, t, cases_plot, caseOpt, plotOpt);
    else
        error("Invalid choice")
    end


end
function [f,ax] = plot_three_components(dataTot, t, cases_plot, caseOpt, plotOpt)

path_struct = plotOpt.path;
variable_name = plotOpt.variable_name;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;

if isfield(plotOpt,"components")
    directions = plotOpt.components;
else
    directions = ["x","y","z"];
end

if isfield(plotOpt,"show_cov")
    show_cov = plotOpt.show_cov; 
else
    show_cov = true;
end

if isfield(plotOpt,"figure")
    f = plotOpt.figure; 
else
    f = figure();
end

if isfield(plotOpt,"axes")
    ax = plotOpt.axes; 
else
    ax = zeros(3,1);
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;
    end
end

for i = 1:3
    leg = nan(2,size(cases_plot,1));
    for i_c = 1:size(cases_plot,1)
        case_path = cases_plot(i_c, :);
        c = case_path(end);
        path_struct_i = [case_path, path_struct];
        path_struct_i_label = [case_path, path_struct(1:end-1),"label"];
        try 
            path_struct_i_struct = num2cell(path_struct_i);
            res = getfield(dataTot, path_struct_i_struct{:});
        catch 
            warning("%s is not in %s",c, join([path_struct{:}],"/"))
            continue
        end
        try
            path_struct_i_label_struct = num2cell(path_struct_i_label);
            label_2 = getfield(dataTot, path_struct_i_label_struct{:});
        catch 
            warning("label %s is not in %s",c, join([path_struct{:}],"/"))
            label_2 = "";
        end    
        try
            label_1 = dataTot.(case_path(1)).("label");
        catch
            label_1 = "";
        end
        label = label_2 + newline + label_1;
        

        if isfield(caseOpt,c) && isfield(caseOpt.(c),"mean")
            resLineSpecs = get_ls_specs(caseOpt.(c).mean);
        else
            resLineSpecs = get_ls_specs(struct);
        end
        data = f_trans(res.("mean").(variable_name));
        leg(1,i_c) = plot(ax(i), t, data(i,:), ...
            "DisplayName", label, ...
            "Color", colorOrder(i_c,:), ...
            "Linestyle", resLineSpecs.ls_c, ...
            "Marker", resLineSpecs.m_c, ...
            "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

        if show_cov && isfield(res.("std"), variable_name)
            data_sig = f_trans(res.("std").(variable_name));
            if isfield(caseOpt,c) && isfield(caseOpt.(c),"cov")
                resLineSpecs = get_ls_specs(caseOpt.(c).cov);
            else
                resLineSpecs = get_ls_specs(struct);
                resLineSpecs.ls_c = "--";
            end

            leg(2,i_c) = plot(ax(i), t, data(i,:) + 3*data_sig(i,:), ...
                "DisplayName", sprintf("%s $3\\sigma$", label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

            plot(ax(i), t, data(i,:) - 3*data_sig(i,:), ...
                "DisplayName", label, ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        end
    end
    grid(ax(i),"on")

    title(ax(i),title_label)
    if isfield(plotOpt,"legend")
        leg_opt = plotOpt.legend;
    else
        leg_opt = {"Location","bestoutside","NumColumns",1, 'interpreter',"latex"};
    end
    legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
    xlabel(ax(i),"Time [s]")
    ylabel(ax(i),sprintf("%s %s", directions(i), y_unit))
end

end
function [f,ax] = plot_three_components_error(dataTot, t, cases_plot, caseOpt, plotOpt)

path_struct = plotOpt.path;
variable_name = plotOpt.variable_name;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;

directions = ["x","y","z"];

if isfield(plotOpt,"figure")
    f = plotOpt.figure; 
else
    f = figure();
end

if isfield(plotOpt,"axes")
    ax = plotOpt.axes; 
else
    ax = zeros(3,1);
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;
    end
end

if isfield(plotOpt,"show_cov")
    show_cov = plotOpt.show_cov; 
else
    show_cov = true;
end

if isfield(plotOpt,"show_mean")
    show_mean = plotOpt.show_mean; 
else
    show_mean = true;
end

for i = 1:3
    leg = nan(2,size(cases_plot,1));
    for i_c = 1:size(cases_plot,1)
        case_path = cases_plot(i_c, :);
        c = case_path(end);
        path_struct_i = [case_path, path_struct];
        path_struct_i_label = [case_path, path_struct(1:end-1),"label"];
        try 
            path_struct_i_struct = num2cell(path_struct_i);
            res = getfield(dataTot, path_struct_i_struct{:});
        catch 
            warning("%s is not in %s",c, join([path_struct{:}],"/"))
            continue
        end    
        
        try
            path_struct_i_label_struct = num2cell(path_struct_i_label);
            label_2 = getfield(dataTot, path_struct_i_label_struct{:});
        catch 
            warning("label %s is not in %s",c, join([path_struct{:}],"/"))
            label_2 = "";
        end    
        try
            label_1 = dataTot.(case_path(1)).("label");
        catch
            label_1 = "";
        end
        label = label_2 + newline + label_1;
        

        if isfield(caseOpt,c) && isfield(caseOpt.(c),"mean")
            resLineSpecs = get_ls_specs(caseOpt.(c).mean);
        else
            resLineSpecs = get_ls_specs(struct);
        end
        
        if show_mean && isfield(res, "mean") && isfield(res.("mean"), variable_name)
            data = f_trans(res.("mean").(variable_name));
            leg(1,i_c) = plot(ax(i), t, data(i,:), ...
                "DisplayName", label, ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        end        

        if show_cov && isfield(res.("std"), variable_name)
            data_sig = f_trans(res.("std").(variable_name));
            if isfield(caseOpt,c) && isfield(caseOpt.(c),"cov")
                resLineSpecs = get_ls_specs(caseOpt.(c).cov);
            else
                resLineSpecs = get_ls_specs(struct);
                resLineSpecs.ls_c = "--";
            end

            leg(2,i_c) = plot(ax(i), t, 3*data_sig(i,:), ...
                "DisplayName", sprintf("%s $3\\sigma$", label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

            plot(ax(i), t, -3*data_sig(i,:), ...
                "DisplayName", sprintf("%s %3\\sigma$", label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        end
    end
    grid(ax(i),"on")

    title(ax(i), title_label)
    if isfield(plotOpt,"legend")
        leg_opt = plotOpt.legend;
    else
        leg_opt = {"Location","bestoutside","NumColumns",1, 'interpreter',"latex"};
    end
    legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
    xlabel(ax(i), "Time [s]")
    ylabel(ax(i), sprintf("%s %s", directions(i), y_unit))
end

end
function [f,ax] = plot_three_components_error_log(dataTot, t, cases_plot, caseOpt, plotOpt)

path_struct = plotOpt.path;
variable_name = plotOpt.variable_name;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;

directions = ["x","y","z"];

if isfield(plotOpt,"figure")
    f = plotOpt.figure; 
else
    f = figure();
end

if isfield(plotOpt,"axes")
    ax = plotOpt.axes; 
else
    ax = zeros(3,1);
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;
    end
end

if isfield(plotOpt,"show_cov")
    show_cov = plotOpt.show_cov; 
else
    show_cov = true;
end

if isfield(plotOpt,"show_mean")
    show_mean = plotOpt.show_mean; 
else
    show_mean = true;
end

for i = 1:3
    leg = nan(2,size(cases_plot,1));
    for i_c = 1:size(cases_plot,1)
        case_path = cases_plot(i_c, :);
        c = case_path(end);
        path_struct_i = [case_path, path_struct];
        path_struct_i_label = [case_path, path_struct(1:end-1),"label"];
        try 
            path_struct_i_struct = num2cell(path_struct_i);
            res = getfield(dataTot, path_struct_i_struct{:});
        catch 
            warning("%s is not in %s",c, join([path_struct{:}],"/"))
            continue
        end    
        
        try
            path_struct_i_label_struct = num2cell(path_struct_i_label);
            label_2 = getfield(dataTot, path_struct_i_label_struct{:});
        catch 
            warning("label %s is not in %s",c, join([path_struct{:}],"/"))
            label_2 = "";
        end    
        try
            label_1 = dataTot.(case_path(1)).("label");
        catch
            label_1 = "";
        end
        label = label_2 + newline + label_1;
        

        if isfield(caseOpt,c) && isfield(caseOpt.(c),"mean")
            resLineSpecs = get_ls_specs(caseOpt.(c).mean);
        else
            resLineSpecs = get_ls_specs(struct);
        end
        
        if show_mean && isfield(res, "mean") && isfield(res.("mean"), variable_name)
            data = abs(f_trans(res.("mean").(variable_name)));
            leg(1,i_c) = plot(ax(i), t, data(i,:), ...
                "DisplayName", label, ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        end        

        if show_cov && isfield(res.("std"), variable_name)
            data_sig = f_trans(res.("std").(variable_name));
            if isfield(caseOpt,c) && isfield(caseOpt.(c),"cov")
                resLineSpecs = get_ls_specs(caseOpt.(c).cov);
            else
                resLineSpecs = get_ls_specs(struct);
                resLineSpecs.ls_c = "--";
            end

            leg(2,i_c) = plot(ax(i), t, 3*data_sig(i,:), ...
                "DisplayName", sprintf("%s $3\\sigma$", label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);


        end
    end
    grid(ax(i),"on")
    set(ax(i), "YScale", "log")
    title(ax(i), title_label)
    if isfield(plotOpt,"legend")
        leg_opt = plotOpt.legend;
    else
        leg_opt = {"Location","bestoutside","NumColumns",1, 'interpreter',"latex"};
    end
    legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
    xlabel(ax(i), "Time [s]")
    ylabel(ax(i), sprintf("%s %s", directions(i), y_unit))
end

end
function [f,ax] = plot_three_components_error_v2(dataTot, t, caseOpt, plotOpt)

pred_type = plotOpt.pred_type;
v_mean = plotOpt.variable_name_mean;
v_var = plotOpt.variable_name_var;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;

directions = ["x","y","z"];

if isfield(plotOpt,"figure")
    f = plotOpt.figure; 
else
    f = figure();
end

if isfield(plotOpt,"axes")
    ax = plotOpt.axes; 
else
    ax = zeros(3,1);
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;
    end
end

if isfield(plotOpt,"show_cov")
    show_cov = plotOpt.show_cov; 
else
    show_cov = true;
end
for i = 1:3
    leg = nan(2,numel(dataTot));
    for i_c = 1:numel(dataTot)
        res = dataTot(i_c);
        res = res{1};
        caseOpt_i = caseOpt(i_c);
        if ~isfield(res,pred_type)
            warning("%s not in %s",pred_type, c)
            continue
        end
        res_pred = res.(pred_type);

        if isstruct(caseOpt_i) && isfield(caseOpt_i,"mean")
            resLineSpecs = get_ls_specs(caseOpt_i.mean);
        else
            resLineSpecs = get_ls_specs(struct);
        end
        if isfield(res,"err") && isfield(res.err,v_mean)
            data = f_trans(res.err.(v_mean));
            leg(1,i_c) = plot(ax(i), t, data(i,:), ...
                "DisplayName", res.label, ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        end

        if show_cov && isfield(res_pred, v_var)
            data_sig = f_trans(res_pred.(v_var));
            if isstruct(caseOpt_i) && isfield(caseOpt_i,"cov")
                resLineSpecs = get_ls_specs(caseOpt_i.cov);
            else
                resLineSpecs = get_ls_specs(struct);
                resLineSpecs.ls_c = "--";
            end

            leg(2,i_c) = plot(ax(i), t, 3*data_sig(i,:), ...
                "DisplayName", sprintf("%s $3\\sigma$", res.label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

            plot(ax(i), t, -3*data_sig(i,:), ...
                "DisplayName", sprintf("%s %3\\sigma$", res.label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        end
    end
    grid(ax(i),"on")

    title(ax(i), title_label)
    if isfield(plotOpt,"legend")
        leg_opt = plotOpt.legend;
    else
        leg_opt = {"Location","bestoutside","NumColumns",1, 'interpreter',"latex"};
    end
    legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
    xlabel(ax(i), "Time [s]")
    ylabel(ax(i), sprintf("%s %s", directions(i), y_unit))
end

end
function [f,ax] = plot_three_components_filt_pred(dataTot, t, cases_plot, caseOpt, plotOpt)

pred_types = ["pred", "filt"];
v_mean = plotOpt.variable_name_mean;
v_var = plotOpt.variable_name_var;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;
get_color = @(i) colorOrder(mod1(i,size(colorOrder,1)), :);

if isfield(plotOpt,"components")
    directions = plotOpt.components;
else
    directions = ["x","y","z"];
end

if isfield(plotOpt,"show_cov")
    show_cov = plotOpt.show_cov; 
else
    show_cov = true;
end

if isfield(plotOpt,"figure")
    f = plotOpt.figure; 
else
    f = figure();
end

if isfield(plotOpt,"axes")
    ax = plotOpt.axes; 
else
    ax = zeros(3,1);
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;
    end
end

N_c = length(cases_plot);
N_p = length(pred_types);
for i = 1:3

    leg = nan(2,2,length(cases_plot));
    for i_c = 1:N_c
        c = cases_plot(i_c);
        res = dataTot.(c);
        for i_p = 1:N_p
            ii = sub2ind([N_p, N_c], i_p, i_c);
            pred_type = pred_types(i_p);
            
            if ~isfield(dataTot.(c),pred_type)
                continue
            end
                
            res_pred = dataTot.(c).(pred_type);
            
            if isfield(caseOpt,c) && isfield(caseOpt,pred_type) && isfield(caseOpt.(pred_type).(c),"mean")
                resLineSpecs = get_ls_specs(caseOpt.(pred_type).(c).mean);
            else
                resLineSpecs = get_ls_specs(struct);
            end
            if isfield(res_pred, v_mean)
                data = f_trans(res_pred.(v_mean));
                leg(1, i_p, i_c) = plot(ax(i), t, data(i,:), ...
                    "DisplayName", sprintf("%s %s", res.label, pred_type), ...
                    "Color", get_color(ii), ...
                    "Linestyle", resLineSpecs.ls_c, ...
                    "Marker", resLineSpecs.m_c, ...
                    "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
            end
            
            if show_cov && isfield(res_pred, v_var)
                data_sig = f_trans(res_pred.(v_var));
                if isfield(caseOpt,c) && isfield(caseOpt,pred_type) && isfield(caseOpt.(pred_type).(c),"cov")
                    resLineSpecs = get_ls_specs(caseOpt.(pred_type).(c).cov);
                else
                    resLineSpecs = get_ls_specs(struct);
                    resLineSpecs.ls_c = "--";
                end
                
                leg(2,i_p,i_c) = plot(ax(i), t, data(i,:) + 3*data_sig(i,:), ...
                    "DisplayName", sprintf("%s %s $3\\sigma$", res.label, pred_type), ...
                    "Color", get_color(ii), ...
                    "Linestyle", resLineSpecs.ls_c, ...
                    "Marker", resLineSpecs.m_c, ...
                    "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
                
                plot(ax(i), t, data(i,:) - 3*data_sig(i,:), ...
                    "DisplayName", res.label, ...
                    "Color", get_color(ii), ...
                    "Linestyle", resLineSpecs.ls_c, ...
                    "Marker", resLineSpecs.m_c, ...
                    "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
            end
        end
    end
    grid(ax(i),"on")

    title(ax(i),title_label)
    if isfield(plotOpt,"legend")
        leg_opt = plotOpt.legend;
    else
        leg_opt = {"Location","bestoutside","NumColumns",1, 'interpreter',"latex"};
    end
    legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
    xlabel(ax(i),"Time [s]")
    ylabel(ax(i),sprintf("%s %s", directions(i), y_unit))
end

end
function [f,ax] = plot_three_components_norm_error_log(dataTot, t, cases_plot, caseOpt, plotOpt)

pred_type = plotOpt.pred_type;
v_mean = plotOpt.variable_name_mean;
v_var = plotOpt.variable_name_var;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;


f = figure();
ax = gca();
hold on;

leg = nan(2,length(cases_plot));
for i_c = 1:length(cases_plot)
    c = cases_plot(i_c);
    res = dataTot.(c);
    res_pred = dataTot.(c).(pred_type);

    if isfield(caseOpt,c) && isfield(caseOpt.(c),"mean")
        resLineSpecs = get_ls_specs(caseOpt.(c).mean);
    else
        resLineSpecs = get_ls_specs(struct);
    end
    if ~isfield(res,"err")
        warning("%s not in %s","err", c)
        continue
    end
    if ~isfield(res.err,v_mean)
        warning("%s not in %s",v_mean, c)
        continue
    end
    data = f_trans(res.err.(v_mean));
    leg(1,i_c) = plot(t, norm_time(data), ...
        "DisplayName", res.label, ...
        "Color", colorOrder(i_c,:), ...
        "Linestyle", resLineSpecs.ls_c, ...
        "Marker", resLineSpecs.m_c, ...
        "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

    if isfield(res_pred, v_var)
        data_sig = f_trans(res_pred.(v_var));
        if isfield(caseOpt,c) && isfield(caseOpt.(c),"cov")
            resLineSpecs = get_ls_specs(caseOpt.(c).cov);
        else
            resLineSpecs = get_ls_specs(struct);
            resLineSpecs.ls_c = "--";
        end

        leg(2,i_c) = plot(t, 3*norm_time(data_sig), ...
            "DisplayName", sprintf("%s%s%s", res.label, newline, "$3 \sqrt{\sum \sigma_i }$"), ...
            "Color", colorOrder(i_c,:), ...
            "Linestyle", resLineSpecs.ls_c, ...
            "Marker", resLineSpecs.m_c, ...
            "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
    end
end
grid on
set(ax, "YScale","log")
title(title_label)
if isfield(plotOpt,"legend")
    leg_opt = plotOpt.legend;
else
    leg_opt = {"Location","bestoutside","NumColumns",1, 'interpreter',"latex"};
end
legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
xlabel("Time [s]")
ylabel(sprintf("%s", y_unit))


end
function [ax] = plot_three_components_update(ax, dataTot, t, cases_plot, caseOpt, plotOpt)

pred_type = plotOpt.pred_type;
v_mean = plotOpt.variable_name_mean;
v_var = plotOpt.variable_name_var;
f_trans = plotOpt.transformation;
title_label = plotOpt.title_label;
y_unit = plotOpt.y_unit;

N_end = length(t);
colorOrder = colororder;

directions = ["x","y","z"];

for i = 1:3
    hold(ax(i),"on");

    leg = zeros(2,length(cases_plot));
    for i_c = 1:length(cases_plot)
        c = cases_plot(i_c);
        res = dataTot.(c);
        res_pred = dataTot.(c).(pred_type);

        if isfield(caseOpt,c) && isfield(caseOpt.(c),"mean")
            resLineSpecs = get_ls_specs(caseOpt.(c).mean);
        else
            resLineSpecs = get_ls_specs(struct);
        end
        data = f_trans(res_pred.(v_mean));
        leg(1,i_c) = plot(ax(i),t, data(i,:), ...
            "DisplayName", res.label, ...
            "Color", colorOrder(i_c,:), ...
            "Linestyle", resLineSpecs.ls_c, ...
            "Marker", resLineSpecs.m_c, ...
            "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

        if isfield(res_pred, v_var)
            data_sig = f_trans(res_pred.(v_var));
            if isfield(caseOpt,c) && isfield(caseOpt.(c),"cov")
                resLineSpecs = get_ls_specs(caseOpt.(c).cov);
            else
                resLineSpecs = get_ls_specs(struct);
            end

            leg(2,i_c) = plot(ax(i), t, data(i,:) + 3*data_sig(i,:), ...
                "DisplayName", sprintf("%s $3\\sigma$", res.label), ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);

            plot(ax(i), t, data(i,:) - 3*data_sig(i,:), ...
                "DisplayName", res.label, ...
                "Color", colorOrder(i_c,:), ...
                "Linestyle", resLineSpecs.ls_c, ...
                "Marker", resLineSpecs.m_c, ...
                "MarkerIndices", resLineSpecs.m_offset:resLineSpecs.m_space:N_end);
        else
            leg(2,i_c) = nan;
        end
    end
    grid on

    title(ax(i), title_label)
    if isfield(plotOpt,"legend")
        leg_opt = plotOpt.legend;
    else
        leg_opt = {"Location","bestoutside","NumColumns",1, 'interpreter',"latex"};
    end
    legend(reshape(leg(~isnan(leg)),[],1), leg_opt{:})
    xlabel(ax(i), "Time [s]")
    ylabel(ax(i), sprintf("%s %s", directions(i), y_unit))
end

end
