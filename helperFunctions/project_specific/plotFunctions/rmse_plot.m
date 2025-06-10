function [f, ax] = IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)

v_name = figurePlotOpts.variable_name;
f_trans = figurePlotOpts.transformation;
title_label = figurePlotOpts.title_label;
y_unit = figurePlotOpts.y_unit;

f = figure();
colorOrder = colororder;
ax = zeros(3,1);
dirs = ["x", "y","z"];
for i = 1:3
    ax(i) = subplot(3,1,i);
    hold on
    leg = nan(size(cases_plot,1),1);
    for i_c = 1:size(cases_plot,1)
        case_path = cases_plot(i_c, :);
        c = case_path(end);
        path_struct_i = case_path;
        path_struct_i_label = [case_path(1), "label"];
        try 
            path_struct_i_struct = num2cell(path_struct_i);
            res = getfield(data, path_struct_i_struct{:});
            
            path_struct_i_label_struct = num2cell(path_struct_i_label);
            label_1 = getfield(data, path_struct_i_label_struct{:});
            label_2 = res.("label");
            label = label_2 + newline + label_1;
        catch ME
            warning("%s is not in %s",c, join([path_struct_i_struct{:}],"/"))
            continue
        end

        if isfield(casePlotOpts,c)
            resLineSpecs = get_ls_specs(casePlotOpts.(c));
        else
            resLineSpecs = get_ls_specs(struct);
        end
        resData = f_trans(res.(v_name));
        temp = plot(IN_time_array, squeeze(resData(:,i,:)), ...
            "Color", colorOrder(i_c,:),...
            "LineStyle",resLineSpecs.ls_c,....
            "DisplayName", label);
        leg(i_c) = temp(1);
    end
    grid on
    legend(reshape(leg(~isnan(leg)),[],1),"Location","best", "Interpreter","latex")
    ylabel(sprintf("%s %s", dirs(i),y_unit))
    xlabel("Time [s]")
    title(title_label)

end

end
function [f, ax] = IN_error_position(data, IN_time_array, cases_plot, casePlotOpts)
    figurePlotOpts = struct;
    figurePlotOpts.variable_name = "p";
    figurePlotOpts.transformation = @(x) x;
    figurePlotOpts.title_label = "Position Error";
    figurePlotOpts.y_unit = "[m]";

    [f,ax] = IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);

end
function [f, ax] = IN_error_rotation(data, IN_time_array, cases_plot, casePlotOpts)

    figurePlotOpts = struct;
    figurePlotOpts.variable_name = "R";
    figurePlotOpts.transformation = @rad2deg;
    figurePlotOpts.title_label = "Rotation Error";
    figurePlotOpts.y_unit = "[deg]";

    [f,ax] = IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);

end
function [f,a] = plot_statistics_position(data, IN_time_array, cases_plot, casePlotOpts, plotType, extras)

    if nargin < 5
        plotType = "error";
    end
    if nargin < 6
        extras = struct;
    end
    figurePlotOpts = extras;
    figurePlotOpts.variable_name = "p";
    figurePlotOpts.transformation = @(x) x;    
    figurePlotOpts.y_unit = "[m]";
    
    if isfield(extras, "title")
        figurePlotOpts.title_label = sprintf("%s ", extras.title);
    else
        figurePlotOpts.title_label = "";
    end
    
    if strcmp(plotType, "error")
        figurePlotOpts.title_label = figurePlotOpts.title_label + "Position Error";
        [f,a] = IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);
    elseif strcmp(plotType, "rmse-all")
        figurePlotOpts.title_label = figurePlotOpts.title_label + "RMSE Position";
        [f,a] = rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);
    elseif strcmp(plotType, "rmse-comp")
        figurePlotOpts.title_label = figurePlotOpts.title_label + "RMSE component-wise Position";
        [f,a] = rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);
    else
        error("Invalid plot type %s", plotType)        
    end
end
function [f,a] = plot_statistics_position_v2(data, IN_time_array, casePlotOpts, plotType, extras)

    if nargin < 4
        plotType = "error";
    end
    if nargin < 5
        extras = struct;
    end
    figurePlotOpts = extras;
    figurePlotOpts.variable_name = "p";
    figurePlotOpts.transformation = @(x) x;    
    figurePlotOpts.y_unit = "[m]";
    
    if isfield(extras, "title")
        figurePlotOpts.title_label = sprintf("%s ", extras.title);
    else
        figurePlotOpts.title_label = "";
    end
    
    if strcmp(plotType, "rmse-all")
        figurePlotOpts.title_label = figurePlotOpts.title_label + "RMSE Position";
        [f,a] = rmse_all_components_general_v2(data, IN_time_array, casePlotOpts, figurePlotOpts);
    else
        error("Invalid plot type %s", plotType)        
    end
end
function [f,a] = plot_statistics_rotation(data, IN_time_array, cases_plot, casePlotOpts, plotType)
   
    if nargin < 5
        plotType = "error";
    end
    if nargin < 6
        extras = struct;
    end
    figurePlotOpts = extras;
    figurePlotOpts.variable_name = "R";
    figurePlotOpts.transformation = @rad2deg;    
    figurePlotOpts.y_unit = "[deg]";
    
    if isfield(extras, "title")
        figurePlotOpts.title_label = sprintf("%s ", extras.title);
    else
        figurePlotOpts.title_label = "";
    end
    
    if strcmp(plotType, "error")
        figurePlotOpts.title_label = figurePlotOpts.title_label + "Position Error";
        [f,a] = IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);
    elseif strcmp(plotType, "rmse-all")
        figurePlotOpts.title_label = figurePlotOpts.title_label + "RMSE Position";
        [f,a] = rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);
    elseif strcmp(plotType, "rmse-comp")
        figurePlotOpts.title_label = figurePlotOpts.title_label + "RMSE component-wise Rotation";
        [f,a] = rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);
    else
        error("Invalid plot type %s", plotType)        
    end
end
function [f,a] = rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)

v_name = figurePlotOpts.variable_name;
f_trans = figurePlotOpts.transformation;
title_label = figurePlotOpts.title_label;
y_unit = figurePlotOpts.y_unit;

f = figure();
a = gca;
colorOrder = colororder;
Nc_max = size(colorOrder,1);
N_end = length(IN_time_array);

hold on
for i_c = 1:size(cases_plot,1)
    case_path = cases_plot(i_c, :);
    c = case_path(end);
    path_struct_i = case_path;
    path_struct_i_label = [case_path(1), "label"];
    try
        path_struct_i_struct = num2cell(path_struct_i);
        res = getfield(data, path_struct_i_struct{:});
        
        path_struct_i_label_struct = num2cell(path_struct_i_label);
        label_1 = getfield(data, path_struct_i_label_struct{:});
        label_2 = res.("label");
        label = label_2 + newline + label_1;
    catch ME
        warning("%s is not in %s",c, join([path_struct_i_struct{:}],"/"))
        continue
    end
    
    
    resData = f_trans(res.(v_name));
    res1 = reshape(resData,size(resData,1),[]);

    if isfield(casePlotOpts,c)
        lineSpecs = get_ls_specs(casePlotOpts.(c));
    else
        lineSpecs = get_ls_specs(struct);
    end

    plot(IN_time_array, sqrt(mean(res1.^2,2)), ...
        "Color", colorOrder(mod1(i_c,Nc_max),:),...
        "LineStyle", lineSpecs.ls_c,....
        "DisplayName", label,...
        "Marker", lineSpecs.m_c, ...
        "MarkerIndices",lineSpecs.m_offset:lineSpecs.m_space:N_end);

end
if isfield(figurePlotOpts,"legend")
    leg_opt = figurePlotOpts.legend;
else
    leg_opt = {"Location","best","NumColumns",1, 'Interpreter',"latex"};
end
legend(leg_opt{:})
grid on
ylabel(sprintf("%s",y_unit))
xlabel("Time [s]")
if ~(isfield(figurePlotOpts,"omit_title") && figurePlotOpts.omit_title)
    title(title_label)
end


end
function [f,a] = rmse_all_components_general_v2(data, IN_time_array, casePlotOptsTot, figurePlotOpts)

v_name = figurePlotOpts.variable_name;
f_trans = figurePlotOpts.transformation;
title_label = figurePlotOpts.title_label;
y_unit = figurePlotOpts.y_unit;

f = figure();
a = gca;
colorOrder = colororder;
Nc_max = size(colorOrder,1);
N_end = length(IN_time_array);

hold on
for i_c = 1:numel(data)
    res = data(i_c);
    casePlotOpts = casePlotOptsTot(i_c);
    resData = f_trans(res.(v_name));
    res1 = reshape(resData,size(resData,1),[]);

    if isstruct(casePlotOpts)
        lineSpecs = get_ls_specs(casePlotOpts);
    else
        lineSpecs = get_ls_specs(struct);
    end

    plot(IN_time_array, sqrt(mean(res1.^2,2)), ...
        "Color", colorOrder(mod1(i_c,Nc_max),:),...
        "LineStyle", lineSpecs.ls_c,....
        "DisplayName", res.label,...
        "Marker", lineSpecs.m_c, ...
        "MarkerIndices",lineSpecs.m_offset:lineSpecs.m_space:N_end);

end
if isfield(figurePlotOpts,"legend")
    leg_opt = figurePlotOpts.legend;
else
    leg_opt = {"Location","best","NumColumns",1, 'Interpreter',"latex"};
end
legend(leg_opt{:})
grid on
ylabel(sprintf("%s",y_unit))
xlabel("Time [s]")
if ~(isfield(figurePlotOpts,"omit_title") && figurePlotOpts.omit_title)
    title(title_label)
end


end
function [f,a] = rmse_all_components_position(data, IN_time_array, cases_plot, casePlotOpts)

    figurePlotOpts = struct;
    figurePlotOpts.variable_name = "p";
    figurePlotOpts.transformation = @(x) x;
    figurePlotOpts.title_label = "RMSE Position";
    figurePlotOpts.y_unit = "[m]";

    [f,a] = rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);

end
function [f,a] = rmse_all_components_rotation(data, IN_time_array, cases_plot, casePlotOpts)

    figurePlotOpts = struct;
    figurePlotOpts.variable_name = "R";
    figurePlotOpts.transformation = @rad2deg;
    figurePlotOpts.title_label = "RMSE Rotation";
    figurePlotOpts.y_unit = "[deg]";

    [f,a] = rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);

end
function [f,ax] = rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)

v_name = figurePlotOpts.variable_name;
f_trans = figurePlotOpts.transformation;
title_label = figurePlotOpts.title_label;
y_unit = figurePlotOpts.y_unit;

f = figure();
dirs = ["x", "y", "z"];
colorOrder = colororder;
N_end = length(IN_time_array);
ax = zeros(3,1);
for i = 1:3
    ax(i) = subplot(3,1,i);
    hold on
    for i_c = 1:length(cases_plot)
        c = cases_plot(i_c);
        res = data.(c);

        if isfield(casePlotOpts,c)
            lineSpecs = get_ls_specs(casePlotOpts.(c));
        else
            lineSpecs = get_ls_specs(struct);
        end
        resData = f_trans(res.(v_name));
        plot(IN_time_array, sqrt(mean(squeeze(resData(:,i,:)).^2,2)), ...
            "Color", colorOrder(i_c,:),...
            "LineStyle", lineSpecs.ls_c,...
            "DisplayName", res.label,...
            "Marker", lineSpecs.m_c, ...
            "MarkerIndices",lineSpecs.m_offset:lineSpecs.m_space:N_end);
    end
    legend("Location", "best", "Interpreter","latex")
    grid on
    ylabel(sprintf("%s %s", dirs(i),y_unit))
    xlabel("Time [s]")
    title(title_label)

end

end
function [f,ax] = rmse_components_position(data, IN_time_array, cases_plot, casePlotOpts)

    figurePlotOpts = struct;
    figurePlotOpts.variable_name = "p";
    figurePlotOpts.transformation = @(x) x;
    figurePlotOpts.title_label = "RMSE component-wise Position";
    figurePlotOpts.y_unit = "[m]";

    [f,ax] = rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);

end
function [f,ax] = rmse_components_rotation(data, IN_time_array, cases_plot, casePlotOpts)

    figurePlotOpts = struct;
    figurePlotOpts.variable_name = "R";
    figurePlotOpts.transformation = @rad2deg;
    figurePlotOpts.title_label = "RMSE component-wise Rotation";
    figurePlotOpts.y_unit = "[deg]";


    [f,ax] = rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts);

end
