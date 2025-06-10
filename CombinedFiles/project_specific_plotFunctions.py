import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.spatial.transform import Rotation as R
import warnings


def mod1(x, m):
    """MATLAB-style mod1 function (1-indexed modulo)"""
    return ((x - 1) % m) + 1


def my_rotm2eul(rotm):
    """Convert rotation matrix to Euler angles (roll, pitch, yaw)"""
    # Assuming rotm is a 3x3xN array where N is the number of time steps
    if rotm.ndim == 3:
        euler_angles = np.zeros((3, rotm.shape[2]))
        for i in range(rotm.shape[2]):
            r = R.from_matrix(rotm[:, :, i])
            euler_angles[:, i] = r.as_euler('xyz')  # roll, pitch, yaw
        return euler_angles
    else:
        r = R.from_matrix(rotm)
        return r.as_euler('xyz')


def norm_time(data):
    """Calculate norm along first axis (equivalent to MATLAB's norm along columns)"""
    return np.linalg.norm(data, axis=0)


def get_ls_specs(opt):
    """
    Get line specifications from options structure
    """
    res = {}
    
    if "ls" in opt:
        res["ls_c"] = opt["ls"]
    else:
        res["ls_c"] = "-"
    
    if "m" in opt:
        res["m_c"] = opt["m"]
    else:
        res["m_c"] = "none"
    
    if "m_space" in opt:
        res["m_space"] = opt["m_space"]
    else:
        res["m_space"] = 1
    
    if "m_offset" in opt:
        res["m_offset"] = opt["m_offset"]
    else:
        res["m_offset"] = 1
    
    return res


def IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts):
    """
    General function for plotting IN errors
    """
    v_name = figurePlotOpts["variable_name"]
    f_trans = figurePlotOpts["transformation"]
    title_label = figurePlotOpts["title_label"]
    y_unit = figurePlotOpts["y_unit"]
    
    f, ax = plt.subplots(3, 1, figsize=(10, 12))
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    dirs = ["x", "y", "z"]
    
    for i in range(3):
        ax[i].hold = True  # Enable hold
        leg = [None] * cases_plot.shape[0]
        
        for i_c in range(cases_plot.shape[0]):
            case_path = cases_plot[i_c, :]
            c = case_path[-1]
            path_struct_i = case_path
            path_struct_i_label = np.concatenate([case_path[:-1], ["label"]])
            
            try:
                # Navigate through nested dictionary structure
                res = data
                for key in path_struct_i:
                    res = res[key]
                
                # Get label
                label_res = data
                for key in path_struct_i_label:
                    label_res = label_res[key]
                label_1 = label_res
                label_2 = res["label"]
                label = f"{label_2}\n{label_1}"
                
            except KeyError:
                warnings.warn(f"{c} is not in {'/'.join(path_struct_i)}")
                continue
            
            if c in casePlotOpts:
                resLineSpecs = get_ls_specs(casePlotOpts[c])
            else:
                resLineSpecs = get_ls_specs({})
            
            resData = f_trans(res[v_name])
            temp = ax[i].plot(IN_time_array, resData[i, :], 
                            color=color_order[i_c],
                            linestyle=resLineSpecs["ls_c"],
                            label=label)
            leg[i_c] = temp[0]
        
        ax[i].grid(True)
        # Filter out None values for legend
        valid_legs = [l for l in leg if l is not None]
        if valid_legs:
            ax[i].legend(valid_legs, loc="best")
        ax[i].set_ylabel(f"{dirs[i]} {y_unit}")
        ax[i].set_xlabel("Time [s]")
        ax[i].set_title(title_label)
    
    return f, ax


def IN_error_position(data, IN_time_array, cases_plot, casePlotOpts):
    """
    Plot position error
    """
    figurePlotOpts = {
        "variable_name": "p",
        "transformation": lambda x: x,
        "title_label": "Position Error",
        "y_unit": "[m]"
    }
    
    return IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)


def IN_error_rotation(data, IN_time_array, cases_plot, casePlotOpts):
    """
    Plot rotation error
    """
    figurePlotOpts = {
        "variable_name": "R",
        "transformation": np.rad2deg,
        "title_label": "Rotation Error",
        "y_unit": "[deg]"
    }
    
    return IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)

# ... existing code ...


def plot_statistics_position(data, IN_time_array, cases_plot, casePlotOpts, plotType="error", extras=None):
    """
    Plot position statistics with different plot types
    """
    if extras is None:
        extras = {}
    
    figurePlotOpts = extras.copy()
    figurePlotOpts["variable_name"] = "p"
    figurePlotOpts["transformation"] = lambda x: x
    figurePlotOpts["y_unit"] = "[m]"
    
    if "title" in extras:
        figurePlotOpts["title_label"] = f"{extras['title']} "
    else:
        figurePlotOpts["title_label"] = ""
    
    if plotType == "error":
        figurePlotOpts["title_label"] += "Position Error"
        return IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)
    elif plotType == "rmse-all":
        figurePlotOpts["title_label"] += "RMSE Position"
        return rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)
    elif plotType == "rmse-comp":
        figurePlotOpts["title_label"] += "RMSE component-wise Position"
        return rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)
    else:
        raise ValueError(f"Invalid plot type {plotType}")


def plot_statistics_position_v2(data, IN_time_array, casePlotOpts, plotType="error", extras=None):
    """
    Plot position statistics version 2
    """
    if extras is None:
        extras = {}
    
    figurePlotOpts = extras.copy()
    figurePlotOpts["variable_name"] = "p"
    figurePlotOpts["transformation"] = lambda x: x
    figurePlotOpts["y_unit"] = "[m]"
    
    if "title" in extras:
        figurePlotOpts["title_label"] = f"{extras['title']} "
    else:
        figurePlotOpts["title_label"] = ""
    
    if plotType == "rmse-all":
        figurePlotOpts["title_label"] += "RMSE Position"
        return rmse_all_components_general_v2(data, IN_time_array, casePlotOpts, figurePlotOpts)
    else:
        raise ValueError(f"Invalid plot type {plotType}")


def plot_statistics_rotation(data, IN_time_array, cases_plot, casePlotOpts, plotType="error", extras=None):
    """
    Plot rotation statistics
    """
    if extras is None:
        extras = {}
    
    figurePlotOpts = extras.copy()
    figurePlotOpts["variable_name"] = "R"
    figurePlotOpts["transformation"] = np.rad2deg
    figurePlotOpts["y_unit"] = "[deg]"
    
    if "title" in extras:
        figurePlotOpts["title_label"] = f"{extras['title']} "
    else:
        figurePlotOpts["title_label"] = ""
    
    if plotType == "error":
        figurePlotOpts["title_label"] += "Position Error"  # Note: This seems to be a bug in original MATLAB code
        return IN_error_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)
    elif plotType == "rmse-all":
        figurePlotOpts["title_label"] += "RMSE Position"  # Note: This seems to be a bug in original MATLAB code
        return rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)
    elif plotType == "rmse-comp":
        figurePlotOpts["title_label"] += "RMSE component-wise Rotation"
        return rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)
    else:
        raise ValueError(f"Invalid plot type {plotType}")


def rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts):
    """
    General RMSE plotting for all components
    """
    v_name = figurePlotOpts["variable_name"]
    f_trans = figurePlotOpts["transformation"]
    title_label = figurePlotOpts["title_label"]
    y_unit = figurePlotOpts["y_unit"]
    
    f, a = plt.subplots(figsize=(10, 6))
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    Nc_max = len(color_order)
    N_end = len(IN_time_array)
    
    for i_c in range(cases_plot.shape[0]):
        case_path = cases_plot[i_c, :]
        c = case_path[-1]
        path_struct_i = case_path
        path_struct_i_label = np.concatenate([case_path[:-1], ["label"]])
        
        try:
            # Navigate through nested dictionary structure
            res = data
            for key in path_struct_i:
                res = res[key]
            
            # Get label
            label_res = data
            for key in path_struct_i_label:
                label_res = label_res[key]
            label_1 = label_res
            label_2 = res["label"]
            label = f"{label_2}\n{label_1}"
            
        except KeyError:
            warnings.warn(f"{c} is not in {'/'.join(path_struct_i)}")
            continue
        
        resData = f_trans(res[v_name])
        res1 = resData.reshape(resData.shape[0], -1)
        
        if c in casePlotOpts:
            lineSpecs = get_ls_specs(casePlotOpts[c])
        else:
            lineSpecs = get_ls_specs({})
        
        # Create marker indices
        marker_indices = range(lineSpecs["m_offset"] - 1, N_end, lineSpecs["m_space"])
        
        a.plot(IN_time_array, np.sqrt(np.mean(res1**2, axis=1)),
               color=color_order[mod1(i_c + 1, Nc_max) - 1],
               linestyle=lineSpecs["ls_c"],
               label=label,
               marker=lineSpecs["m_c"],
               markevery=marker_indices)
    
    if "legend" in figurePlotOpts:
        leg_opt = figurePlotOpts["legend"]
        a.legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
    else:
        a.legend(loc="best", ncol=1)
    
    a.grid(True)
    a.set_ylabel(y_unit)
    a.set_xlabel("Time [s]")
    
    if not ("omit_title" in figurePlotOpts and figurePlotOpts["omit_title"]):
        a.set_title(title_label)
    
    return f, a

# ... existing code ...


def rmse_all_components_general_v2(data, IN_time_array, casePlotOptsTot, figurePlotOpts):
    """
    General RMSE plotting for all components version 2
    """
    v_name = figurePlotOpts["variable_name"]
    f_trans = figurePlotOpts["transformation"]
    title_label = figurePlotOpts["title_label"]
    y_unit = figurePlotOpts["y_unit"]
    
    f, a = plt.subplots(figsize=(10, 6))
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    Nc_max = len(color_order)
    N_end = len(IN_time_array)
    
    for i_c in range(len(data)):
        res = data[i_c]
        casePlotOpts = casePlotOptsTot[i_c]
        resData = f_trans(res[v_name])
        res1 = resData.reshape(resData.shape[0], -1)
        
        if isinstance(casePlotOpts, dict):
            lineSpecs = get_ls_specs(casePlotOpts)
        else:
            lineSpecs = get_ls_specs({})
        
        # Create marker indices
        marker_indices = range(lineSpecs["m_offset"] - 1, N_end, lineSpecs["m_space"])
        
        a.plot(IN_time_array, np.sqrt(np.mean(res1**2, axis=1)),
               color=color_order[mod1(i_c + 1, Nc_max) - 1],
               linestyle=lineSpecs["ls_c"],
               label=res["label"],
               marker=lineSpecs["m_c"],
               markevery=marker_indices)
    
    if "legend" in figurePlotOpts:
        leg_opt = figurePlotOpts["legend"]
        a.legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
    else:
        a.legend(loc="best", ncol=1)
    
    a.grid(True)
    a.set_ylabel(y_unit)
    a.set_xlabel("Time [s]")
    
    if not ("omit_title" in figurePlotOpts and figurePlotOpts["omit_title"]):
        a.set_title(title_label)
    
    return f, a


def rmse_all_components_position(data, IN_time_array, cases_plot, casePlotOpts):
    """
    RMSE for all position components
    """
    figurePlotOpts = {
        "variable_name": "p",
        "transformation": lambda x: x,
        "title_label": "RMSE Position",
        "y_unit": "[m]"
    }
    
    return rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)


def rmse_all_components_rotation(data, IN_time_array, cases_plot, casePlotOpts):
    """
    RMSE for all rotation components
    """
    figurePlotOpts = {
        "variable_name": "R",
        "transformation": np.rad2deg,
        "title_label": "RMSE Rotation",
        "y_unit": "[deg]"
    }
    
    return rmse_all_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)


def rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts):
    """
    General RMSE plotting for individual components
    """
    v_name = figurePlotOpts["variable_name"]
    f_trans = figurePlotOpts["transformation"]
    title_label = figurePlotOpts["title_label"]
    y_unit = figurePlotOpts["y_unit"]
    
    f, ax = plt.subplots(3, 1, figsize=(10, 12))
    dirs = ["x", "y", "z"]
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    N_end = len(IN_time_array)
    
    for i in range(3):
        for i_c in range(len(cases_plot)):
            c = cases_plot[i_c]
            res = data[c]
            
            if c in casePlotOpts:
                lineSpecs = get_ls_specs(casePlotOpts[c])
            else:
                lineSpecs = get_ls_specs({})
            
            resData = f_trans(res[v_name])
            
            # Create marker indices
            marker_indices = range(lineSpecs["m_offset"] - 1, N_end, lineSpecs["m_space"])
            
            ax[i].plot(IN_time_array, np.sqrt(np.mean(resData[i, :, :]**2, axis=1)),
                      color=color_order[i_c],
                      linestyle=lineSpecs["ls_c"],
                      label=res["label"],
                      marker=lineSpecs["m_c"],
                      markevery=marker_indices)
        
        ax[i].legend(loc="best")
        ax[i].grid(True)
        ax[i].set_ylabel(f"{dirs[i]} {y_unit}")
        ax[i].set_xlabel("Time [s]")
        ax[i].set_title(title_label)
    
    return f, ax


def rmse_components_position(data, IN_time_array, cases_plot, casePlotOpts):
    """
    RMSE for position components
    """
    figurePlotOpts = {
        "variable_name": "p",
        "transformation": lambda x: x,
        "title_label": "RMSE component-wise Position",
        "y_unit": "[m]"
    }
    
    return rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)


def rmse_components_rotation(data, IN_time_array, cases_plot, casePlotOpts):
    """
    RMSE for rotation components
    """
    figurePlotOpts = {
        "variable_name": "R",
        "transformation": np.rad2deg,
        "title_label": "RMSE component-wise Rotation",
        "y_unit": "[deg]"
    }
    
    return rmse_components_general(data, IN_time_array, cases_plot, casePlotOpts, figurePlotOpts)

# ... existing code ...


def plot_acc_pos(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot accelerometer position
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name_mean"] = "r"
    plotOpt["variable_name_var"] = "r_std"
    plotOpt["transformation"] = lambda x: x * 1e3
    plotOpt["title_label"] = "Accelerometer pos. 1"
    plotOpt["y_unit"] = "[mm]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_angular_acceleration(data, t, cases_plot, caseOpt, plotType="normal"):
    """
    Plot angular acceleration
    """
    plotOpt = {
        "variable_name_mean": "omega_dot",
        "variable_name_var": "omega_dot_std",
        "transformation": np.rad2deg,
        "title_label": "Angular Acceleration",
        "y_unit": "[deg/s^2]",
        "pred_type": "pred"
    }
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_angular_velocity(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot angular velocity
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name_mean"] = "w"
    plotOpt["variable_name_var"] = "w_std"
    plotOpt["transformation"] = np.rad2deg
    plotOpt["title_label"] = "Angular Velocity"
    plotOpt["y_unit"] = "[deg/s]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "filt_and_pred":
        return plot_three_components_filt_pred(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_angular_velocity_filt_and_pred(data, t, cases_plot, caseOpt):
    """
    Plot angular velocity with both filtered and predicted values
    """
    plotOpt = {
        "variable_name_mean": "w",
        "variable_name_var": "w_std",
        "transformation": np.rad2deg,
        "title_label": "Angular Velocity",
        "y_unit": "[deg/s]",
        "pred_type": "filt"
    }
    
    f, ax = plt.subplots(3, 1, figsize=(10, 12))
    for i in range(3):
        ax[i].grid(True)
    
    # Plot filtered data
    plot_three_components_update(ax, data, t, cases_plot, caseOpt, plotOpt)
    
    # Update for predicted data
    plotOpt["pred_type"] = "pred"
    plot_three_components_update(ax, data, t, cases_plot, caseOpt, plotOpt)
    
    return f, ax


def plot_bias_accelerometer(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot accelerometer bias
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name_mean"] = "b_a"
    plotOpt["variable_name_var"] = "b_a_std"
    plotOpt["transformation"] = lambda x: x
    plotOpt["title_label"] = "Bias accelerometer"
    plotOpt["y_unit"] = "[m/s^2]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_bias_gyroscopes(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot gyroscope bias
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name"] = "b_g"
    plotOpt["transformation"] = np.rad2deg
    
    if "title" in plotOpt:
        plotOpt["title_label"] = f"{plotOpt['title']}, Bias gyroscopes"
    else:
        plotOpt["title_label"] = "Bias gyroscopes"
    
    plotOpt["y_unit"] = "[deg/s]"
    
    if "path" not in plotOpt:
        plotOpt["path"] = ["filt"]
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_bias_omega_dot(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot omega dot bias
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name"] = "b_omega_dot"
    plotOpt["transformation"] = np.rad2deg
    
    if "title" in plotOpt:
        plotOpt["title_label"] = f"{plotOpt['title']}, Bias omega dot"
    else:
        plotOpt["title_label"] = "Bias omega dot"
    
    plotOpt["y_unit"] = "[deg/s^2]"
    
    if "path" not in plotOpt:
        plotOpt["path"] = ["filt"]
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_bias_s(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot specific force bias
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name"] = "b_s"
    plotOpt["transformation"] = lambda x: x
    
    if "title" in plotOpt:
        plotOpt["title_label"] = f"{plotOpt['title']}, Bias specific force"
    else:
        plotOpt["title_label"] = "Bias specific force"
    
    plotOpt["y_unit"] = "[m/s^2]"
    
    if "path" not in plotOpt:
        plotOpt["path"] = ["filt"]
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)
    
# ... existing code ...


def plot_navigation_acceleration(data, t, cases_plot, caseOpt, plotType="normal"):
    """
    Plot navigation acceleration
    """
    plotOpt = {
        "variable_name_mean": "v_dot",
        "variable_name_var": "v_dot_std",
        "transformation": lambda x: x,
        "title_label": "Navigation Acceleration",
        "y_unit": "[m/s^2]",
        "pred_type": "pred"
    }
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_navigation_position(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot navigation position
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name"] = "p"
    plotOpt["transformation"] = lambda x: x
    
    if "title" in plotOpt:
        plotOpt["title_label"] = f"{plotOpt['title']} "
    else:
        plotOpt["title_label"] = ""
    
    plotOpt["title_label"] += "Navigation position"
    plotOpt["y_unit"] = "[m]"
    
    if "path" not in plotOpt:
        plotOpt["path"] = ["filt"]
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_navigation_position_v2(data, t, caseOpt, plotType="normal", extras=None):
    """
    Plot navigation position version 2
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
        print("Using Extras")
    
    plotOpt["variable_name_mean"] = "p"
    plotOpt["variable_name_var"] = "p_std"
    plotOpt["transformation"] = lambda x: x
    
    if "title" in plotOpt:
        plotOpt["title_label"] = f"{plotOpt['title']} "
    else:
        plotOpt["title_label"] = ""
    
    plotOpt["title_label"] += "Navigation position"
    plotOpt["y_unit"] = "[m]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "error":
        return plot_three_components_error_v2(data, t, caseOpt, plotOpt)
    else:
        raise ValueError("wrong")


def plot_navigation_velocity(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot navigation velocity
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name_mean"] = "v"
    plotOpt["variable_name_var"] = "v_std"
    plotOpt["transformation"] = lambda x: x
    plotOpt["title_label"] = "Navigation Velocity"
    plotOpt["y_unit"] = "[m/s]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_nine_components(dataTot, t, cases_plot, caseOpt, plotOpt):
    """
    Plot nine components (3x3 matrix elements)
    """
    pred_type = plotOpt["pred_type"]
    v_name = plotOpt["variable_name_mean"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    
    if "components" in plotOpt:
        directions = plotOpt["components"]
    else:
        directions = ["(1,1)", "(1,2)", "(1,3)",
                      "(2,1)", "(2,2)", "(2,3)",
                      "(3,1)", "(3,2)", "(3,3)"]
    
    if "show_cov" in plotOpt:
        show_cov = plotOpt["show_cov"]
    else:
        show_cov = True
    
    if "figure" in plotOpt:
        f = plotOpt["figure"]
    else:
        f = plt.figure(figsize=(15, 12))
    
    if "axes" in plotOpt:
        ax = plotOpt["axes"]
    else:
        ax = []
        for i in range(9):
            ax.append(plt.subplot(3, 3, i + 1))
    
    for i in range(9):
        leg = []
        for i_c in range(len(cases_plot)):
            c = cases_plot[i_c]
            res = dataTot[c]
            res_pred = dataTot[c][pred_type]
            
            if c in caseOpt and "mean" in caseOpt[c]:
                resLineSpecs = get_ls_specs(caseOpt[c]["mean"])
            else:
                resLineSpecs = get_ls_specs({})
            
            if v_name not in res_pred["mean"]:
                raise KeyError(f"{v_name} is not in {c}")
            
            data = f_trans(res_pred["mean"][v_name])
            
            # Create marker indices
            marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
            
            line1 = ax[i].plot(t, data[i, :],
                              label=res["label"],
                              color=color_order[i_c],
                              linestyle=resLineSpecs["ls_c"],
                              marker=resLineSpecs["m_c"],
                              markevery=marker_indices)
            leg.extend(line1)
            
            if show_cov and v_name in res_pred["std"]:
                data_sig = f_trans(res_pred["std"][v_name])
                
                if c in caseOpt and "cov" in caseOpt[c]:
                    resLineSpecs = get_ls_specs(caseOpt[c]["cov"])
                else:
                    resLineSpecs = get_ls_specs({})
                    resLineSpecs["ls_c"] = "--"
                
                line2 = ax[i].plot(t, data[i, :] + 3 * data_sig[i, :],
                                  label=f"{res['label']} $3\\sigma$",
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line2)
                
                ax[i].plot(t, data[i, :] - 3 * data_sig[i, :],
                          color=color_order[i_c],
                          linestyle=resLineSpecs["ls_c"],
                          marker=resLineSpecs["m_c"],
                          markevery=marker_indices)
        
        ax[i].grid(True)
        ax[i].set_title(title_label)
        
        if "legend" in plotOpt:
            leg_opt = plotOpt["legend"]
            ax[i].legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
        else:
            ax[i].legend(loc="best", ncol=1)
        
        ax[i].set_xlabel("Time [s]")
        ax[i].set_ylabel(f"{directions[i]} {y_unit}")
    
    return f, ax


def plot_rotation(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot rotation (only works for error types)
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name"] = "R"
    plotOpt["transformation"] = np.rad2deg
    plotOpt["y_unit"] = "[deg]"
    
    if "path" not in plotOpt:
        plotOpt["path"] = ["filt"]
    
    if "title" in plotOpt:
        plotOpt["title_label"] = f"{plotOpt['title']} "
    else:
        plotOpt["title_label"] = " "
    
    if plotType == "error":
        plotOpt["title_label"] += "Rotation Error"
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        plotOpt["title_label"] += "Rotation Error"
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        plotOpt["title_label"] += "Rotation Error"
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        plotOpt["title_label"] += "Rotation"
        plotOpt["components"] = ["roll", "pitch", "yaw"]
        plotOpt["transformation"] = lambda x: np.rad2deg(my_rotm2eul(x))
        plotOpt["show_cov"] = False
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_rotation_matrix(dataTot, t, cases_plot, caseOpt):
    """
    Plot rotation matrix elements
    """
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    
    f = plt.figure(figsize=(15, 12))
    ax = []
    
    for i_c in range(len(cases_plot)):
        c = cases_plot[i_c]
        res = dataTot[c]
        R = dataTot[c]["filt"]["R"]
        R_vec = R.reshape(9, -1)
        
        if c in caseOpt and "mean" in caseOpt[c]:
            resLineSpecs = get_ls_specs(caseOpt[c]["mean"])
        else:
            resLineSpecs = get_ls_specs({})
        
        for i in range(9):
            if i_c == 0:  # Create subplots only once
                ax.append(plt.subplot(3, 3, i + 1))
            
            # Create marker indices
            marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
            
            ax[i].plot(t, R_vec[i, :],
                      label=res["label"],
                      color=color_order[i_c],
                      linestyle=resLineSpecs["ls_c"],
                      marker=resLineSpecs["m_c"],
                      markevery=marker_indices)
            ax[i].grid(True)
    
    return f, ax

# ... existing code ...


def plot_rotation_v2(data, t, caseOpt, plotType="normal", extras=None):
    """
    Plot rotation version 2 (only works for error)
    """
    if extras is None:
        plotOpt = {}
        print("Using Extras")
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name_mean"] = "R"
    plotOpt["variable_name_var"] = "R_std"
    plotOpt["transformation"] = np.rad2deg
    
    if "title" in plotOpt:
        plotOpt["title_label"] = f"{plotOpt['title']} "
    else:
        plotOpt["title_label"] = ""
    
    plotOpt["title_label"] += "Rotation"
    plotOpt["y_unit"] = "[deg]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "error":
        return plot_three_components_error_v2(data, t, caseOpt, plotOpt)
    else:
        raise ValueError("wrong")


def plot_specific_force(data, t, cases_plot, caseOpt, plotType="normal"):
    """
    Plot specific force
    """
    plotOpt = {
        "variable_name_mean": "s",
        "variable_name_var": "s_std",
        "transformation": lambda x: x,
        "title_label": "Specific force",
        "y_unit": "[m/s^2]",
        "pred_type": "pred"
    }
    
    if plotType == "error":
        return plot_three_components_error(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "error_log":
        return plot_three_components_error_log(data, t, cases_plot, caseOpt, plotOpt)
    elif plotType == "norm_error_log":
        return plot_three_components_norm_error_log(data, t, cases_plot, caseOpt, plotOpt)
    else:
        return plot_three_components(data, t, cases_plot, caseOpt, plotOpt)


def plot_T_a(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot T_a matrix
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name_mean"] = "T_a"
    plotOpt["transformation"] = lambda x: x
    plotOpt["title_label"] = "T_a"
    plotOpt["y_unit"] = "[-]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "normal":
        return plot_nine_components(data, t, cases_plot, caseOpt, plotOpt)
    else:
        raise ValueError("Invalid choice")


def plot_T_g(data, t, cases_plot, caseOpt, plotType="normal", extras=None):
    """
    Plot T_g matrix
    """
    if extras is None:
        plotOpt = {}
    else:
        plotOpt = extras.copy()
    
    plotOpt["variable_name_mean"] = "T_g"
    plotOpt["transformation"] = lambda x: x
    plotOpt["title_label"] = "T_g"
    plotOpt["y_unit"] = "[-]"
    plotOpt["pred_type"] = "filt"
    
    if plotType == "normal":
        return plot_nine_components(data, t, cases_plot, caseOpt, plotOpt)
    else:
        raise ValueError("Invalid choice")


def plot_three_components(dataTot, t, cases_plot, caseOpt, plotOpt):
    """
    Core function to plot three components (x, y, z)
    """
    path_struct = plotOpt["path"]
    variable_name = plotOpt["variable_name"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    
    if "components" in plotOpt:
        directions = plotOpt["components"]
    else:
        directions = ["x", "y", "z"]
    
    if "show_cov" in plotOpt:
        show_cov = plotOpt["show_cov"]
    else:
        show_cov = True
    
    if "figure" in plotOpt:
        f = plotOpt["figure"]
    else:
        f = plt.figure(figsize=(10, 12))
    
    if "axes" in plotOpt:
        ax = plotOpt["axes"]
    else:
        ax = []
        for i in range(3):
            ax.append(plt.subplot(3, 1, i + 1))
    
    for i in range(3):
        leg = []
        for i_c in range(cases_plot.shape[0]):
            case_path = cases_plot[i_c, :]
            c = case_path[-1]
            path_struct_i = np.concatenate([case_path, path_struct])
            path_struct_i_label = np.concatenate([case_path, path_struct[:-1], ["label"]])
            
            try:
                # Navigate through nested dictionary structure
                res = dataTot
                for key in path_struct_i:
                    res = res[key]
            except KeyError:
                warnings.warn(f"{c} is not in {'/'.join(path_struct)}")
                continue
            
            try:
                # Get label
                label_res = dataTot
                for key in path_struct_i_label:
                    label_res = label_res[key]
                label_2 = label_res
            except KeyError:
                warnings.warn(f"label {c} is not in {'/'.join(path_struct)}")
                label_2 = ""
            
            try:
                label_1 = dataTot[case_path[0]]["label"]
            except KeyError:
                label_1 = ""
            
            label = f"{label_2}\n{label_1}"
            
            if c in caseOpt and "mean" in caseOpt[c]:
                resLineSpecs = get_ls_specs(caseOpt[c]["mean"])
            else:
                resLineSpecs = get_ls_specs({})
            
            data = f_trans(res["mean"][variable_name])
            
            # Create marker indices
            marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
            
            line1 = ax[i].plot(t, data[i, :],
                              label=label,
                              color=color_order[i_c],
                              linestyle=resLineSpecs["ls_c"],
                              marker=resLineSpecs["m_c"],
                              markevery=marker_indices)
            leg.extend(line1)
            
            if show_cov and variable_name in res["std"]:
                data_sig = f_trans(res["std"][variable_name])
                
                if c in caseOpt and "cov" in caseOpt[c]:
                    resLineSpecs = get_ls_specs(caseOpt[c]["cov"])
                else:
                    resLineSpecs = get_ls_specs({})
                    resLineSpecs["ls_c"] = "--"
                
                line2 = ax[i].plot(t, data[i, :] + 3 * data_sig[i, :],
                                  label=f"{label} $3\\sigma$",
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line2)
                
                ax[i].plot(t, data[i, :] - 3 * data_sig[i, :],
                          color=color_order[i_c],
                          linestyle=resLineSpecs["ls_c"],
                          marker=resLineSpecs["m_c"],
                          markevery=marker_indices)
        
        ax[i].grid(True)
        ax[i].set_title(title_label)
        
        if "legend" in plotOpt:
            leg_opt = plotOpt["legend"]
            ax[i].legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
        else:
            ax[i].legend(loc="best", ncol=1)
        
        ax[i].set_xlabel("Time [s]")
        ax[i].set_ylabel(f"{directions[i]} {y_unit}")
    
    return f, ax


def plot_three_components_error(dataTot, t, cases_plot, caseOpt, plotOpt):
    """
    Core function to plot three component errors
    """
    path_struct = plotOpt["path"]
    variable_name = plotOpt["variable_name"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    directions = ["x", "y", "z"]
    
    if "figure" in plotOpt:
        f = plotOpt["figure"]
    else:
        f = plt.figure(figsize=(10, 12))
    
    if "axes" in plotOpt:
        ax = plotOpt["axes"]
    else:
        ax = []
        for i in range(3):
            ax.append(plt.subplot(3, 1, i + 1))
    
    if "show_cov" in plotOpt:
        show_cov = plotOpt["show_cov"]
    else:
        show_cov = True
    
    if "show_mean" in plotOpt:
        show_mean = plotOpt["show_mean"]
    else:
        show_mean = True
    
    for i in range(3):
        leg = []
        for i_c in range(cases_plot.shape[0]):
            case_path = cases_plot[i_c, :]
            c = case_path[-1]
            path_struct_i = np.concatenate([case_path, path_struct])
            path_struct_i_label = np.concatenate([case_path, path_struct[:-1], ["label"]])
            
            try:
                # Navigate through nested dictionary structure
                res = dataTot
                for key in path_struct_i:
                    res = res[key]
            except KeyError:
                warnings.warn(f"{c} is not in {'/'.join(path_struct)}")
                continue
            
            try:
                # Get label
                label_res = dataTot
                for key in path_struct_i_label:
                    label_res = label_res[key]
                label_2 = label_res
            except KeyError:
                warnings.warn(f"label {c} is not in {'/'.join(path_struct)}")
                label_2 = ""
            
            try:
                label_1 = dataTot[case_path[0]]["label"]
            except KeyError:
                label_1 = ""
            
            label = f"{label_2}\n{label_1}"
            
            if c in caseOpt and "mean" in caseOpt[c]:
                resLineSpecs = get_ls_specs(caseOpt[c]["mean"])
            else:
                resLineSpecs = get_ls_specs({})
            
            if show_mean and "mean" in res and variable_name in res["mean"]:
                data = f_trans(res["mean"][variable_name])
                
                # Create marker indices
                marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                
                line1 = ax[i].plot(t, data[i, :],
                                  label=label,
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line1)
            
            if show_cov and variable_name in res["std"]:
                data_sig = f_trans(res["std"][variable_name])
                
                if c in caseOpt and "cov" in caseOpt[c]:
                    resLineSpecs = get_ls_specs(caseOpt[c]["cov"])
                else:
                    resLineSpecs = get_ls_specs({})
                    resLineSpecs["ls_c"] = "--"
                
                # Create marker indices
                marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                
                line2 = ax[i].plot(t, 3 * data_sig[i, :],
                                  label=f"{label} $3\\sigma$",
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line2)
                
                ax[i].plot(t, -3 * data_sig[i, :],
                          color=color_order[i_c],
                          linestyle=resLineSpecs["ls_c"],
                          marker=resLineSpecs["m_c"],
                          markevery=marker_indices)
        
        ax[i].grid(True)
        ax[i].set_title(title_label)
        
        if "legend" in plotOpt:
            leg_opt = plotOpt["legend"]
            ax[i].legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
        else:
            ax[i].legend(loc="best", ncol=1)
        
        ax[i].set_xlabel("Time [s]")
        ax[i].set_ylabel(f"{directions[i]} {y_unit}")
    
    return f, ax

# ... existing code ...


def plot_three_components_error_log(dataTot, t, cases_plot, caseOpt, plotOpt):
    """
    Plot three component errors on log scale
    """
    path_struct = plotOpt["path"]
    variable_name = plotOpt["variable_name"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    directions = ["x", "y", "z"]
    
    if "figure" in plotOpt:
        f = plotOpt["figure"]
    else:
        f = plt.figure(figsize=(10, 12))
    
    if "axes" in plotOpt:
        ax = plotOpt["axes"]
    else:
        ax = []
        for i in range(3):
            ax.append(plt.subplot(3, 1, i + 1))
    
    if "show_cov" in plotOpt:
        show_cov = plotOpt["show_cov"]
    else:
        show_cov = True
    
    if "show_mean" in plotOpt:
        show_mean = plotOpt["show_mean"]
    else:
        show_mean = True
    
    for i in range(3):
        leg = []
        for i_c in range(cases_plot.shape[0]):
            case_path = cases_plot[i_c, :]
            c = case_path[-1]
            path_struct_i = np.concatenate([case_path, path_struct])
            path_struct_i_label = np.concatenate([case_path, path_struct[:-1], ["label"]])
            
            try:
                # Navigate through nested dictionary structure
                res = dataTot
                for key in path_struct_i:
                    res = res[key]
            except KeyError:
                warnings.warn(f"{c} is not in {'/'.join(path_struct)}")
                continue
            
            try:
                # Get label
                label_res = dataTot
                for key in path_struct_i_label:
                    label_res = label_res[key]
                label_2 = label_res
            except KeyError:
                warnings.warn(f"label {c} is not in {'/'.join(path_struct)}")
                label_2 = ""
            
            try:
                label_1 = dataTot[case_path[0]]["label"]
            except KeyError:
                label_1 = ""
            
            label = f"{label_2}\n{label_1}"
            
            if c in caseOpt and "mean" in caseOpt[c]:
                resLineSpecs = get_ls_specs(caseOpt[c]["mean"])
            else:
                resLineSpecs = get_ls_specs({})
            
            if show_mean and "mean" in res and variable_name in res["mean"]:
                data = np.abs(f_trans(res["mean"][variable_name]))
                
                # Create marker indices
                marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                
                line1 = ax[i].plot(t, data[i, :],
                                  label=label,
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line1)
            
            if show_cov and variable_name in res["std"]:
                data_sig = f_trans(res["std"][variable_name])
                
                if c in caseOpt and "cov" in caseOpt[c]:
                    resLineSpecs = get_ls_specs(caseOpt[c]["cov"])
                else:
                    resLineSpecs = get_ls_specs({})
                    resLineSpecs["ls_c"] = "--"
                
                # Create marker indices
                marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                
                line2 = ax[i].plot(t, 3 * data_sig[i, :],
                                  label=f"{label} $3\\sigma$",
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line2)
        
        ax[i].grid(True)
        ax[i].set_yscale('log')
        ax[i].set_title(title_label)
        
        if "legend" in plotOpt:
            leg_opt = plotOpt["legend"]
            ax[i].legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
        else:
            ax[i].legend(loc="best", ncol=1)
        
        ax[i].set_xlabel("Time [s]")
        ax[i].set_ylabel(f"{directions[i]} {y_unit}")
    
    return f, ax


def plot_three_components_error_v2(dataTot, t, caseOpt, plotOpt):
    """
    Plot three component errors version 2
    """
    pred_type = plotOpt["pred_type"]
    v_mean = plotOpt["variable_name_mean"]
    v_var = plotOpt["variable_name_var"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    directions = ["x", "y", "z"]
    
    if "figure" in plotOpt:
        f = plotOpt["figure"]
    else:
        f = plt.figure(figsize=(10, 12))
    
    if "axes" in plotOpt:
        ax = plotOpt["axes"]
    else:
        ax = []
        for i in range(3):
            ax.append(plt.subplot(3, 1, i + 1))
    
    if "show_cov" in plotOpt:
        show_cov = plotOpt["show_cov"]
    else:
        show_cov = True
    
    for i in range(3):
        leg = []
        for i_c in range(len(dataTot)):
            res = dataTot[i_c][0] if isinstance(dataTot[i_c], list) else dataTot[i_c]
            caseOpt_i = caseOpt[i_c]
            
            if pred_type not in res:
                warnings.warn(f"{pred_type} not in case {i_c}")
                continue
            
            res_pred = res[pred_type]
            
            if isinstance(caseOpt_i, dict) and "mean" in caseOpt_i:
                resLineSpecs = get_ls_specs(caseOpt_i["mean"])
            else:
                resLineSpecs = get_ls_specs({})
            
            if "err" in res and v_mean in res["err"]:
                data = f_trans(res["err"][v_mean])
                
                # Create marker indices
                marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                
                line1 = ax[i].plot(t, data[i, :],
                                  label=res["label"],
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line1)
            
            if show_cov and v_var in res_pred:
                data_sig = f_trans(res_pred[v_var])
                
                if isinstance(caseOpt_i, dict) and "cov" in caseOpt_i:
                    resLineSpecs = get_ls_specs(caseOpt_i["cov"])
                else:
                    resLineSpecs = get_ls_specs({})
                    resLineSpecs["ls_c"] = "--"
                
                # Create marker indices
                marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                
                line2 = ax[i].plot(t, 3 * data_sig[i, :],
                                  label=f"{res['label']} $3\\sigma$",
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line2)
                
                ax[i].plot(t, -3 * data_sig[i, :],
                          color=color_order[i_c],
                          linestyle=resLineSpecs["ls_c"],
                          marker=resLineSpecs["m_c"],
                          markevery=marker_indices)
        
        ax[i].grid(True)
        ax[i].set_title(title_label)
        
        if "legend" in plotOpt:
            leg_opt = plotOpt["legend"]
            ax[i].legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
        else:
            ax[i].legend(loc="best", ncol=1)
        
        ax[i].set_xlabel("Time [s]")
        ax[i].set_ylabel(f"{directions[i]} {y_unit}")
    
    return f, ax


def plot_three_components_filt_pred(dataTot, t, cases_plot, caseOpt, plotOpt):
    """
    Plot three components with both filtered and predicted data
    """
    pred_types = ["pred", "filt"]
    v_mean = plotOpt["variable_name_mean"]
    v_var = plotOpt["variable_name_var"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    
    def get_color(i):
        return color_order[mod1(i, len(color_order)) - 1]
    
    if "components" in plotOpt:
        directions = plotOpt["components"]
    else:
        directions = ["x", "y", "z"]
    
    if "show_cov" in plotOpt:
        show_cov = plotOpt["show_cov"]
    else:
        show_cov = True
    
    if "figure" in plotOpt:
        f = plotOpt["figure"]
    else:
        f = plt.figure(figsize=(10, 12))
    
    if "axes" in plotOpt:
        ax = plotOpt["axes"]
    else:
        ax = []
        for i in range(3):
            ax.append(plt.subplot(3, 1, i + 1))
    
    N_c = len(cases_plot)
    N_p = len(pred_types)
    
    for i in range(3):
        leg = []
        for i_c in range(N_c):
            c = cases_plot[i_c]
            res = dataTot[c]
            for i_p in range(N_p):
                ii = i_p * N_c + i_c
                pred_type = pred_types[i_p]
                
                if pred_type not in dataTot[c]:
                    continue
                
                res_pred = dataTot[c][pred_type]
                
                if (c in caseOpt and pred_type in caseOpt and 
                    c in caseOpt[pred_type] and "mean" in caseOpt[pred_type][c]):
                    resLineSpecs = get_ls_specs(caseOpt[pred_type][c]["mean"])
                else:
                    resLineSpecs = get_ls_specs({})
                
                if v_mean in res_pred:
                    data = f_trans(res_pred[v_mean])
                    
                    # Create marker indices
                    marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                    
                    line1 = ax[i].plot(t, data[i, :],
                                      label=f"{res['label']} {pred_type}",
                                      color=get_color(ii + 1),
                                      linestyle=resLineSpecs["ls_c"],
                                      marker=resLineSpecs["m_c"],
                                      markevery=marker_indices)
                    leg.extend(line1)
                
                if show_cov and v_var in res_pred:
                    data_sig = f_trans(res_pred[v_var])
                    
                    if (c in caseOpt and pred_type in caseOpt and 
                        c in caseOpt[pred_type] and "cov" in caseOpt[pred_type][c]):
                        resLineSpecs = get_ls_specs(caseOpt[pred_type][c]["cov"])
                    else:
                        resLineSpecs = get_ls_specs({})
                        resLineSpecs["ls_c"] = "--"
                    
                    # Create marker indices
                    marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                    
                    line2 = ax[i].plot(t, data[i, :] + 3 * data_sig[i, :],
                                      label=f"{res['label']} {pred_type} $3\\sigma$",
                                      color=get_color(ii + 1),
                                      linestyle=resLineSpecs["ls_c"],
                                      marker=resLineSpecs["m_c"],
                                      markevery=marker_indices)
                    leg.extend(line2)
                    
                    ax[i].plot(t, data[i, :] - 3 * data_sig[i, :],
                              color=get_color(ii + 1),
                              linestyle=resLineSpecs["ls_c"],
                              marker=resLineSpecs["m_c"],
                              markevery=marker_indices)
        
        ax[i].grid(True)
        ax[i].set_title(title_label)
        
        if "legend" in plotOpt:
            leg_opt = plotOpt["legend"]
            ax[i].legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
        else:
            ax[i].legend(loc="best", ncol=1)
        
        ax[i].set_xlabel("Time [s]")
        ax[i].set_ylabel(f"{directions[i]} {y_unit}")
    
    return f, ax


def plot_three_components_norm_error_log(dataTot, t, cases_plot, caseOpt, plotOpt):
    """
    Plot norm of three component errors on log scale
    """
    pred_type = plotOpt["pred_type"]
    v_mean = plotOpt["variable_name_mean"]
    v_var = plotOpt["variable_name_var"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    
    f, ax = plt.subplots(figsize=(10, 6))
    
    leg = []
    for i_c in range(len(cases_plot)):
        c = cases_plot[i_c]
        res = dataTot[c]
        res_pred = dataTot[c][pred_type]
        
        if c in caseOpt and "mean" in caseOpt[c]:
            resLineSpecs = get_ls_specs(caseOpt[c]["mean"])
        else:
            resLineSpecs = get_ls_specs({})
        
        if "err" not in res:
            warnings.warn(f"err not in {c}")
            continue
        
        if v_mean not in res["err"]:
            warnings.warn(f"{v_mean} not in {c}")
            continue
        
        data = f_trans(res["err"][v_mean])
        
        # Create marker indices
        marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
        
        line1 = ax.plot(t, norm_time(data),
                       label=res["label"],
                       color=color_order[i_c],
                       linestyle=resLineSpecs["ls_c"],
                       marker=resLineSpecs["m_c"],
                       markevery=marker_indices)
        leg.extend(line1)
        
        if v_var in res_pred:
            data_sig = f_trans(res_pred[v_var])
            
            if c in caseOpt and "cov" in caseOpt[c]:
                resLineSpecs = get_ls_specs(caseOpt[c]["cov"])
            else:
                resLineSpecs = get_ls_specs({})
                resLineSpecs["ls_c"] = "--"
            
            # Create marker indices
            marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
            
            line2 = ax.plot(t, 3 * norm_time(data_sig),
                           label=f"{res['label']}\n$3 \\sqrt{{\\sum \\sigma_i }}$",
                           color=color_order[i_c],
                           linestyle=resLineSpecs["ls_c"],
                           marker=resLineSpecs["m_c"],
                           markevery=marker_indices)
            leg.extend(line2)
    
    ax.grid(True)
    ax.set_yscale('log')
    ax.set_title(title_label)
    
    if "legend" in plotOpt:
        leg_opt = plotOpt["legend"]
        ax.legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
    else:
        ax.legend(loc="best", ncol=1)
    
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(y_unit)
    
    return f, ax


def plot_three_components_update(ax, dataTot, t, cases_plot, caseOpt, plotOpt):
    """
    Update existing axes with three component plots
    """
    pred_type = plotOpt["pred_type"]
    v_mean = plotOpt["variable_name_mean"]
    v_var = plotOpt["variable_name_var"]
    f_trans = plotOpt["transformation"]
    title_label = plotOpt["title_label"]
    y_unit = plotOpt["y_unit"]
    
    N_end = len(t)
    color_order = plt.cm.tab10(np.linspace(0, 1, 10))
    directions = ["x", "y", "z"]
    
    for i in range(3):
        leg = []
        for i_c in range(len(cases_plot)):
            c = cases_plot[i_c]
            res = dataTot[c]
            res_pred = dataTot[c][pred_type]
            
            if c in caseOpt and "mean" in caseOpt[c]:
                resLineSpecs = get_ls_specs(caseOpt[c]["mean"])
            else:
                resLineSpecs = get_ls_specs({})
            
            data = f_trans(res_pred[v_mean])
            
            # Create marker indices
            marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
            
            line1 = ax[i].plot(t, data[i, :],
                              label=res["label"],
                              color=color_order[i_c],
                              linestyle=resLineSpecs["ls_c"],
                              marker=resLineSpecs["m_c"],
                              markevery=marker_indices)
            leg.extend(line1)
            
            if v_var in res_pred:
                data_sig = f_trans(res_pred[v_var])
                
                if c in caseOpt and "cov" in caseOpt[c]:
                    resLineSpecs = get_ls_specs(caseOpt[c]["cov"])
                else:
                    resLineSpecs = get_ls_specs({})
                
                # Create marker indices
                marker_indices = range(resLineSpecs["m_offset"] - 1, N_end, resLineSpecs["m_space"])
                
                line2 = ax[i].plot(t, data[i, :] + 3 * data_sig[i, :],
                                  label=f"{res['label']} $3\\sigma$",
                                  color=color_order[i_c],
                                  linestyle=resLineSpecs["ls_c"],
                                  marker=resLineSpecs["m_c"],
                                  markevery=marker_indices)
                leg.extend(line2)
                
                ax[i].plot(t, data[i, :] - 3 * data_sig[i, :],
                          color=color_order[i_c],
                          linestyle=resLineSpecs["ls_c"],
                          marker=resLineSpecs["m_c"],
                          markevery=marker_indices)
        
        ax[i].grid(True)
        ax[i].set_title(title_label)
        
        if "legend" in plotOpt:
            leg_opt = plotOpt["legend"]
            ax[i].legend(**dict(zip(["loc", "ncol"], leg_opt[:2])))
        else:
            ax[i].legend(loc="best", ncol=1)
        
        ax[i].set_xlabel("Time [s]")
        ax[i].set_ylabel(f"{directions[i]} {y_unit}")
    
    return ax