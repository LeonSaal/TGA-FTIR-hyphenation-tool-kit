import logging

import matplotlib
import PySimpleGUI as sg

from ..config import cfg, fmt, lang

logger = logging.getLogger(__name__)


def setting_window():
    tabs = []
    params = ["hwhm", "height"]
    for section, vals in cfg.items():
        if section == "DEFAULT":
            continue
        elif section == "paths":
            tab = [
                [
                    sg.T(name, expand_x=True),
                    sg.Input(val, k=f"-SET_{section}::{name}"),
                    sg.B(lang.browse, k=f"-B_{section}::{name}"),
                ]
                for name, val in vals.items()
            ]

        elif section == "correction":
            tab = []
            for name, value in vals.items():
                key = f"-SET_{section}::{name}"
                if key in ['tga_plot', 'dry_weight_plot']:
                    tab.append([sg.CBox(name, k=key)])
                elif key == 'ir_plot':
                    tab.append([sg.CBox(gas, k=gas) for gas in gases])
       
        elif section == "plotting":
            tab = []
            for name, val in vals.items():
                key = f"-SET_{section}::{name}"
                if name == 'mpl-style':
                    tab.append([
                    sg.T(name, expand_x=True),
                    sg.Combo(
                        matplotlib.style.available + ["default"],
                        default_value=val,
                        k=key,
                        ),
                    ])
                elif val in ['True', 'False']:
                    tab.append([sg.CBox(name,default=eval(val), k=key)])
                elif name.endswith('axis'):
                    if name.startswith('x'):
                        vals = ['sample_temp', 'time']
                    elif name.startswith('y'):
                        vals = ['orig', 'rel']
                    tab.append([sg.T(name),sg.Push(),sg.Combo(vals, default_value=val, k=key)])
                elif name== 'xlim':
                    tab.append([sg.T(name),sg.Push(),sg.Input(default_text=val, k=key)])
                elif name == 'ylim':
                    tab.append([sg.T(name),sg.Push(),sg.Combo([None, 'auto'], default_value=val, k=key)])


        elif section == "fitting":
            name, val = "tol_center", vals["tol_center"]
            tol_center = [
                sg.T(name),
                sg.Slider(
                    (0, float(val) * 2),
                    default_value=float(val),
                    orientation="h",
                    resolution=0.5,
                    k=f"-SET_{section}::tol_center",
                ),
            ]
            rows = [
                [
                    [sg.Push()],
                    [sg.T("min")],
                    [sg.Push(), sg.T("0"), sg.Push()],
                    [sg.T("max")],
                ]
            ]
            resolution = {"hwhm": 0.5, "height": 0.05}
            for param in params:
                val = [vals[f"{param}_{suff}"] for suff in ["min", "0", "max"]]
                rows.append(
                    [
                        [sg.T(param)],
                        [
                            sg.Input(
                                default_text=val[0],
                                size=4,
                                enable_events=True,
                                k=f"-SET_{section}::{param}_min",
                            )
                        ],
                        [
                            sg.Slider(
                                (val[0], val[2]),
                                default_value=val[1],
                                orientation="h",
                                resolution=resolution[param],
                                k=f"-SET_{section}::{param}_0",
                            )
                        ],
                        [
                            sg.Input(
                                default_text=val[2],
                                size=4,
                                enable_events=True,
                                k=f"-SET_{section}::{param}_max",
                            )
                        ],
                    ]
                )

            tab = [tol_center, [sg.Column([row[i] for row in rows]) for i in range(4)]]
        elif section == "savgol":
            resolution = {"window_length": 2, "polyorder": 1}
            ranges = {"window_length": (1, 301), "polyorder": (1, 5)}
            tab = [
                [
                    sg.T(name, expand_x=True),
                    sg.Slider(
                        range=ranges[name],
                        default_value=val,
                        k=f"-SET_{section}::{name}",
                        resolution=resolution[name],
                        orientation="h",
                    ),
                ]
                for name, val in vals.items()
            ]
        else:
            tab = [
                [sg.T(name, expand_x=True), sg.Input(val, k=f"-SET_{section}::{name}")]
                for name, val in vals.items()
            ]

        tabs.append([sg.Tab(section, tab)])

    layout = [
        [sg.TabGroup(tabs, k="-SET_TAB-", enable_events=True)],
        [
            sg.B(lang.OK, key="-SET_K-"),
            sg.B(lang.cancel, key="-SET_X-"),
            sg.Push(),
            sg.B(lang.apply, k="-SET_APP-"),
            sg.B(lang.reset, k="-SET_RESET-"),
        ],
    ]

    window = sg.Window(lang.settings, layout, modal=True)
    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "-SET_X-"]:
            break
        if type(event) == str:
            if event.startswith("-B_"):
                _, name = event.split("::")
                key = f'-SET_{values["-SET_TAB-"]}::{name}'
                folder = sg.popup_get_folder("", no_window=True)
                if folder:
                    window[key].update(folder)
        if event in ["-SET_K-", "-SET_APP-"]:
            for key, value in values.items():
                if key.startswith("-SET_") and "::" in key:
                    section, name = key.removeprefix("-SET_").split("::")
                    cfg[section][name] = str(value)
            matplotlib.style.use("default")
            matplotlib.style.use(cfg["plotting"]["mpl-style"])

            if event == "-SET_K-":
                break

        if event == "-SET_RESET-":
            for key, value in values.items():
                if key.startswith("-SET_") and "::" in key:
                    section, name = key.removeprefix("-SET_").split("::")
                    window[key].update(cfg[section][name])
        if event.endswith("_min") or event.endswith("_max"):
            for param in params:
                if (
                    not (lb := values[f"-SET_fitting::{param}_min"])
                    or not lb.replace(".", "", 1).isdigit()
                ):
                    l = cfg["fitting"][f"{param}_min"]
                else:
                    l = lb

                if (
                    not (ub := values[f"-SET_fitting::{param}_max"])
                    or not ub.replace(".", "", 1).isdigit()
                ):
                    u = cfg["fitting"][f"{param}_max"]
                else:
                    u = ub

                if l > u:
                    l, u = (
                        cfg["fitting"][f"{param}_min"],
                        cfg["fitting"][f"{param}_max"],
                    )
                    logger.warn("Upper bound must be bigger than lower bound.")
                window[f"-SET_fitting::{param}_0"].update(range=(l, u))
                window[f"-SET_fitting::{param}_min"].update(l)
                window[f"-SET_fitting::{param}_max"].update(u)

    window.close()
