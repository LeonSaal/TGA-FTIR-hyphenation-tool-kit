import configparser
import logging
import os
import re
import webbrowser
from heapq import merge
from pathlib import Path
from typing import Iterable, List, Mapping

import matplotlib.pyplot as plt
import PySimpleGUI as sg

from ..calibration import calibrate
from ..classes import Sample, Worklist
from ..config import cfg, fmt
from ..input_output import samplelog
from ..links import LINKS
from .lang import EN
from .popups import (corr_window, fit_set_window, get_value_window,
                     mass_step_window, plot_set_window, samplelog_window,
                     setting_window)

plt.ion()

logger = logging.getLogger(__name__)
width = 50

lang = EN

opts = {
    lang.file: [lang.open, lang.load],
    lang.wl: [lang.load,  lang.from_samplelog, "-", lang.rename, lang.cl,"-", lang.save],
    lang.settings: [lang.edit, lang.load, lang.save],
    lang.cali: [lang.load, lang.recali],
    lang.recali: [lang.max, lang.iter, lang.co_oxi, lang.co_oxi_iter, lang.mlr],
    lang.info: [lang.about, lang.help]
}

opts_sample = {
    lang.plot: [lang.tga, lang.ir, lang.dir, lang.irdtg, lang.heatflow, lang.fit],
    lang.samples: [
        lang.plot,
        lang.corr,
        "-",
        lang.fit,
        "-",
        lang.get_value,
        lang.mass_step,
        "-",
        lang.delete,
        lang.save_wl,
    ],
}

opts_wl = {
    lang.plot: [lang.tga, lang.ir, lang.dir, lang.heatflow, lang.fit],
    lang.samples: [lang.plot,lang.corr,"-",lang.fit, lang.rob,"-",lang.delete,lang.save_wl],
}


def update_samples(wl):
    data = [
        [
            sample.name,
            sample.alias,
            sample.sample,
            sample.run,
            ", ".join(sample.info.gases),
            lang.heatflow in sample.tga,
        ]
        for sample in wl.samples
    ]
    return data


def make_menu(top: str, keys: Iterable, opts: Mapping):
    menu = []
    for key in keys:
        if key in opts:
            menu.extend([key, make_menu(key, opts[key], opts)])
        elif key == "-":
            menu.append("---")
        else:
            menu.append(f"{key}::{top}")
    return menu

def disable_menu(menu: list, disable: list):
    out =[]
    for i, item in enumerate(menu):
        if isinstance(item, list):
            out.append(disable_menu(item, disable))
        else:
            if item not in [f'{key}::{parent}' for key, parent in disable]:
                out.append(item)
    return out

menu_definition = [
    make_menu("", [opt], opts) for opt in [lang.file, lang.wl, lang.settings, lang.cali, lang.info]
]


def gui():
    # SETUP
    global cfg
    cols_wl = [lang.name, lang.alias, lang.sample, lang.run, lang.gases, lang.heatflow]
    worklist_frame = [
        [
            sg.Table(
                [],
                headings=cols_wl,
                key="-SAMPLES-",
                num_rows=6,
                justification="left",
                col_widths=[5, 5, width - 5],
                select_mode="extended",
                enable_events=True,
                enable_click_events=True,
                expand_x=True,
                tooltip=lang.tip_itable,
                right_click_menu=sg.MENU_RIGHT_CLICK_DISABLED,
            )
        ],
    ]

    layout = [
        [sg.MenuBar(menu_definition, key="-MENU-")],
        [sg.VPush()],
        [sg.Frame(lang.wl, worklist_frame, expand_x=True, k="-SAMPLELIST-")],
        [sg.VPush()],
        [
            sg.Multiline(
                k="-OUT-",
                #reroute_stdout=True,
                #reroute_stderr=True,
                autoscroll=True,
                disabled=True,
                size=(None, 10),
            )
        ],
        [sg.VPush()],
    ]

    window = sg.Window(lang.prog_name, layout)
    logging.basicConfig(level=logging.INFO, format=fmt, style="{", force=True)
    log = samplelog()

    wl = Worklist(lang.wl)
    subset = None
    import matplotlib.pyplot as plt

    # EVENT LOOP
    while True:
        # Converter
        event, values = window.read()

        if event == sg.WINDOW_CLOSED:
            break

        subset = wl[values["-SAMPLES-"]]

        if event:
            window["-OUT-"].update()

        if type(event) == str:
            event = str(event)
            # WORKLIST OPTIONS
            if event.endswith(lang.wl):
                if not os.path.exists(cfg["paths"]["output"]):
                    os.mkdir(cfg["paths"]["output"])
                if event.startswith(lang.load):
                    path = sg.popup_get_file(
                        "",
                        no_window=True,
                        initial_folder=cfg["paths"]["output"],
                        file_types=((".wkl", ["*.wkl"]),),
                    )
                    if path:
                        fname, _ = os.path.splitext(Path(path).name)
                        wl.load(fname=fname)
                        window["-SAMPLELIST-"].update(wl.name)

                if event.startswith(lang.save):
                    path = sg.popup_get_file(
                        "",
                        no_window=True,
                        initial_folder=cfg["paths"]["output"],
                        save_as=True,
                        file_types=((".wkl", ["*.wkl"]),),
                    )
                    if path:
                        fname, _ = os.path.splitext(Path(path).name)
                        wl.save(fname=fname)
                if event.startswith(lang.rename):
                    text = sg.popup_get_text(
                        f"Rename '{wl.name}' to:", default_text=wl.name
                    )
                    if text:
                        wl.name = text

                if event.startswith(lang.from_samplelog):
                    if names := samplelog_window(log):
                        for name in names:
                            wl.append(Sample(name))
                if event.startswith(lang.cl):
                    wl = Worklist(lang.wl)

                window["-SAMPLES-"].update(values=update_samples(wl))
                window["-SAMPLELIST-"].update(wl.name)

            # FILE OPTIONS
            if event.endswith(lang.file):
                if event.startswith(lang.open):
                    selection = sg.popup_get_file(
                        "",
                        initial_folder=cfg["paths"]["data"],
                        no_window=True,
                        multiple_files=True,
                        file_types=((".txt, .csv", ["*.txt", "*.csv"]),),
                    )
                    files = list(
                        set(
                            [
                                res.group()
                                for file in selection
                                if (
                                    res := re.search(
                                        r"\w{2}_\d\.\d_\d{5}_\d{3}", Path(file).name
                                    )
                                )
                            ]
                        )
                    )
                    for file in files:
                        wl.append(Sample(file))

                    window["-SAMPLES-"].update(values=update_samples(wl))

                if event.startswith(lang.load):
                    if not (path := os.path.exists(cfg["paths"]["output"])):
                        os.makedirs(path)
                    selection = sg.popup_get_file(
                        "",
                        initial_folder=cfg["paths"]["output"],
                        no_window=True,
                        multiple_files=True,
                        file_types=((".pkl", ["*.pkl"]),),
                    )
                    if selection:
                        for file in selection:
                            name, _ = os.path.splitext(Path(file).name)
                            wl.append(Sample(name, mode="pickle"))
                        window["-SAMPLES-"].update(values=update_samples(wl))

            # SETTINGS
            if event.endswith(lang.settings):
                fname = "settings.ini"
                if event.startswith(lang.load):
                    file = sg.popup_get_file(
                        "",
                        no_window=True,
                        file_types=((".ini", ["*.ini"]),),
                        initial_folder=cfg["paths"]["home"],
                        default_path=fname,
                    )
                    if file:
                        cfg = configparser.ConfigParser()
                        cfg.read(file)
                elif event.startswith(lang.edit):
                    setting_window()
                elif event.startswith(lang.save):
                    file = sg.popup_get_file(
                        "",
                        no_window=True,
                        file_types=((".ini", [".ini"]),),
                        initial_folder=cfg["paths"]["home"],
                        default_path=fname,
                        save_as=True,
                    )
                    if file:
                        with open(file, "w", encoding='latin-1') as configfile:
                            cfg.write(configfile)
            # CALIBRATION
            if event.endswith(lang.cali):
                if event.startswith(lang.load):
                    sg.popup_get_file(
                        "", no_window=True, initial_folder=cfg["paths"]["calibration"]
                    )

            if event.endswith(lang.recali):
                method, _ = event.split("::")
                calibrate(mode="recalibrate", method=method)


            # INFO
            if event.endswith(lang.info):
                if event.startswith(lang.about):
                    webbrowser.open(LINKS.REPO)
                if event.startswith(lang.help):
                    webbrowser.open(LINKS.WIKI)

            if subset:
                if len(subset) == 1:
                    sample = subset
                    menu = make_menu("", [lang.samples], opts_sample)

                    disable = []
                    if lang.fit not in sample.results:
                        disable.append((lang.fit, lang.plot))
                    for att in ['ir', 'tga']:
                        if att not in sample.__dict__:
                            disable.append((lang.__dict__[att], lang.plot))
                    if lang.heatflow not in sample.tga:
                        disable.append((lang.heatflow,lang.plot))
                    menu = disable_menu(menu, disable)    

                else:
                    menu = make_menu("", [lang.samples], opts_wl)
                    menu = disable_menu(menu, [(lang.fit, lang.plot)])

                window["-SAMPLES-"].set_right_click_menu(menu)
            else:
                continue

            # SAMPLELIST
            if event.endswith(lang.samples) and values["-SAMPLES-"]:
                if event.startswith(lang.fit):
                    reference = fit_set_window(lang.fit)
                    if reference:
                        subset.fit(**reference)

                elif event.startswith(lang.rob):
                    reference = fit_set_window(lang.rob)
                    if reference:
                        subset.robustness(**reference)

                elif event.startswith(lang.delete):
                    for i in reversed(values["-SAMPLES-"]):
                        wl.pop(i)
                    window["-SAMPLES-"].update(values=update_samples(wl))

                elif event.startswith(lang.corr):
                    refs = all(log.reference.loc[[sample.name for sample in subset]].notna())
                    if not refs:
                        path = sg.popup_get_file(
                            "", no_window=True, initial_folder=cfg["paths"]["data"]
                        )
                        if path:
                            ref, _ = os.path.splitext(Path(path).name)
                    else: 
                        ref = None
                    

                    if type(subset) == Worklist:
                        gases = list(
                            eval(
                                "&".join(
                                    [
                                        repr(set(sample.info.gases))
                                        for sample in subset
                                    ]
                                )
                            )
                        )
                    else:
                        gases = subset.info.gases

                    plot = corr_window(gases)   
                    print(plot)
                    wl.corr(references=ref, plot=plot)

                elif event.startswith(lang.save_wl):
                    subset = wl[values["-SAMPLES-"]]
                    if type(subset) == Worklist:
                        file_type = ".wkl"
                    else:
                        file_type = ".pkl"
                    path = sg.popup_get_file(
                        "",
                        default_path=subset.name,
                        no_window=True,
                        initial_folder=cfg["paths"]["output"],
                        save_as=True,
                        file_types=((f"{file_type}", [f"*{file_type}"]),),
                    )
                    if path:
                        fname, _ = os.path.splitext(Path(path).name)
                        subset.save(fname=fname, how="pickle")

                elif event.startswith(lang.mass_step):
                    sample = wl[values["-SAMPLES-"]]
                    settings = mass_step_window(max(sample.tga.sample_temp), 1)
                    if settings:
                        sample.mass_step(**settings)

                elif event.startswith(lang.get_value):
                    get_value_window(subset)
                    pass

            # PLOTTING
            if event.endswith(lang.plot) and values["-SAMPLES-"]:
                plot = event.split("::")[0]
                if type(subset) == Worklist:
                    gases = list(
                        eval(
                            "&".join(
                                [
                                    repr(set(sample.info.gases))
                                    for sample in subset
                                ]
                            )
                        )
                    )

                    one_gas = True
                else:
                    gases = subset.info.gases
                    one_gas = False
                settings = plot_set_window(plot, gases, one_gas=one_gas)
                if settings:
                    subset.plot(plot, **settings)
