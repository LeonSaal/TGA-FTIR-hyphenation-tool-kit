from imp import reload
import PySimpleGUI as sg
import matplotlib
from traitlets import default

# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:16:15 2022

@author: Leon
"""
import PySimpleGUI as sg
from dotmap import DotMap
from ..config import cfg, config, fmt
from ..classes import Sample, Worklist
import os
import copy
import re
from pathlib import Path
from ..fitting import get_presets
import logging
from .lang import EN
import matplotlib.pyplot as plt
import configparser
plt.ion()

logger = logging.getLogger(__name__)

lang = EN

fonts = DotMap()
fonts.head = "Calibri 15 bold"
fonts.sub = "Calibri 10 bold"

default_folder = ""  # r'C:\Users\Leon\tubCloud2\Shared\LC-OCD\Rohdaten'
width = 50


opts = {lang.plot:[lang.tga, lang.ir, lang.dir, lang.irdtg, lang.heatflow, lang.fit],
lang.file:[lang.open, lang.load, lang.save],
lang.wl: [lang.open, lang.load, lang.save],
lang.settings:[lang.edit, lang.load, lang.save],
lang.cali: [lang.load, lang.recali],
lang.samples:[lang.fit, lang.rob, lang.plot, lang.delete]
}

def update_samples(wl):
    data = [
        [
            sample.name,
            sample.alias,
            sample.sample,
            sample.run,
            ', '.join(sample.info.gases),
            lang.heatflow in sample.tga,
        ]
        for sample in wl.samples
    ]
    return data
    

def make_menu(top, keys,opts):
    menu = []
    for key in keys:
        if key in opts:
            menu.append([key, make_menu(key, opts[key],opts)])
        else:
            menu.append(f'{key}::{top}')
    return menu

def make_rcl_menu(top, keys,opts):
    menu = []
    for key in keys:
        if key in opts:
            menu.extend([key, make_rcl_menu(key, opts[key],opts)])
        else:
            menu.append(f'{key}::{top}')
    return menu
    

menu_definition = make_menu("",[lang.file,lang.wl, lang.settings, lang.cali],opts)
right_click_menu =make_rcl_menu("",[lang.samples], opts)

def setting_window():
    tabs = []
    params = ['hwhm', 'height']
    for section, vals in cfg.items():
        if section == "DEFAULT":
            continue
        elif section == "paths":
            tab = [[sg.T(name, expand_x=True),sg.Input(val, k=f"-SET_{section}::{name}"),sg.B(lang.browse, k=f"-B_{section}::{name}"),] for name, val in vals.items()]
        elif section == 'plotting':
            tab = [[sg.T(name, expand_x=True), sg.Combo(matplotlib.style.available+['default'], default_value=val, k=f"-SET_{section}::{name}")] for name, val in vals.items()]
        elif section == 'fitting':
            name, val= 'tol_center', vals['tol_center']
            tol_center = [sg.T(name), sg.Slider((0,float(val)*2), default_value=float(val), orientation='h', resolution=0.5, k=f'-SET_{section}::tol_center')]
            rows = [[[sg.Push()],[sg.T('min')],[sg.Push(),  sg.T('0'), sg.Push()],[ sg.T('max')]]]
            resolution={'hwhm':0.5,'height': 0.05}
            for param in params:
                val = [vals[f'{param}_{suff}'] for suff in ['min', '0', 'max']]
                rows.append([[sg.T(param)], [sg.Input(default_text=val[0], size=4, enable_events=True,k=f'-SET_{section}::{param}_min')], [sg.Slider((val[0], val[2]), default_value=val[1], orientation='h', resolution=resolution[param], k=f'-SET_{section}::{param}_0')], [sg.Input(default_text=val[2], size=4, enable_events=True, k=f'-SET_{section}::{param}_max')]])
            
            tab = [tol_center, [sg.Column([row[i] for row in rows]) for i in range(4) ]]
        elif section =='savgol':
            resolution = {'window_length':2, 'polyorder':1}
            ranges = {'window_length':(1,301), 'polyorder':(1,5)}
            tab=[[sg.T(name, expand_x=True), sg.Slider(range= ranges[name],default_value=val, k=f'-SET_{section}::{name}', resolution=resolution[name], orientation='h')] for name, val in vals.items()]
        else:
            tab = [[sg.T(name, expand_x=True), sg.Input(val, k=f'-SET_{section}::{name}')] for name, val in vals.items()]

        tabs.append([sg.Tab(section, tab)])

    layout = [
        [sg.TabGroup(tabs, k="-SET_TAB-", enable_events=True)],
        [sg.B(lang.OK, key="-SET_K-"), sg.B(lang.cancel, key="-SET_X-"),sg.Push(), sg.B(lang.apply, k='-SET_APP-'), sg.B(lang.reset, k='-SET_RESET-')],
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
        if event in ["-SET_K-",'-SET_APP-']:
            for key, value in values.items():
                if key.startswith("-SET_") and "::" in key:
                    section, name = key.removeprefix("-SET_").split("::")
                    cfg[section][name] = str(value)
            matplotlib.style.use(cfg['plotting']['mpl-style'])

            if event == '-SET_K-':
                break
                
        if event == '-SET_RESET-':
            for key, value in values.items():
                if key.startswith("-SET_") and "::" in key:
                    section, name = key.removeprefix("-SET_").split("::")
                    window[key].update(cfg[section][name])
        if event.endswith('_min') or event.endswith('_max'):
            for param in params:
                d, u = float(values[f'{param}_min']), float(values[f'{param}_max'])
                if d > u:
                    d, u = cfg['fitting'][f'{param}_min'], cfg['fitting'][f'{param}_min']
                    logger.warn('Upper bound must be bigger than lower bound.')
                window[f'{param}_slider'].update(range=(d, u))


    window.close()

def plot_set_window(plot, gases):
    cboxes = ['save', 'title', 'legend']
    settings = {
    'x_axis': sg.Column([[sg.T('x_axis'),sg.Combo(['sample_temp', 'time'], default_value='sample_temp',k='x_axis')]]),
    'y_axis': sg.Column([[sg.T('y_axis'),sg.Combo(['rel','orig'], default_value='orig',k='y_axis')]]),
    'xlim': sg.Column([[sg.T('xlim'), sg.Input(default_text= 'None, None',k='xlim')]]),
    'ylim': sg.Column([[sg.T('ylim'),sg.Combo([None, 'auto'], default_value='auto', k='ylim')]]),
    'gases' :sg.Column([[sg.T('show gases')]+[sg.CBox(gas, k=gas) for gas in gases]])
    }

    settings.update({name: sg.CBox(name, k=name) for name in cboxes})
    base_opts = ['xlim','x_axis','y_axis', 'title', 'legend','save']
    if plot  in [lang.tga, lang.heatflow]:
        options = ['ylim']
        
    if plot in [lang.ir, lang.irdtg]:
        options = ['gases']
    
    if plot in ['fit', 'robustness']:
        options=[]

    layout = [[sg.T(f'{lang.plot}: {plot}', )]]+[[settings[opt]] for opt in options+base_opts]+[[sg.B(lang.OK, k=lang.OK), sg.Cancel(k='-X-')]]
    window = sg.Window(lang.plot, layout, modal=True)

    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "-X-"]:
            break 
        if event==lang.OK:
            window.close()
            out = {name: val if name != 'xlim' else eval(val) for name, val in values.items() if not name.startswith('-')}
            if plot in [lang.ir, lang.irdtg]:
                out['gases']=[gas for gas in gases if values[gas]]
            return out
    window.close()

def fit_set_window(name):
    layout = [
        [sg.T(lang.ref),sg.Combo(get_presets(None), k="-REF-")],
        [sg.CBox('show plot', k='-PLOT-', default=True)],
        [sg.B(name, k="-GO-"), sg.Cancel(k="-X-")],
    ]
    window = sg.Window(lang.fit, layout, modal=True)

    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "-X-"]:
            break
        if event == "-GO-":
            window.close()
            return {'reference': values["-REF-"], 'plot': values["-PLOT-"]}
    window.close()


def gui():
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
                enable_click_events=True,
                expand_x=True,
                tooltip=lang.tip_itable,
                right_click_menu=right_click_menu,
            )
        ],
    ]

    # figure_frame = [
    #     [
    #         sg.T(lang.incl_signal, font=fonts.sub),
    #         sg.CBox("OC", default=True, key="-OC_P-"),
    #         sg.CBox("UV", default=True, key="-UV_P-"),
    #         sg.CBox("UV2", default=True, key="-UV2_P-"),
    #         sg.CBox("t", default=False, key="-T_P-"),
    #         sg.CBox(lang.bounds_int, default=False, k='-BOUNDS_INT-', disabled=True),
    #         sg.Push(),
    #         sg.Canvas(key="-CONTROLS-"),
    #         sg.B(lang.cl, key="-FIG_CLEAR-", visible=False),
    #     ],
    #     [
    #         sg.Column(
    #             layout=[
    #                 [
    #                     sg.Canvas(
    #                         key="-FIGURE-",
    #                         # it's important that you set this size
    #                         size=(400 * 2, 400),
    #                     )
    #                 ]
    #             ],
    #             background_color="#DAE0E6",
    #             pad=(0, 0),
    #         )
    #     ],
    # ]

    layout = [
        [sg.MenuBar(menu_definition, key="-MENU-")],
        [sg.VPush()],
        [sg.Frame(lang.wl, worklist_frame, expand_x=True)],
        [sg.VPush()],
        [
            sg.Multiline(
                k="-OUT-",
                # reroute_stdout=True,
                # reroute_stderr=True,
                autoscroll=True,
                disabled=True,
                size=(None, 10),
            )
        ],
        [sg.VPush()],
        [sg.Button(lang.convert, key="-RUN-"), sg.Push()],
    ]

    window = sg.Window(lang.prog_name, layout)
    logging.basicConfig(level=logging.INFO, format=fmt, style="{", force=True)

    wl = Worklist(lang.wl)
    while True:
        # Converter
        event, values = window.read()
        print(event)
        if event:
            window["-OUT-"].update()

        if event == sg.WINDOW_CLOSED:
            break

        if type(event) == str:
            if event.endswith(lang.wl):
                if event.startswith(lang.open):
                    sg.popup_get_file("", no_window=True)
                if event.startswith(lang.load):
                    sg.popup_get_file("", no_window=True)

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
                    # if not (path:=os.path.exists(config["paths"]["output"])):
                    #     os.mkdir(path)
                    selection = sg.popup_get_file(
                        "",
                        initial_folder=cfg["paths"]["output"],
                        no_window=True,
                        multiple_files=True,
                        file_types=((".pkl", ["*.pkl"]),),
                    )
                    if selection:
                        for file in selection:
                            name = os.path.splitext(Path(file).name)[0]
                            wl.append(Sample(name))
                        window["-SAMPLES-"].update(values=update_samples(wl))

            if event.endswith(lang.settings):
                fname = 'settings.ini'
                if event.startswith(lang.load):
                    file = sg.popup_get_file("", no_window=True, file_types=((".ini", ["*.ini"]),), initial_folder=cfg['paths']['home'], default_path=fname)
                    if file:
                        cfg = configparser.ConfigParser()
                        cfg.read(file)
                if event.startswith(lang.edit):
                    setting_window()
                if event.startswith(lang.save):
                    file = sg.popup_get_file("", no_window=True, file_types=((".ini", [".ini"]),), initial_folder=cfg['paths']['home'], default_path=fname, save_as=True)
                    if file: 
                        with open(file, "w") as configfile:
                            cfg.write(configfile)

            if event.endswith(lang.cali):
                if event.startswith(lang.load):
                    sg.popup_get_file("", no_window=True)

            if event.endswith(lang.samples) and values["-SAMPLES-"]:
                if event.startswith(lang.fit):
                    reference = fit_set_window(lang.fit)
                    if reference:
                        wl[values["-SAMPLES-"]].fit(**reference)
                if event.startswith(lang.rob):
                    reference = fit_set_window(lang.rob)
                    if reference:
                        wl[values["-SAMPLES-"]].robustness(**reference)
                if event.startswith(lang.delete):
                    for i in reversed(values["-SAMPLES-"]):
                        wl.pop(i)

                    window["-SAMPLES-"].update(values=update_samples(wl))
                
            if event.endswith(lang.plot) and values["-SAMPLES-"]:
                plot = event.split("::")[0]
                gases = wl[values["-SAMPLES-"]].info.gases
                settings = plot_set_window(plot,gases)
                wl[values["-SAMPLES-"]].plot(plot, **settings)
                

        if type(event) == tuple:
            if event[0] == "-SAMPLES-":
                print(values["-SAMPLES-"])
            


        
