from typing import Iterable, Mapping

import PySimpleGUI as sg

from ...config import PLOTTING
from ..config import lang


def plot_set_window(plot: str, gases: Iterable, one_gas=False) -> Mapping:
    cboxes = ["save", "title", "legend"]
    settings = {
        "x_axis": sg.Column(
            [
                [
                    sg.T("x_axis"),
                    sg.Combo(
                        ["sample_temp", "time"], k="x_axis"
                    ),
                ]
            ]
        ),
        "y_axis": sg.Column(
            [
                [
                    sg.T("y_axis"),
                    sg.Combo(["rel", "orig"], k="y_axis"),
                ]
            ]
        ),
        "xlim": sg.Column(
            [[sg.T("xlim"), sg.Input(k="xlim")]]
        ),
        "ylim": sg.Column(
            [[sg.T("ylim"), sg.Combo([None, "auto"], k="ylim")]]
        ),
    }
    if one_gas:
        settings["gases"] = (
            sg.Column(
                [
                    [sg.T("show gases")]
                    + [sg.Radio(gas, 0, k=f"-{gas}", default=(i==0)) for i, gas in enumerate(gases)]
                ]
            ),
        )
    else:
        settings["gases"] = (
            sg.Column(
                [[sg.T("show gases")] + [sg.CBox(gas, k=f"-{gas}") for gas in gases]]
            ),
        )

    settings.update({name: sg.CBox(name, k=name) for name in cboxes})
    base_opts = ["xlim", "x_axis", "y_axis", "title", "legend", "save"]
    if plot in [lang.tga, lang.heatflow]:
        options = ["ylim"]

    if plot in [lang.ir, lang.irdtg, lang.dir]:
        options = ["gases"]

    if plot in ["fit", "robustness"]:
        options = []

    options += base_opts
    layout = (
        [
            [
                sg.T(
                    f"{lang.plot}: {plot}",
                )
            ]
        ]
        + [[settings[opt]] for opt in options]
        + [[sg.B(lang.OK, k='-K-'), sg.Cancel(k="-X-")]]
    )

    if PLOTTING['dialog'] =='False':
        out = {key: PLOTTING[key] if (PLOTTING[key] not in ['True', 'False'] and key != 'xlim') else eval(PLOTTING[key]) for key in options}
        return out

    window = sg.Window(lang.plot, layout, modal=True, finalize=True)

    for elem  in window.__dict__['AllKeysDict'].keys():
        if not elem.startswith('-'):
            if (val:=PLOTTING[elem]) in ['True', 'False']:
                default = eval(val)
            else:
                default = val
            window[elem].update(default)


    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "-X-"]:
            break
        if event == '-K-':
            window.close()
            out = {
                name: val if name != "xlim" else eval(val)
                for name, val in values.items()
                if not name.startswith("-")
            }
            if plot in [lang.ir, lang.irdtg]:
                out["gases"] = [gas for gas in gases if values[f"-{gas}"]]
                if one_gas:
                    out["gas"], = out["gases"]
                    del out["gases"]
            return out
    window.close()
