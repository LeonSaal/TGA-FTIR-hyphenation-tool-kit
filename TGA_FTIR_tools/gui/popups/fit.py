import PySimpleGUI as sg
from ..config import lang
from typing import Mapping
from ...fitting import get_presets


def fit_set_window(name) -> Mapping:
    layout = [
        [sg.T(lang.ref), sg.Combo(get_presets(None), k="-REF-")],
        [sg.CBox("show plot", k="-PLOT-", default=True)],
        [sg.B(name, k="-GO-"), sg.Cancel(k="-X-")],
    ]
    window = sg.Window(lang.fit, layout, modal=True)
    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "-X-"]:
            break
        if event == "-GO-":
            window.close()
            return {"reference": values["-REF-"], "plot": values["-PLOT-"]}
    window.close()
