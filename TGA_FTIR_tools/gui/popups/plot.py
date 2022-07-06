import PySimpleGUI as sg
from ..config import lang
from typing import Mapping, Iterable


def plot_set_window(plot: str, gases: Iterable, one_gas=False) -> Mapping:
    cboxes = ["save", "title", "legend"]
    settings = {
        "x_axis": sg.Column(
            [
                [
                    sg.T("x_axis"),
                    sg.Combo(
                        ["sample_temp", "time"], default_value="sample_temp", k="x_axis"
                    ),
                ]
            ]
        ),
        "y_axis": sg.Column(
            [
                [
                    sg.T("y_axis"),
                    sg.Combo(["rel", "orig"], default_value="orig", k="y_axis"),
                ]
            ]
        ),
        "xlim": sg.Column(
            [[sg.T("xlim"), sg.Input(default_text="None, None", k="xlim")]]
        ),
        "ylim": sg.Column(
            [[sg.T("ylim"), sg.Combo([None, "auto"], default_value="auto", k="ylim")]]
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

    layout = (
        [
            [
                sg.T(
                    f"{lang.plot}: {plot}",
                )
            ]
        ]
        + [[settings[opt]] for opt in options + base_opts]
        + [[sg.B(lang.OK, k=lang.OK), sg.Cancel(k="-X-")]]
    )
    window = sg.Window(lang.plot, layout, modal=True)

    while True:
        event, values = window.read()
        print(event, values)
        if event in [sg.WIN_CLOSED, "-X-"]:
            break
        if event == lang.OK:
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
