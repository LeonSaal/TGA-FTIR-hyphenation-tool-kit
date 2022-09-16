from typing import Iterable

import PySimpleGUI as sg

from ..config import lang


def mass_step_window(T_max, I_max):
    sliders = ["rel_height", "prominence", "height", "width"]
    bounds = [I_max, I_max, I_max, T_max]
    init = [0.9, 0, 0, 100]
    resolution = [0.001, 0.001, 0.001, 1]
    layout = [
        *[
            [
                sg.T(name),
                sg.Push(),
                sg.Slider(
                    range=(0, bound),
                    default_value=init,
                    resolution=res,
                    k=name,
                    orientation="h",
                ),
            ]
            for name, bound, init, res in zip(sliders, bounds, init, resolution)
        ],
        [sg.CBox("plot", default=True, k="plot")],
        [sg.OK(k="-K-"), sg.Push(), sg.Cancel(k="-X-")],
    ]

    window = sg.Window(lang.mass_step, layout=layout)
    while True:
        event, values = window.read()

        if event in [sg.WIN_CLOSED, "-X-"]:
            break

        if event == "-K-":
            window.close()
            return {key: val for key, val in values.items() if not key.startswith("-")}

    window.close()


def get_value_window(sample):
    layout = [
        [
            sg.T("which"),
            sg.Combo(
                sample.tga.columns.to_list(), default_value="sample_mass", k="which"
            ),
        ],
        [
            sg.T("at"),
            sg.Combo(
                sample.tga.columns.to_list(), default_value="reference_temp", k="at"
            ),
        ],
        [sg.T("values"), sg.Input("", k="-values")],
        [sg.OK(k="-K-")],
        [
            sg.Table(
                [],
                headings=["---", "---"],
                col_widths=[12, 12],
                auto_size_columns=False,
                justification="left",
                k="-VALUES-",
            )
        ],
    ]

    window = sg.Window(lang.get_value, layout=layout)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break
        if event == "-K-":
            if values["-values"]:
                args = {
                    key: val for key, val in values.items() if not key.startswith("-")
                }
                data = sample.get_value(*eval(f"{values['-values']},"), **args)
                data =  [
                    [idx, val] for idx, val in data.T.itertuples()
                ] #[[values["which"], values["at"]]] +

                print(data)
                window["-VALUES-"].update(values=data)

    window.close()


def corr_window(gases: Iterable):
    layout = [[sg.T(f'{lang.tga}:'),sg.CBox(lang.tga, k='tga'), sg.CBox(lang.dry_weight, k='dry_weight')], [sg.T(f'{lang.ir}:')]+[sg.CBox(gas.upper(), k=gas) for gas in gases], [sg.OK(k='-K-')]]

    window = sg.Window(lang.corr_plot, layout=layout)

    while True:
        event, values= window.read()

        if event == sg.WIN_CLOSED:
            break
        
        if event == '-K-':
            out = {}
            out['tga'] = values['tga']
            out['dry_weight'] = values['dry_weight']
            out['ir'] = [gas for gas in gases if values[gas]]
            window.close()
            return out

    window.close()
    return False
