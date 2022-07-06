import PySimpleGUI as sg
from ..config import lang
import pandas as pd


def samplelog_window(data: pd.DataFrame):
    layout = [
        [sg.T(lang.sel_rows)],
        [
            sg.Table(
                [[*vals] for vals in data[["alias", "reference"]].itertuples()],
                headings=[lang.sample, lang.alias, lang.ref],
                select_mode="extended",
                enable_click_events=True,
                expand_x=True,
                justification="left",
                k="-LOG-",
            )
        ],
        [sg.OK(k="-K-"), sg.Push(), sg.Cancel(k="-X-")],
    ]

    window = sg.Window(lang.samplelog, layout)

    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "-X-"]:
            break
        elif event == "-K-":
            window.close()
            data = window["-LOG-"].get()
            return [data[i][0] for i in values["-LOG-"]]

    window.close()
    return
