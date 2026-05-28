import ipywidgets as widgets
from typing import Any
from dataclasses import Field

def widget_for_field(f: Field, value: Any):
    default = f.default
    meta = f.metadata
    label = meta.get("label", f.name)

    if label == 'Dead Extremum':
        widget =  widgets.Dropdown(
            options=[('Min', 1), ('Max', 2)],
            value=1
        )
    elif f.type is bool:
        if default == False:
            value = 2
        else: 
            value = 1
        widget =  widgets.Dropdown(
            options=[('True', 1), ('False', 2)],
            value=value
        )

    elif f.type in (int, float) or "Optional" in str(f.type):
        widget =  widgets.FloatText(
            value=value
        )

    else:
        widget =  widgets.Text(
        value="" if value is None else str(value),
    )

    return labeled(widget, label)

def labeled(widget, label, help_text=None):
    label_widget = widgets.HTML(
        f"<b>{label}</b>" + (f"<br><small>{help_text}</small>" if help_text else "")
    )
    return widgets.VBox([label_widget, widget])
