from dataclasses import fields, replace
from rheoscale.config import RheoscaleConfig
from .base import widget_for_field
import ipywidgets as widgets

class ConfigWidget:
    def __init__(self, config: RheoscaleConfig):
        self.config = config
        self.widgets = {}

        for f in fields(config):
            if f.name.startswith("_"):
                continue

            value = getattr(config, f.name)
            self.widgets[f.name] = widget_for_field(f, value)

        self.ui = widgets.VBox(list(self.widgets.values()))

    def to_config(self) -> RheoscaleConfig:
        updates = {
            name: w.value
            for name, w in self.widgets.items()
        }
        return replace(self.config, **updates)
