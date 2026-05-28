import ipywidgets as widgets
from IPython.display import display
from rheoscale.config import RheoscaleConfig
from .config_widgets import ConfigWidget
from rheoscale.rheoscale_runner import RheoscaleRunner

def launch():
    cw = ConfigWidget(
        RheoscaleConfig(protein_name="")
    )

    output = widgets.Output()
    run_button = widgets.Button(description="Run RheoScale")

    def on_run(_):
        output.clear_output()
        with output:
            try:
                config = cw.to_config()
                runner =RheoscaleRunner(config)
                runner.run()
                print("Analysis complete")
            except Exception as e:
                print("Error:", e)

    run_button.on_click(on_run)

    display(widgets.VBox([cw.ui, run_button, output]))
