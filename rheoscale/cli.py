import argparse
from pathlib import Path
from rheoscale.config import RheoscaleConfig  
from rheoscale.rheoscale_runner import RheoscaleRunner


#this is the cli version of the script
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rheoscale",
        description="Run RheoScale analysis with configurable parameters"
    )

    # --- Core ---
    parser.add_argument(
        "protein_name",
        help="Name of the protein being analyzed"
    )

    parser.add_argument(
        "-i", "--input-file",
        dest="input_file_name",
        type=str,
        help="Path to input CSV file"
    )

    parser.add_argument(
        "-n", "--num-positions",
        dest="number_of_positions",
        type=int,
        help="Number of positions to analyze"
    )

    # --- WT options ---
    parser.add_argument("--wt-val", type=float, help="WT value")
    parser.add_argument("--wt-error", type=float, help="WT error")
    parser.add_argument("--wt-name", type=str, default="WT")

    # --- Value ranges ---
    parser.add_argument("--min-val", type=float)
    parser.add_argument("--max-val", type=float)
    parser.add_argument("--error-val", type=float)

    # --- Binning ---
    parser.add_argument("--bins", dest="number_of_bins", type=int)
    parser.add_argument(
        "--dead-extremum",
        choices=["Min", "Max"],
        default="Min"
    )
    parser.add_argument("--neutral-binsize", type=float)

    # --- Flags ---
    parser.add_argument("--log-scale", action="store_true")
    parser.add_argument("--even-bins", action="store_true")
    parser.add_argument("--output-histograms", action="store_true")

    # --- Output ---
    parser.add_argument(
        "-o", "--output-dir",
        default="Rheoscale_analysis",
        help="Output directory"
    )

    # --- Config I/O ---
    parser.add_argument(
        "--save-config",
        type=Path,
        help="Save config to JSON file and exit"
    )

    parser.add_argument(
        "--load-config",
        type=Path,
        help="Load config from JSON file"
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # ---- Load from JSON if requested ----
    if args.load_config:
        config = RheoscaleConfig.from_json(args.load_config)
    else:
        config = RheoscaleConfig(
            protein_name=args.protein_name,
            input_file_name=args.input_file_name,
            number_of_positions=args.number_of_positions,
            log_scale=args.log_scale,
            WT_val=args.wt_val,
            WT_error=args.wt_error,
            WT_name=args.wt_name,
            min_val=args.min_val,
            max_val=args.max_val,
            error_val=args.error_val,
            number_of_bins=args.number_of_bins,
            dead_extremum=args.dead_extremum,
            neutral_binsize=args.neutral_binsize,
            output_dir=args.output_dir,
            output_histogram_plots=args.output_histograms,
            even_bins=args.even_bins,
        )

    # ---- Save config if requested ----
    if args.save_config:
        config.to_json(args.save_config)
        print(f"Config saved to {args.save_config}")
        return

    # ---- Run validations / setup ----
    config._validate_and_make_output()

    # ---- Hand off to analysis ----
    print("Running RheoScale with config:")
    print(config)
    if config.input_file_name is None:
        raise ValueError('For CLI you must have a file that you can upload to the config must have --input_file_name')
    runner = RheoscaleRunner(config)
    runner.run()


if __name__ == "__main__":
    main()
