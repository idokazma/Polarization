"""CLI dispatcher: ``python -m polarization {single,multiple,blockcrystals}``.

Replaces the MATLAB ``main.m`` entry point. Overrides are passed as
``--override key=value`` where the value is parsed as JSON, so lists, floats,
booleans, and strings all round-trip naturally:

    python -m polarization single \\
        --override lambda_=8e-7 \\
        --override omega_factors='[-3e-5, 0, 3e-5]'
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict


def _parse_override(value: str) -> tuple[str, object]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"override must be key=value, got {value!r}")
    key, raw = value.split("=", 1)
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        parsed = raw
    return key, parsed


def _overrides_to_dict(items) -> Dict[str, object]:
    return dict(items or [])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="polarization")
    sub = parser.add_subparsers(dest="target", required=True)

    for name in ("single", "multiple", "blockcrystals"):
        sp = sub.add_parser(name)
        sp.add_argument(
            "--override", action="append", type=_parse_override, default=[],
            help="key=value override (value parsed as JSON)",
        )
        sp.add_argument("--no-plot", action="store_true", help="skip plotting")
        if name == "blockcrystals":
            sp.add_argument("--data", default=None, help="path to arraypoints.mat")

    args = parser.parse_args(argv)
    overrides = _overrides_to_dict(getattr(args, "override", []))

    if args.target == "single":
        from polarization.examples.single_scatterer import single_scatterer
        single_scatterer(overrides=overrides, plot=not args.no_plot)
    elif args.target == "multiple":
        from polarization.examples.multiple_scatterers import multiple_scatterers
        multiple_scatterers(overrides=overrides, plot=not args.no_plot)
    elif args.target == "blockcrystals":
        from polarization.examples.block_crystals import block_crystals
        block_crystals(overrides=overrides, plot=not args.no_plot, data_path=args.data)
    else:
        parser.error(f"unknown target {args.target!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
