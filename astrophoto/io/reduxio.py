import argparse as ap
import json

COMBINE_MEDIAN = 1
COMBINE_MEAN = 2

defaults = {
    "root": ".",
    "fpattern": "*.fit",

    "bias_fout": "bias.fits",
    "bias_combine": COMBINE_MEDIAN,

    "dark_fout": "dark.fits",
    "dark_texp": 0.0,
    "dark_combine": COMBINE_MEDIAN,
    "dark_is_current": 0,
}

is_quiet = False



def init_parser():

    parser = ap.ArgumentParser(
        "redux",
        usage="%(prog)s [options]",
        description="""
        Reduces astronomical imaging data.
        """,
        formatter_class=ap.HelpFormatter
    )

    parser.add_argument("-q", "--quiet", action="store_const", const=True, default=False)
    parser.add_argument("process", nargs="?", default="data", type=str)

    return parser


def init_defaults(fname: str=None):
    global defaults
    if fname is not None and fname != "":
        with open(fname, "r") as f:
            json_data = json.load(f)
            defaults.update(**json_data)

def store_defaults(fname):
    if fname is not None and fname != "":
        with open(fname, "w") as f:
            json.dump(defaults, f, indent="\t")
    else:
        printf("Responses not saved.")

def set_quiet(q):
    global is_quiet
    is_quiet = q

def printf(msg):
    if not is_quiet:
        print(msg)

def println():
    printf("="*80)

def prompt(msg: str, key: str, cast_fn=str):
    val = defaults.get(key, "")
    resp = input(f"{msg} [{val}] ")

    if resp == "":
        resp = cast_fn(val)
    else:
        resp = cast_fn(resp.strip())
        defaults[key] = resp

    return resp

def prompt_choices(msg: str, choices: str, key: str, cast_fn=str):
    val = defaults.get(key, "")
    resp = input(f"{msg} [{val}]\n{choices}\n")

    if resp == "":
        resp = cast_fn(val)
    else:
        resp = cast_fn(resp.strip())
        defaults[key] = resp

    return resp
