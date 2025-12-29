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

    parser.add_argument("process", nargs="?", default=None, type=str)
    parser.add_argument("-c", "--config", type=str, default=None)
    parser.add_argument("-q", "--quiet", action="store_const", const=True, default=False)



    return parser


def init_defaults(fname: str=None):
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

def save_and_quit(code: int):
    cfg_outfile = input("Save responses to config file (path)? ")
    store_defaults(cfg_outfile)
    exit(code)

def update_default(key: str, val):
    defaults[key] = val

def get_default(key: str):
    return defaults.get(key)

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
    if is_quiet:
        return val

    resp = input(f"{msg} [{val}] ")

    try:

        if resp == "":
            resp = cast_fn(val)
        else:
            resp = cast_fn(resp.strip())
            defaults[key] = resp

    except ValueError:
        print("Bad input - quitting.")
        save_and_quit(40)

    return resp

def prompt_choices(msg: str, choices: str, key: str, cast_fn=str):
    val = defaults.get(key, "")
    if is_quiet:
        return val

    resp = input(f"{msg} [{val}]\n{choices}\n")

    try:
        if resp == "":
            resp = cast_fn(val)
        else:
            resp = cast_fn(resp.strip())
            defaults[key] = resp
    except ValueError:
        print("Bad input - quitting.")
        save_and_quit(40)

    return resp


