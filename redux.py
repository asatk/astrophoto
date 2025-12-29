from astrophoto.image import calib_image
from astrophoto.image.image import stack_image
from astrophoto.io import println, prompt, init_parser, init_defaults, store_defaults, set_quiet
from astrophoto.bias import calib_bias
from astrophoto.dark import calib_dark
from astrophoto.flat import calib_flat
from astrophoto.io.reduxio import save_and_quit

procedures = [
        "bias",
        "dark",
        "flat",
        "image",
        "stack"
    ]

if __name__ == "__main__":

    parser = init_parser()
    args = vars(parser.parse_args())

    set_quiet(args["quiet"])

    cfg_infile = args["config"]
    init_defaults(cfg_infile)

    proc = args["process"]
    if proc is None:
        proc = input(f"Data reduction process ({procedures}): ")

    println()

    root_dir = prompt("Project's root directory (path)", "root")

    println()

    if proc == "bias":
        calib_bias(root_dir)
    elif proc == "dark":
        calib_dark(root_dir)
    elif proc == "flat":
        calib_flat(root_dir)
    elif proc == "image":
        calib_image(root_dir)
    elif proc == "stack":
        stack_image(root_dir)

    println()

    save_and_quit(0)
