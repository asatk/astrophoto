from astrophoto.image import calib_image
from astrophoto.image.image import stack_image
from astrophoto.io import println, prompt, init_parser, init_defaults, store_defaults, set_quiet
from astrophoto.bias import calib_bias
from astrophoto.dark import calib_dark
from astrophoto.flat import calib_flat

if __name__ == "__main__":

    parser = init_parser()
    args = vars(parser.parse_args())

    proc = args["process"]
    set_quiet(args["quiet"])

    procedures = [
        "bias",
        "dark",
        "flat",
        "image",
        "stack"
    ]

    cfg_infile = input("Config file (path)? ")
    init_defaults(cfg_infile)

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

    cfg_outfile = input("Save responses to config file (path)? ")
    store_defaults(cfg_outfile)
