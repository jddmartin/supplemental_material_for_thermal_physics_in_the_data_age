"""A Python program to create a plot of molecular heat capacities using
the data at:  http://webbook.nist.gov/chemistry/fluid/

Updated versions will be available at:
  https://github.com/jddmartin/supplemental_material_for_thermal_physics_in_the_data_age

Requires Python 3.8.10 or greater and the additional packages:
  numpy
  scipy
  matplotlib

These packages are normally installed using "pip".  Consult the documentation
relevant to your Python installation.

If run as a main program, a *.pdf figure file is produced:
  molecular_heat_capacities_generated.pdf
The matplotlib style sheet file: "molecular_heat_capacities.mplstyle" should
be in the same directory as this file, as this file controls the formatting
of the *.pdf figure file.

Alternately, this *.py file can be imported as a module, and the function
"make_plot" run to display the plot on the screen.

When the program is first run, the relevant data will be downloaded from
NIST.  But this data is stored (in the "downloaded_data" subdirectory) so that
subsequent runs do not require downloading the data.
"""

import os.path
from collections import namedtuple
from collections import OrderedDict
import pickle
import urllib.request, urllib.error, urllib.parse
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

"""Through inspection of the html page:
  http://webbook.nist.gov/chemistry/fluid/
we can determine the "key" for each species that we can obtain data for.
Here are keys and corresponding species:
"""
species = {
    "C7732185": "Water",
    "C7727379": "Nitrogen",
    "C1333740": "Hydrogen",
    "B5000001": "Parahydrogen",
    "B5000007": "Orthohydrogen",
    "C7782390": "Deuterium",
    "C7782447": "Oxygen",
    "C7782414": "Fluorine",
    "C630080": "Carbon monoxide",
    "C124389": "Carbon dioxide",
    "C10024972": "Dinitrogen monoxide",
    "C7789200": "Deuterium oxide",
    "C67561": "Methanol",
    "C74828": "Methane",
    "C74840": "Ethane",
    "C74851": "Ethene",
    "C74986": "Propane",
    "C115071": "Propene",
    "C74997": "Propyne",
    "C75194": "Cyclopropane",
    "C106978": "Butane",
    "C75285": "Isobutane",
    "C109660": "Pentane",
    "C78784": "2-Methylbutane",
    "C463821": "2,2-Dimethylpropane",
    "C110543": "Hexane",
    "C107835": "2-Methylpentane",
    "C110827": "Cyclohexane",
    "C142825": "Heptane",
    "C111659": "Octane",
    "C111842": "Nonane",
    "C124185": "Decane",
    "C112403": "Dodecane",
    "C7440597": "Helium",
    "C7440019": "Neon",
    "C7440371": "Argon",
    "C7439909": "Krypton",
    "C7440633": "Xenon",
    "C7664417": "Ammonia",
    "C7783542": "Nitrogen trifluoride",
    "C75694": "Trichlorofluoromethane (R11)",
    "C75718": "Dichlorodifluoromethane (R12)",
    "C75729": "Chlorotrifluoromethane (R13)",
    "C75730": "Tetrafluoromethane (R14)",
    "C75434": "Dichlorofluoromethane (R21)",
    "C75456": "Methane, chlorodifluoro- (R22)",
    "C75467": "Trifluoromethane (R23)",
    "C75105": "Methane, difluoro- (R32)",
    "C593533": "Fluoromethane (R41)",
    "C76131": "1,1,2-Trichloro-1,2,2-trifluoroethane (R113)",
    "C76142": "1,2-Dichloro-1,1,2,2-tetrafluoroethane (R114)",
    "C76153": "Chloropentafluoroethane (R115)",
    "C76164": "Hexafluoroethane (R116)",
    "C306832": "Ethane, 2,2-dichloro-1,1,1-trifluoro- (R123)",
    "C2837890": "Ethane, 1-chloro-1,2,2,2-tetrafluoro- (R124)",
    "C354336": "Ethane, pentafluoro- (R125)",
    "C811972": "Ethane, 1,1,1,2-tetrafluoro- (R134a)",
    "C1717006": "1,1-Dichloro-1-fluoroethane (R141b)",
    "C75683": "1-Chloro-1,1-difluoroethane (R142b)",
    "C420462": "Ethane, 1,1,1-trifluoro- (R143a)",
    "C75376": "Ethane, 1,1-difluoro- (R152a)",
    "C76197": "Octafluoropropane (R218)",
    "C431890": "1,1,1,2,3,3,3-Heptafluoropropane (R227ea)",
    "C431630": "1,1,1,2,3,3-Hexafluoropropane (R236ea)",
    "C690391": "1,1,1,3,3,3-Hexafluoropropane (R236fa)",
    "C679867": "1,1,2,2,3-Pentafluoropropane (R245ca)",
    "C460731": "1,1,1,3,3-Pentafluoropropane (R245fa)",
    "C115253": "Octafluorocyclobutane (RC318)",
    "C71432": "Benzene",
    "C108883": "Toluene",
    "C355259": "Decafluorobutane",
    "C678262": "Dodecafluoropentane",
    "C7446095": "Sulfur dioxide",
    "C7783064": "Hydrogen sulfide",
    "C2551624": "Sulfur hexafluoride",
    "C463581": "Carbonyl sulfide",
}


def download_data():
    """Download thermodynamic fluid data from:
      http://webbook.nist.gov/chemistry/fluid/
    for each species, putting the downloaded data files in the subdirectory:
      "./downloaded_data"
    """
    # The code requests data for each molecule over a range of temperatures
    # using a single URL (with a part that changes for each molecule).  The
    # form for the required URLs was determined by interactively generating
    # a manual request for the data of a single specific molecule using a
    # web-browser and then examining the specific URL.

    # for each species, request data and write it to file:
    keys = list(species.keys())
    for i, akey in enumerate(keys):
        print("Downloading %d of %d, key: %s" % (i + 1, len(keys), akey))
        raw = (
            r"http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID="
            + akey
            + r"&Type=IsoBar&Digits=5&P=0&THigh=1000&TLow=0&TInc=1&"
            + r"RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fl&HUnit=kJ%2Fmol&W"
            + r"Unit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
        )
        response = urllib.request.urlopen(raw)
        page_text = response.read()  # these pages are just plain text --- not html
        apath = os.path.join(os.path.dirname(__file__),
            "downloaded_data")

        if not os.path.exists(apath):
            os.makedirs(apath)

        fname_with_path = os.path.join(apath, akey + "_generated.txt")
        with open(fname_with_path, "w") as f:
            f.write(page_text.decode())


def create_pickle_file():
    """Creates a *.pickle file from downloaded thermodynamic data in subdirectory:
     downloaded_data/
    The data in downloaded_data should be downloaded prior to running this script,
    by running
     0_download_data.py
    """

    C = namedtuple("C", "pos full_name units short_name")

    amap = [
        C(0, "Temperature", "K", "temperature"),
        C(1, "Pressure", "MPa", "pressure"),
        C(2, "Density", "mol/l", "density"),
        C(3, "Volume", "l/mol", "volume"),
        C(4, "Internal Energy", "kJ/mol", "internal_energy"),
        C(5, "Enthalpy", "kJ/mol", "enthalpy"),
        C(6, "Entropy", "J/mol*K", "entropy"),
        C(7, "Cv", "J/mol*K", "cv"),
        C(8, "Cp", "J/mol*K", "cp"),
        C(9, "Sound Spd.", "m/s", "sound_speed"),
        C(10, "Joule-Thomson", "K/MPa", "joule_thomson"),
        C(11, "Viscosity", "uPa*s", "viscosity"),
        C(12, "Therm. Cond.", "W/m*K", "thermal_conductivity"),
        C(13, "Phase", "", "phase"),
    ]

    # read in data for each species and load it into a dictionary:
    d = OrderedDict()
    for akey, name in list(species.items()):
        d[akey] = OrderedDict()
        sd = d[akey]
        sd["name"] = np.string_(name)
        for column in amap:
            sd[column.short_name] = []

        apath = os.path.join(os.path.dirname(__file__),
            "downloaded_data")

        with open(os.path.join(apath, akey + "_generated.txt"), "r") as f:
            f.readline()
            for i, line in enumerate(f.readlines()):
                data = line.split()
                for column in amap:
                    try:
                        num = float(data[column.pos])
                    except ValueError:
                        num = float("nan")
                    sd[column.short_name].append(num)
        for column in amap:
            sd[column.short_name] = np.array(sd[column.short_name])

    with open(os.path.join(apath, "nist_fluid_data_generated.pickle"),
              "wb") as f:
        pickle.dump(d, f)

def load_pickled_data():
    apath = os.path.join(os.path.dirname(__file__),
        "downloaded_data")

    with open(os.path.join(apath, "nist_fluid_data_generated.pickle"),
              "rb") as f:
        d = pickle.load(f)
    return d

def make_plot():
    """Generate plot, without producing output file.
    In ipython it can be useful to:
      %run -n "name_of_this_script".py
    to avoid triggering "if __name__ == "__main__": code
    """

    d = load_pickled_data()

    # for clarity, don't plot:
    dont_plot = [
        "B5000001",  # "Parahydrogen"
        "B5000007",  # "Orthohydrogen"
        "C7789200",  # "Deuterium oxide"
        "C630080",  # "Carbon monoxide"
    ]

    # define species that will have customized labels on plot:
    ls = namedtuple(
        "ls",
        "nist_key nist_name line_color line_style alpha label " "position offset ha va",
    )
    specials = [
        ls(
            "C7440597",
            "Helium",
            "r",
            "-",
            0.75,
            "He, Ne, Ar, Kr, Xe",
            (700, 1.5*2),
            (0, 10),
            "left",
            "bottom",
        ),
        ls(
            "C7440019",
            "Neon",
            "r",
            "-",
            0.75,
            "",
            (0, 0),
            (0, 0),
            "center",
            "center",
        ),
        ls(
            "C7440371",
            "Argon",
            "r",
            "-",
            0.75,
            "",
            (0, 0),
            (0, 0),
            "center",
            "center",
        ),
        ls(
            "C7439909",
            "Krypton",
            "r",
            "-",
            0.75,
            "",
            (0, 0),
            (0, 0),
            "center",
            "center",
        ),
        ls(
            "C7440633",
            "Xenon",
            "r",
            "-",
            0.75,
            "",
            (0, 0),
            (0, 0),
            "center",
            "center",
        ),
        ls(
            "C1333740",
            "Hydrogen",
            "b",
            "-",
            0.75,
            "H$_2$",
            (131.0, 1.94*2),
            (15, -5),
            "left",
            "center",
        ),
        ls(
            "C7782390",
            "Deuterium",
            "g",
            "-",
            0.75,
            "D$_2$",
            (58, 2.23*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C7727379",
            "Nitrogen",
            "r",
            ":",
            0.75,
            "N$_2$",
            (980, 2.93*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C7782447",
            "Oxygen",
            "y",
            ":",
            0.75,
            "O$_2$",
            (980.0, 3.2*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C7782414",
            "Fluorine",
            "m",
            ":",
            0.75,
            "F$_2$",
            (300, 2.8*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C7732185",
            "Water",
            "b",
            "-",
            0.75,
            "H$_2$O",
            (1000, 3.97*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C7783064",
            "Hydrogen sulfide",
            "g",
            "-",
            0.75,
            "H$_2$S",
            (760, 4*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C74828",
            "Methane",
            "r",
            "-",
            0.75,
            "CH$_4$",
            (625, 5.5*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C124389",
            "Carbon dioxide",
            "m",
            "-",
            0.75,
            "CO$_2$",
            (1000, 5.54*2),
            (15, 0),
            "left",
            "center",
        ),
        ls(
            "C811972",
            "Ethane, 1,1,1,2-tetrafluoro- (R134a)",
            "b",
            "-",
            1.0,
            "R-134a",
            (173, 6.22*2),
            (-20, 0),
            "right",
            "center",
        ),
    ]

    plt.ion()
    fig = plt.figure(1)
    axes = fig.add_subplot(111, autoscale_on=False,
        xlim=(0, 1125), ylim=(0, 8*2))

    axes.set_xlabel("temperature (K)")
    yl = axes.set_ylabel(r"$\dfrac{C_V}{Nk_B/2}$",
        labelpad=20, rotation=0)

    plt.minorticks_on()
    axes.tick_params(axis="y", which="minor", left="off")

    plotted = []  # to keep track of what we have plotted

    # define conversion factor for heat capacities:
    cv_conversion_factor = 2.0 / (sc.N_A * sc.k)

    # plot everything that has a special label:
    for special in specials:
        plotted.append(special.nist_key)
        sd = d[special.nist_key]
        axes.plot(
            sd["temperature"],
            sd["cv"] * cv_conversion_factor,
            label=special.nist_key,
            color=special.line_color,
            linestyle=special.line_style,
            alpha=0.75,
        )

    # plot everything that is not being deliberately excluded that
    # has not been plotted yet:
    for k in list(d.keys()):
        sd = d[k]
        if (k not in dont_plot) and (k not in plotted):
            axes.plot(
                sd["temperature"],
                sd["cv"] * cv_conversion_factor,
                label=k,
                alpha=0.25,
                color="k",
            )

    # place labels for species:
    def place(special):
        sd = d[special.nist_key]
        axes = plt.gca()
        if (special.position[0] != 0) and (special.position[1] != 0):
            ann = axes.annotate(
                special.label,
                xy=special.position,
                xytext=special.offset,
                textcoords="offset points",
                bbox=None,
                arrowprops=dict(
                    arrowstyle="->",
                    connectionstyle="arc3",
                    shrinkB=0,
                    edgecolor="black",
                ),
            )
            ann.set_verticalalignment(special.va)
            ann.set_horizontalalignment(special.ha)

    for special in specials:
        place(special)

    # label room temperature:
    axes.vlines(300, 0.825 * 2, 8, alpha=0.25, color="k", linestyle="dashed")
    axes.text(
        300,
        0.15 * 2,
        "room\n  temperature",
        horizontalalignment="center",
        verticalalignment="bottom",
        multialignment="center",
    )

if __name__ == "__main__":

    style_file = os.path.splitext(__file__)[0] + ".mplstyle"
    output_filename = os.path.splitext(__file__)[0] + "_generated.pdf"

    tries = 0
    max_tries = 2
    while (tries < max_tries):
        tries += 1
        try:
            with plt.style.context(style_file):
                make_plot()
                # use output filename based on script filename:
                plt.savefig(output_filename)
            break
        except FileNotFoundError as e:
            print("downloading data")
            download_data()
            print("creating pickle file")
            create_pickle_file()
    else:
        print(f"Failed after: {max_tries} tries")
