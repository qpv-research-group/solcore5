# This code is heavily based on/copied from the refractiveindex.info-sqlite code
# found here on GitHub: https://github.com/HugoGuillen/refractiveindex.info-sqlite
# However, there are some modifications from the original code

import os
import sqlite3
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional
import logging

import pandas as pd
import xarray as xr

from .. import SOLCORE_ROOT, config
from .nk_database import NK, MaterialNKDatabaseError
from .refractiveindex_info_DB import dboperations as DB


@NK.register_database(name="refractiveindex.info")
def get_refractiveindex_info_nk_data(name: str, pageid: Optional[int] = None, **kwrgs):
    """Gets the nk data from the refractiveindex.info database.

    If pageid is None, then the first record found for the material is provided, logging
    a warning about it.

    The support for this db is very limited, for now. In particular, the existence of
    several options for a given material and the definition of the alloys, appearing
    as the combination of binaries rather than a ternary - eg. AlGaAs appears as
    AlAs-GaAs - makes its integration with the parameteric materials or other nk
    databases pretty cumbersome. In these cases, it might be simpler to download the
    data from refractiveindex.info and loaded into a DataArray.

    Raises:
        MaterialNKDatabaseError if there's no records for that material or the pageid
        provided does not correspond to the given material.

    Args:
        name: Name of the material.
        pageid: id of the page where to find the material.

    Returns:
        DataArray with the complex refractive index as function of wavelength, in m.
    """
    logging.getLogger().warning(
        "The integration of the refractiveindex.info database "
        "is experimental and results should be taken with care."
    )
    available = search_db(name, exact=True)
    if available.empty:
        msg = f"No NK data in 'refractiveindex.info database for material '{name}'. "
        raise MaterialNKDatabaseError(msg)

    if pageid is None:
        record = available.iloc[0]
        logging.getLogger().warning(
            f"No 'pageid' provided for material '{name}'. "
            f"Using first record found:\n{record}"
        )
        pageid = record.pageid

    elif pageid not in available.pageid:
        msg = (
            f"pageid={pageid} does not correspond to a record for material '{name}'. "
            f"Available records are:\n{available}"
        )
        raise MaterialNKDatabaseError(msg)

    try:
        wl, n = nkdb_load_n(pageid)
        _, k = nkdb_load_n(pageid)
        wl *= 1e-6
    except Exception as err:
        raise MaterialNKDatabaseError(err)

    return xr.DataArray(n + 1.0j * k, dims=["wavelength"], coords={"wavelength": wl})


def download_db(url=None, interpolation_points=200, confirm=False):
    """
    This function downloads the refractiveindex.info database and creates on SQLite database at
    the path specified in the user config file.

    :param url: URL from which the zip archive of the database will be downloaded. Default is "https://refractiveindex.info/download/database/rii-database-2017-09-05.zip"
    :param interpolation_points: how many interpolation points to save the data at. Default is 200.
    :param confirm: if True, will not ask if you want to download database again even if it has been downloaded previously
    :return:
    """
    if "nk" not in config.sources("Others"):
        NK_PATH = os.path.join(config.user_folder, "nk", "nk.db")
        if not os.path.exists(NK_PATH):
            os.makedirs(Path(NK_PATH).parent)
        config["Others", "nk"] = NK_PATH
    else:
        NK_PATH = config["Others", "nk"]

    if os.path.isfile(NK_PATH) and not confirm:
        response = input(
            "There is already a downloaded database file. "
            "Do you want to download it again (Y/n)? "
        )

        if response not in "Yy":
            return
    else:
        confirm = True

    if confirm:
        db = DB.Database(NK_PATH)
        with TemporaryDirectory() as output:
            if url is None:
                db.create_database_from_url(interpolation_points, outputfolder=output)
            else:
                db.create_database_from_url(
                    interpolation_points, url, outputfolder=output
                )


def search_db(term="", exact=False, print_list: bool = False):
    """
    Search the downloaded SQLite database.

    :param term: search term (e.g. the name of a material, or a source).
    :param exact: search by exact (True) or approximate (False, default) terms.
    :param print_list: If the results should be printed in the terminal.
    :return: A list of tuples of with one tuple per database entry matching the
    search term.
    The first entry of each tuple is the pageid of the database entry.
    """

    NK_PATH = config["Others", "nk"]

    db = DB.Database(NK_PATH)
    conn = sqlite3.connect(db.db_path)
    c = conn.cursor()
    if not exact:
        c.execute(
            "SELECT * FROM pages WHERE shelf like ? or book like ? or page like ? or filepath like ?",
            ["%" + term + "%" for i in range(4)],
        )
    else:
        c.execute(
            "SELECT * FROM pages WHERE shelf like ? or book like ? or page like ? or filepath like ?",
            [term for i in range(4)],
        )

    results = pd.DataFrame(c.fetchall(), columns=db._get_pages_columns())
    conn.close()

    if print_list:
        print(results)

    return results


def nkdb_load_n(pageid):
    NK_PATH = os.path.abspath(
        config["Others"]["nk"].replace("SOLCORE_ROOT", SOLCORE_ROOT)
    )

    db = DB.Database(NK_PATH)
    res = db.get_material_n_numpy(int(pageid))
    wl = res[:, 0]
    n = res[:, 1]
    return wl, n


def nkdb_load_k(pageid):
    NK_PATH = os.path.abspath(
        config["Others"]["nk"].replace("SOLCORE_ROOT", SOLCORE_ROOT)
    )

    db = DB.Database(NK_PATH)
    res = db.get_material_k_numpy(int(pageid))
    wl = res[:, 0]
    k = res[:, 1]
    return wl, k


def create_nk_txt(pageid, file, folder=""):
    """
    This function creates two files called [file]_n.txt and [file]_k.txt (with [file]
    as specified in the arguments. The format matches that used by create_new_material,
    with the first column being the wavelength in metres and the second column the n or k values.

    :param pageid: pageid (number) of the database entry to be used
    :param file: name of the file to be created: the n and k data will be saved separately as
    [file]_n.txt and [file]_k.txt
    :param folder: folder where the files should be saved
    :return: parameter_source: file with list of other parameters for the new material
    """

    NK_PATH = os.path.abspath(
        config["Others"]["nk"].replace("SOLCORE_ROOT", SOLCORE_ROOT)
    )

    db = DB.Database(NK_PATH)
    if not os.path.exists(folder) and folder != "":
        os.makedirs(folder)
    res = db.get_material_txt(int(pageid), output=file, folder=folder)
