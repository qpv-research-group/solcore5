import os
import sqlite3

from solcore.material_data.refractiveindex_info_DB import dboperations as DB
from solcore import config, SOLCORE_ROOT

NK_PATH = os.path.abspath(config['Others']['nk'].replace('SOLCORE_ROOT', SOLCORE_ROOT))


def download_db(url = None, interpolation_points = 200):
    db = DB.Database(NK_PATH)
    if url is None:
        db.create_database_from_url(interpolation_points)
    else:
        db.create_database_from_url(interpolar_points, url)


def search_db(term="", exact=False):
    db = DB.Database(NK_PATH)
    conn = sqlite3.connect(db.db_path)
    c = conn.cursor()
    if not exact:
        c.execute('SELECT * FROM pages WHERE shelf like ? or book like ? or page like ? or filepath like ?',
                  ["%" + term + "%" for i in range(4)])
    else:
        c.execute('SELECT * FROM pages WHERE shelf like ? or book like ? or page like ? or filepath like ?',
                  [term for i in range(4)])
    results = c.fetchall()
    if len(results) == 0:
        print("No results found.")
    else:
        print(len(results), "results found.")
        columns = db._get_pages_columns()
        print("\t".join(columns))
        for r in results:
            print("\t".join(map(str, r[:])))
    conn.close()
    # return results


def nkdb_load_n(pageid):
    db = DB.Database(NK_PATH)
    res = db.get_material_n_numpy(int(pageid))
    wl = res[:, 0]
    n = res[:, 1]
    return wl, n


def nkdb_load_k(pageid):
    db = DB.Database(NK_PATH)
    res = db.get_material_k_numpy(int(pageid))
    wl = res[:, 0]
    k = res[:, 1]
    return wl, k