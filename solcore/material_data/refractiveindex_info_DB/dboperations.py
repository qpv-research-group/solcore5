# From: https: // github.com / HugoGuillen / refractiveindex.info - sqlite

from collections import namedtuple,OrderedDict
import os
import yaml
import sqlite3
import numpy as np

# from solcore.material_data.refractiveindex_info_DB import dbmaterial
from solcore.material_data.refractiveindex_info_DB import DBMaterial

Shelf = namedtuple('Shelf', ['shelf', 'name'])
Book = namedtuple('Book', ['book', 'name'])
Page = namedtuple('Page', ['page', 'name', 'path'])
Entry = namedtuple('Entry',['id','shelf','book','page'])

# Check latest available database in https://refractiveindex.info/download.php
_riiurl = "https://refractiveindex.info/download/database/rii-database-2018-07-01.zip"

class Database:

    def __init__(self, sqlitedbpath):
        self.db_path = sqlitedbpath
        if not os.path.isfile(sqlitedbpath):
            print("Database file not found.")
        else:
            print("Database file found at",sqlitedbpath)

    def create_database_from_folder(self, yml_database_path, interpolation_points=100):
        create_sqlite_database(yml_database_path, self.db_path,interpolation_points=interpolation_points)

    def create_database_from_url(self,interpolation_points=100,riiurl=_riiurl, outputfolder=""):
        Database.DownloadRIIzip(riiurl=riiurl, outputfolder=outputfolder)
        outputfolder = os.path.join(outputfolder, "database")
        self.create_database_from_folder(outputfolder, interpolation_points=interpolation_points)
        pass

    def check_url_version(self):
        print(_riiurl)

    def search_custom(self,sqlquery):
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute(sqlquery)
        results = c.fetchall()
        if len(results)==0:
            print("No results found.")
        else:
            print(len(results),"results found.")
        conn.close()
        return results

    def search_pages(self,term="",exact=False):
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        if not exact:
            c.execute('SELECT * FROM pages WHERE shelf like ? or book like ? or page like ? or filepath like ?', ["%"+term+"%" for i in range(4)])
        else:
            c.execute('SELECT * FROM pages WHERE shelf like ? or book like ? or page like ? or filepath like ?', [term for i in range(4)])
        results = c.fetchall()
        if len(results)==0:
            print("No results found.")
        else:
            print(len(results),"results found.")
            columns = self._get_pages_columns()
            print("\t".join(columns))
            for r in results:
                print("\t".join(map(str,r[:])))
        conn.close()
        #return results

    def search_id(self,pageid):
        info = self._get_page_info(pageid)
        if info is None:
            print("PageID not found.")
        else:
            print("\t".join(info.keys()))
            print("\t".join(map(str,info.values())))

    def search_n(self,n,delta_n):
        print("*Search n =",n,"delta_n =",delta_n)
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        interval = [n-delta_n,n+delta_n]
        c.execute('''select r.pageid,shelf,book,page,r.wave,r.refindex
                    from refractiveindex r join pages p on r.pageid = p.pageid
                    where refindex between ? and ? ''',interval)
        results = c.fetchall()
        if len(results)==0:
            print("No results found.")
        else:
            print(len(results),"results found.")
            print("pageid|shelf|book|page|wavelength|n")
            for r in results:
                print(r)
        conn.close()

    def search_k(self,k,delta_k):
        print("*Search k =",k,"delta_k =",delta_k)
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        interval = [k-delta_k,k+delta_k]
        c.execute('''select e.pageid,shelf,book,page,e.wave,e.coeff
                    from extcoeff e join pages p on e.pageid = p.pageid
                    where coeff between ? and ?''',interval)
        results = c.fetchall()
        if len(results)==0:
            print("No results found.")
        else:
            print(len(results),"results found.")
            print("pageid|shelf|book|page|wavelength|k")
            for r in results:
                print(r)
        conn.close()

    def search_nk(self,n,delta_n,k,delta_k):
        print("*Search n =",n,"delta_n =",delta_n,"k =",k,"delta_k =",delta_k)
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        interval = [n-delta_n,n+delta_n,k-delta_k,k+delta_k]
        c.execute('''select r.pageid, shelf, book, page, r.wave, r.refindex, e.coeff
                    from refractiveindex r join extcoeff e on r.pageid = e.pageid and r.wave = e.wave
                    join pages p on r.pageid = p.pageid
                    where refindex between ? and ? and coeff between ? and ?''',interval)
        results = c.fetchall()
        if len(results)==0:
            print("No results found.")
        else:
            print(len(results),"results found.")
            print("pageid|shelf|book|page|wavelength|n|k")
            for r in results:
                print(r)
        conn.close()

    def get_material(self, pageid):
        pagedata = self._get_page_info(pageid)
        if pagedata is None:
            print("PageID not found.")
            return None
        else:
            conn = sqlite3.connect(self.db_path)
            c = conn.cursor()
            wavelengths_r = None
            wavelengths_e = None
            refractive = None
            extinction = None
            if pagedata['hasrefractive'] == 1:
                c.execute('''select wave,refindex
                            from refractiveindex
                            where pageid = ?
                            order by wave asc''', [pageid])
                results = c.fetchall()
                wavelengths_r = [r[0] for r in results]
                refractive = [r[1] for r in results]
            if pagedata['hasextinction'] == 1:
                c.execute('''select wave,coeff
                            from extcoeff
                            where pageid = ?
                            order by wave asc''', [pageid])
                results = c.fetchall()
                wavelengths_e = [r[0] for r in results]
                extinction = [r[1] for r in results]
            conn.close()
            print("Material",pagedata['filepath'],"loaded.")
            return DBMaterial.FromLists(pagedata,wavelengths_r=wavelengths_r,refractive=refractive,
                                      wavelengths_e=wavelengths_e,extinction=extinction)

    def get_material_n_numpy(self,pageid):
        mat = self.get_material(pageid)
        if mat is None:
            return None
        n = mat.get_complete_refractive()
        if n is None:
            print("Material has no refractive data.")
            return None
        return np.array(n)

    def get_material_k_numpy(self,pageid):
        mat = self.get_material(pageid)
        if mat is None:
            return None
        k = mat.get_complete_extinction()
        if k is None:
            print("Material has no extinction data.")
            return None
        return np.array(k)

    def get_material_csv(self, pageid, output="", folder=""):
        mat = self.get_material(pageid)
        if mat is None:
            print("PageID not found.")
            return None
        matInfo = mat.get_page_info()
        #print(matInfo)
        if output=="":
            output = ",".join([str(matInfo['pageid']),matInfo['shelf'],matInfo['book'],matInfo['page']])+".csv"
        if folder != "":
            output = folder+os.sep+output
        mat.to_csv(output)

    def get_material_txt(self, pageid, output="", folder=""):
        mat = self.get_material(pageid)
        if mat is None:
            print("PageID not found.")
            return None
        matInfo = mat.get_page_info()
        #print(matInfo)
        if output=="":
            output = "_".join([str(matInfo['pageid']),matInfo['shelf'],matInfo['book'],matInfo['page']])+".txt"
        if folder != "":
            output = folder+os.sep+output

        if '.txt' not in output:
            output = output+'.txt'

        mat.to_txt(output)

    def get_material_csv_all(self, outputfolder):
        allids = self._get_all_pageids()
        for id in allids:
            print("Processing",id)
            self.get_material_csv(pageid=id, output="", folder=outputfolder)

    def _get_pages_columns(self):
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('PRAGMA table_info(pages);')
        results = c.fetchall()
        names = [r[1] for r in results]
        conn.close()
        return names

    def _get_page_info(self,pageid):
        columns = self._get_pages_columns()
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('SELECT * FROM pages WHERE pageid = ?',[pageid])
        results = c.fetchall()
        if len(results) == 0:
            conn.close()
            return None
        else:
            row = results[0]
            data = OrderedDict.fromkeys(columns)
            for idx,c in enumerate(columns):
                data[c] = row[idx]
            #data = {columns[i]:row[i] for i in range(len(columns))}
            conn.close()
            return data

    def _get_all_pageids(self):
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('SELECT pageid FROM pages')
        results = c.fetchall()
        if len(results) == 0:
            conn.close()
            return None
        else:
            pageids = [row[0] for row in results]
            return pageids

    @staticmethod
    def DownloadRIIzip(outputfolder="", riiurl=_riiurl):
        import requests, zipfile, io, tempfile

        with tempfile.TemporaryDirectory() as fp:
            print("Making request to",riiurl)
            r = requests.get(riiurl)
            if r.ok:
                print("Downloaded and extracting...")
                z = zipfile.ZipFile(io.BytesIO(r.content))
                z.extractall(path=outputfolder)
                print("Wrote",outputfolder+"/database","from",riiurl)
                #The destination+database is the result.
                result = True
            else:
                print("There was a problem with the request.")
                result = False
        return result

def extract_entry_list(db_path):
    entries = []
    referencePath = os.path.normpath(db_path)
    idx = 0
    with open(os.path.join(referencePath, os.path.normpath("library.yml")), "r") as f:
        catalog = yaml.safe_load(f)
    for sh in catalog:
        shelf = Shelf(sh['SHELF'], sh['name'])
        for b in sh['content']:
            if 'DIVIDER' not in b:
                book = Book(b['BOOK'],b['name'])
                for p in b['content']:
                    if 'DIVIDER' not in p:
                        page = Page(p['PAGE'],
                                    p['name'],
                                    os.path.join(os.path.join(referencePath, 'data'), os.path.normpath(p['data'])))

                        entries.append(Entry(str(idx),shelf,book,page))
                        idx+=1
    return entries

def print_pretty_entry_list(entries):
    for e in entries:
        print(",".join([e.id,e.shelf.shelf,e.book.book,e.page.page]))

def pretty_entry(entry):
    e = entry
    return ",".join([e.id,e.shelf.shelf,e.book.book,e.page.page])

def create_sqlite_database(refractiveindex_db_path,new_sqlite_db,interpolation_points=100):
    conn = sqlite3.connect(new_sqlite_db)
    c = conn.cursor()
    c.execute('''DROP TABLE IF EXISTS pages;''')
    c.execute('''DROP TABLE IF EXISTS refractiveindex;''')
    c.execute('''DROP TABLE IF EXISTS extcoeff;''')
    c.execute('''CREATE TABLE pages
    (pageid int, shelf text COLLATE NOCASE, book text COLLATE NOCASE, page text COLLATE NOCASE,
    filepath text COLLATE NOCASE,
    hasrefractive integer, hasextinction integer,
    rangeMin real, rangeMax real, points int)''')
    c.execute('''CREATE TABLE refractiveindex (pageid int, wave real, refindex real)''')
    c.execute('''CREATE TABLE extcoeff (pageid int, wave real, coeff real)''')
    conn.commit()
    conn.close()
    _populate_sqlite_database(refractiveindex_db_path,new_sqlite_db,interpolation_points=interpolation_points)

def _populate_sqlite_database(refractiveindex_db_path,new_sqlite_db,interpolation_points=100):
    entries = extract_entry_list(refractiveindex_db_path)
    conn = sqlite3.connect(new_sqlite_db)
    c = conn.cursor()
    for e in entries:
        try:
            mat = DBMaterial(filename=e.page.path,interpolation_points=interpolation_points)
            hasrefractive=0
            hasextinction=0
            if mat.has_refractive():
                refr = mat.get_complete_refractive()
                hasrefractive = 1
                values = [[e.id,r[0],r[1]] for r in refr]
                c.executemany('INSERT INTO refractiveindex VALUES (?,?,?)', values)
            if mat.has_extinction():
                ext = mat.get_complete_extinction()
                hasextinction = 1
                values = [[e.id,ex[0],ex[1]] for ex in ext]
                c.executemany('INSERT INTO extcoeff VALUES (?,?,?)', values)
            c.execute("INSERT INTO pages VALUES (?,?,?,?,?,?,?,?,?,?)",
                      [e.id,
                       e.shelf.shelf,
                       e.book.book,
                       e.page.page,
                       os.sep.join(e.page.path.split(os.sep)[-3:]),
                       hasrefractive,
                       hasextinction,
                       mat.rangeMin,
                       mat.rangeMax,
                       mat.points])
        except Exception as error:
            print("LOG:",pretty_entry(e),":",error)
    conn.commit()
    conn.close()
    print("***Wrote SQLite DB on ",new_sqlite_db)


def pipeline_test():
    #Database.DownloadRIIzip()
    db = Database("../refractive.db")

    #db.check_url_version()
    #db.create_database_from_folder(yml_database_path="database", interpolation_points=200)
    #db.create_database_from_url(interpolation_points=200)
    #db.create_database_from_url(riiurl="http://refractiveindex.info/download/database/rii-database-2015-07-05.zip")

    db.search_pages()
    #db.search_pages("otanicar")
    #db.search_pages("au",exact=True)
    #Id 327 for Formula2 test

    #print(db._get_page_info(1542))
    #db.search_id(1542)

    #mat = db.get_material(1542) #Only extinction
    #mat = db.get_material(372) #Both (formula)
    #mat = db.get_material(1) #Both (tabulated)
    #print("HasRefractive?",mat.has_refractive())
    #print("HasExtinction?",mat.has_extinction())
    #print(mat.get_complete_refractive())
    #print(mat.get_complete_extinction())
    #print(mat.get_page_info())
    #mat.to_csv(output="mat1.csv")

    #db.get_material_csv(1542,output="",folder="all")
    #db.get_material_csv_all(outputfolder="all")

    #db.search_n(n=0.3,delta_n=.001)
    #db.search_k(k=0.3,delta_k=.001)
    #db.search_nk(n=0.3, delta_n=0.1,k=0.3,delta_k=0.1)

    #print(db.search_custom('select * from pages where shelf="main" and book="Ag" and page LIKE "%k%"'))
    #print(db.search_custom('select wave,coeff from extcoeff where pageid = 1 and wave between 0.3 and 0.4'))
    #print(db.search_custom('''select p.filepath, r.wave,refindex,coeff
    #                from refractiveindex r inner join extcoeff e on r.pageid = e.pageid and r.wave = e.wave
    #                inner join pages p on r.pageid = p.pageid
    #                where r.wave = .301'''))

if __name__ == '__main__':
    pipeline_test()
