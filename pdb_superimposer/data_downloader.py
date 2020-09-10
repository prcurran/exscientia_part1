"""

A module for key PDB data collection tools

"""
import json
import os
from urllib.parse import quote_plus
from urllib.request import Request, urlopen


def pdb_search_query(query):
    """
    Wrapper around the `RCSB PDB Search API  <https://search.rcsb.org/index.html#search-api/>`_.

    :param query: Query as a json str
    :type query: str
    :return: Search result
    :rtype: dict
    """

    # encode the query
    url = f'http://search.rcsb.org/rcsbsearch/v1/query?json={quote_plus(query)}'
    # call the url
    response = urlopen(Request(url))

    if response.status == 200:
        # response received in bytes
        decoded = response.read().decode("utf-8")
        # covert json str to dict
        return json.loads(decoded)
    else:
        print(f"Status:{response.status}\n{response.msg}")


def ftp_download(pdb_code, out_dir="pdb"):
    """
    Wrapper around `RCSB PDB FTP <https://www.rcsb.org/pdb/static.do?p=download/ftp/ftp_site_layout.html/>`_.

    :param pdb_code: PDB code for desired file
    :type pdb_code: str
    :param out_dir: Path to output directory
    :type out_dir: str
    :return: Returns PDB code if download failed
    :rtype: None or str
    """

    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    response = urlopen(Request(url))
    f = response.read().decode("utf-8")
    # write out decoded file
    with open(os.path.join(out_dir, f"{pdb_code}.pdb"), "w") as w:
        w.write(f)

