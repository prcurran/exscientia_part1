import os
from multiprocessing import Pool

from tqdm import tqdm

from pdb_superimposer.data_downloader import pdb_search_query, ftp_download
from pdb_superimposer.superimposition import ChainSuperimposer, Helper


def wrap_ftp_download(inputs):
    #  simple wrapper to manage flow of args
    pdb, out_dir = inputs
    return ftp_download(pdb, out_dir)


def main():
    ########################################################
    search_query_path = "data/uniprot_search_query.json"
    out_dir = "pdb"
    ref_id = "2VTA"
    polymer_entity = 0  # PDB speak for chain
    ########################################################

    # Step 1: Run the UniProt Search
    with open(search_query_path, "r") as r:
        query = r.read()  # keep query as str

    results = pdb_search_query(query)

    pdb_entities = {i["identifier"].split("_")[0]: int(i["identifier"].split("_")[1]) - 1
                    for i in results["result_set"][:50]}

    # Step 2: Download PDBs
    args = ((a, out_dir) for a in pdb_entities.keys())
    with Pool(processes=3) as pool:
        fails = [r for r in tqdm(pool.imap_unordered(wrap_ftp_download, args), total=len(pdb_entities)) if r]
        for f in fails:
            pdb_entities.pop(f)

    # Step 3: Overlay
    reference = Helper.protein_from_file(ref_id, os.path.join(out_dir, f"{ref_id}.pdb"))
    ref_chain = [c for c in reference[0]][polymer_entity]  # always take the first model

    for pdb, entity in pdb_entities.items():
        other = Helper.protein_from_file(pdb, os.path.join(out_dir, f"{pdb}.pdb"))
        # always take the first model, the polymer_entity code is the code that matched in the UniProt search
        other_chain = [c for c in other[0]][polymer_entity]

        cs = ChainSuperimposer(reference=ref_chain, other=other_chain, other_struc=other)
        cs.superimpose()

        Helper.protein_to_file(other, os.path.join(out_dir, f"aligned_{pdb}.pdb"))


if __name__ == "__main__":
    main()

