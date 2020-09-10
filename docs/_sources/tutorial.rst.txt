************
Introduction
************

==========================================
Objective
==========================================

The aim of this task was to retrieve and superimpose PDB structures using
open source tools and services and to demonstrate the functionality using CDK2,
UniProt: P24941 as an example.


==========================================
Solution
==========================================

The solution to this problem required two parts:

    #. Search and retrieve data from the PDB. This was implemented by using urllib from the Python standard library to wrap around the RCSB PDB's webservices.
    #. Superimpose the retrieved data. This was implemented using BioPython's PDB parsing and Superimposition functionality

Although this is a small demo task, I have tried to demonstrate coding practices that I
would employ in larger projects:

    * I've used a "tests first" approach.
    * Setup continuous integration with CircleCI.
    * Configured documentation with Sphinx


******************
Test the Code
******************

1. Clone this repository

.. code-block:: shell

    git clone https://github.com/prcurran/pdb_superimposer.git


2. Create a new environment

.. code-block:: shell

    conda env create -f environment.yml


3. Install this package

.. code-block:: shell

    conda activate super
    conda env create -f environment.yml
    pip install .


4. Run the CDK2 example

.. code-block:: shell

    python example.py


5. Remove environment and clean up

.. code-block:: shell

    conda env remove --name super --all
    rm -rf pdb_superimposer


******************
CDK2 Example
******************

This section will run through the code for the CDK2 example

=========================
Step 1: Search the PDB
=========================

The first step is to search the PDB. This is done using the RCSB Search API which is a RESTful API
service that uses URL encoded JSON payloads over HTTP.

The full documentation for the query syntax can be found `here <http://search.rcsb.org/>`_.

For this example, we use 3 of the request building blocks:
    * `query` : This is where the search expression is constructed. For this example only the uniprot code is used but others could have been included (such as structure resolution).
    * `request_options` : This controls various aspects of the search request, the pagination was alter to ensure all results were returned (default = first 10).
    * `return_type`: Specifics the type of returned identifiers, the polymer entity was used over entity since each entity can contain multiple polymers.

.. code-block:: JSON

    {
      "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "operator": "exact_match",
              "value": "P24941",
              "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
            }
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "operator": "exact_match",
              "value": "UniProt",
              "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
            }
          }
        ]
      },
      "request_options": {
        "pager": {
          "start": 0,
          "rows": 500
        }
      },
      "return_type": "polymer_entity"
    }


The query was read from file and `urllib <https://docs.python.org/3/library/urllib.html/>`_. was used to send the request

.. code-block:: python

    search_query_path = "example_search.json"

    with open(search_query_path, "r") as r:
        query = r.read()  # keep query as str

    results = pdb_search_query(query)

    pdb_entities = {i["identifier"].split("_")[0]: int(i["identifier"].split("_")[1]) - 1
                    for i in results["result_set"]}


=================================
Step 2: Download the Search Hits
=================================

Next, the data is download. This is done using the `RCSB FTP service <https://www.rcsb.org/pdb/static.do?p=download/ftp/ftp_site_layout.html/>`_.
Since for CDK2 there are over 400 entries, `multiprocessing <https://docs.python.org/3/library/multiprocessing.html/>`_ has been implemented to download
the files in parallel. I've also incorported a `tqdm <https://github.com/tqdm/tqdm/>`_
progress bar so that the user can see something is happening.

.. code-block:: python

    def wrap_ftp_download(inputs):
    #  simple wrapper to manage flow of args
    pdb, out_dir = inputs
    return ftp_download(pdb, out_dir)

    args = ((a, out_dir) for a in pdb_entities.keys())
    with Pool(processes=processes) as pool:
        list(tqdm(pool.imap_unordered(wrap_ftp_download, args), total=len(pdb_entities)))


=================================
Step 3: The Superimposition
=================================

Finally, the downloaded files can be superimposed. `BioPython's <https://biopython.org/>`_
:class:`Bio.PDB.Superimposer.Superimposer` contains
the functionality to do the actual transformation, minimising the RMS in the solution. However,
the selection of the atoms to be considered in the superimposition had to be implemented in this
package. Each residue in a PDB file contains a sequence identifier which corresponds to that residues
position when aligned with the UniProt reference sequence. I used these identifier to ensure:

    * In both the `reference` and `other` chain there are residues at a given index (some are missing, some are expression tags)
    * At a given index, the residue in `reference` and `other` is the same (some are mutated)
    * For a given residue, all atoms are present (some are partially model)
    * For a given residue, only heavy atoms are considered

I also included functionality to only select binding site residues since this is a useful operation
in some use case. As this task required ALL CDK2 structures, and some structures don't have bound ligands
it wasn't used here.

.. code-block:: python

    reference = Helper.protein_from_file(ref_id, os.path.join(out_dir, f"{ref_id}.pdb"))
    ref_chain = [c for c in reference[0]][polymer_entity]

    for pdb, entity in pdb_entities.items():
        other = Helper.protein_from_file(pdb, os.path.join(out_dir, f"{pdb}.pdb"))
        other_chain = [c for c in other[0]][polymer_entity]

        cs = ChainSuperimposer(reference=ref_chain, other=other_chain, other_struc=other)
        cs.superimpose()

        Helper.protein_to_file(other, <out_path>)
