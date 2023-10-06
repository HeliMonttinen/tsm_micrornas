"""
Scripts to download data from Ensembl REST API.

"""

import requests
import sys


def request_tool(species, seq_region,
                 target_phylum):

    server = "https://rest.ensembl.org"
    ext = "/alignment/region/"+species+"/" + seq_region\
            + "?species_set_group=" + target_phylum +\
            ';method_link_type=EPO_EXTENDED'

    headers = {"Content-Type": "application/json"}
    r = requests.get(server+ext, headers=headers)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    return decoded[0]

