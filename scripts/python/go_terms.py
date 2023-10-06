"""
This script downloads the alignments using
a given list of Ensebl identifiers and locations.
"""
import mygene
from collections import defaultdict


def get_go_terms(identifier):
    """
    Get go terms for a list of Ensemble identifiers.

    Argument
    ========
    idenfier: An ensemble identifier

    Returns
    =======
    ids_dict: A dictionary containing all go identifiers in
              each go category
    ids_terms: Adictionaty containing all go terms in
               each go category for
    go_ids_name_and_term = A dictionary where each go category
                           works as a key for dictionaries where
                           go_identifiers are keays and go terms
                           are values.
    """

    ids_dict = defaultdict(dict)
    ids_terms = defaultdict(dict)
    go_ids_name_and_term = defaultdict(dict)

    mg = mygene.MyGeneInfo()
    try:
        mg_dict = mg.getgene(identifier, fields='go')
    except:
        print(error)
    if mg_dict is not None:

        for category in ["BP", "CC", "MF"]:
            go_identifier = "None"
            go_term = "None"
            ids_dict[category] = list()
            ids_terms[category] = list()
            go_ids_name_and_term[category] = dict()

            if 'go' in mg_dict and category in mg_dict["go"]:

                for item in mg_dict["go"][category]:
                    try:
                        go_identifier = item["id"]
                        go_term = item["term"]
                        ids_dict[category].append(go_identifier)
                        ids_terms[category].append(go_term)
                        go_ids_name_and_term[category][go_identifier] = go_term
                    except:
                        print("error")

            else:
                ids_dict[category].append(go_identifier)
                ids_terms[category].append(go_term)
                go_ids_name_and_term[category][go_identifier] = go_term

    return ids_dict, ids_terms, go_ids_name_and_term
