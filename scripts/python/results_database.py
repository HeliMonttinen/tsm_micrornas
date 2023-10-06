from mongoengine import (Document, StringField,
                         ListField, IntField,
                         BooleanField)


class Results_microRNA(Document):
    """
    A class which defines a structure for a microRNA document in
    a mongo database.
    """

    identifier = StringField(required=True)
    miR_id = StringField()
    position = IntField()
    human_chromosome = StringField()
    human_location = StringField()
    query = StringField()
    parent = StringField()
    sister = StringField()
    ancestor = StringField()
    query_child1 = StringField()
    query_child2 = StringField()
    tsm_source_start = IntField()
    tsm_source_end = IntField()
    tsm_target_start = IntField()
    tsm_target_end = IntField()
    taxonomic_group = ListField()
    sequence_length = IntField()
    TSM_length = IntField()
    mismatches = IntField()
    quality = BooleanField()
    overlap = StringField()
    elements = StringField()
    genomic_context = StringField()

    meta = {
        'db_alias': 'default',
        'collection':'TSMs',
        'indexes': [
            'identifier',
            '$identifier',  # text index
            '#identifier',  # hashed index
        ]
    }


def add_information(result_dict):
    """
    A function which adds information into a database document.
    Requires a dictionary-formatted results as an input.
    """

    result_document = Results_microRNA(
            identifier=result_dict["identifier"],
            miR_id=result_dict["miR_id"],
            position=result_dict["position"],
            human_chromosome=result_dict["human_chromosome"],
            human_location=result_dict["human_location"],
            query=result_dict["query"],
            parent=result_dict["parent"],
            sister=result_dict["sister"],
            ancestor=result_dict["ancestor"],
            query_child1=result_dict["child1"],
            query_child2=result_dict["child2"],
            tsm_source_start=result_dict["tsm_orig_start"],
            tsm_source_end=result_dict["tsm_orig_end"],
            tsm_target_start=result_dict["tsm_target_start"],
            tsm_target_end=result_dict["tsm_target_end"],
            taxonomic_group=result_dict["taxonomic_group"],
            sequence_length=result_dict["sequence_length"],
            TSM_length=result_dict["TSM_length"],
            mismatches=result_dict["mismatches"],
            quality=result_dict["quality"])

    result_document.save()



def add_information_allgen(result_dict):
    """
    A function which adds information into a database document.
    Requires a dictionary-formatted results as an input.
    """

    result_document = Results_microRNA(
            identifier=result_dict["identifier"],
            position=result_dict["position"],
            human_chromosome=result_dict["human_chromosome"],
            human_location=result_dict["human_location"],
            query=result_dict["query"],
            parent=result_dict["parent"],
            sister=result_dict["sister"],
            ancestor=result_dict["ancestor"],
            query_child1=result_dict["child1"],
            query_child2=result_dict["child2"],
            tsm_source_start=result_dict["tsm_orig_start"],
            tsm_source_end=result_dict["tsm_orig_end"],
            tsm_target_start=result_dict["tsm_target_start"],
            tsm_target_end=result_dict["tsm_target_end"],
            taxonomic_group=result_dict["taxonomic_group"],
            TSM_length=result_dict["TSM_length"],
            mismatches=result_dict["mismatches"],
            quality=result_dict["quality"],
            overlap=result_dict["overlap"])

    result_document.save()
