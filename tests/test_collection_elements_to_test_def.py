from planemo.galaxy.workflows import _elements_to_test_def

HDA_ELEMENT = {
    "id": "e9d955ca403027d5",
    "model_class": "DatasetCollectionElement",
    "element_index": 0,
    "element_identifier": "forward",
    "element_type": "hda",
    "object": {
        "id": "dd871f1d42ebeefa",
        "model_class": "HistoryDatasetAssociation",
        "state": "ok",
        "hda_ldda": "hda",
        "history_id": "a2619fc82a1e31e5",
        "create_time": "2022-06-01T10:08:46.155211",
        "uuid": "72c62dcd-ba7e-4644-9ff9-7303425e391f",
        "peek": "bla",
        "hid": 3,
        "genome_build": "?",
        "file_ext": "fastqsanger",
        "metadata_dbkey": "?",
        "metadata_sequences": 2500,
        "deleted": False,
        "misc_info": "uploaded fastqsanger file",
        "history_content_type": "dataset",
        "tags": [],
        "update_time": "2022-06-01T10:54:45.515672",
        "visible": False,
        "file_size": 918886,
        "data_type": "galaxy.datatypes.sequence.FastqSanger",
        "name": "forward",
        "purged": False,
        "validated_state_message": None,
        "misc_blurb": "2,500 sequences",
        "metadata_data_lines": 10000,
        "validated_state": "unknown",
    },
}

COLLECTION_ELEMENT = {
    "id": "3c870daec78336de",
    "model_class": "DatasetCollectionElement",
    "element_index": 0,
    "element_identifier": "ERR3485802",
    "element_type": "dataset_collection",
    "object": {
        "id": "5d2260615043ce1b",
        "model_class": "DatasetCollection",
        "collection_type": "paired",
        "populated": True,
        "element_count": 1,
        "contents_url": None,
        "elements": [
            HDA_ELEMENT,
        ],
    },
}


def test_hda_to_input_test_def():
    element_def = _elements_to_test_def(
        elements=[HDA_ELEMENT], test_data_base_path="test-data/label", download_function=lambda *args, **kwargs: None
    )
    assert element_def == [{"class": "File", "identifier": "forward", "path": "test-data/label_forward.fastqsanger"}]


def test_hda_to_output_test_df():
    element_def = _elements_to_test_def(
        elements=[HDA_ELEMENT],
        test_data_base_path="test-data/label",
        download_function=lambda *args, **kwargs: None,
        definition_style="output",
    )
    assert element_def == {"forward": {"path": "test-data/label_forward.fastqsanger"}}


def test_dataset_collection_element_to_input_test_def():
    element_def = _elements_to_test_def(
        elements=[COLLECTION_ELEMENT],
        test_data_base_path="test-data/label",
        download_function=lambda *args, **kwargs: None,
    )
    assert element_def == [
        {
            "class": "Collection",
            "type": "paired",
            "identifier": "ERR3485802",
            "elements": [{"class": "File", "identifier": "forward", "path": "test-data/label_forward.fastqsanger"}],
        }
    ]


def test_dataset_collection_element_to_output_test_df():
    element_def = _elements_to_test_def(
        elements=[COLLECTION_ELEMENT],
        test_data_base_path="test-data/label",
        download_function=lambda *args, **kwargs: None,
        definition_style="output",
    )
    assert element_def == {"ERR3485802": {"elements": {"forward": {"path": "test-data/label_forward.fastqsanger"}}}}
