from unittest import mock

from planemo.galaxy.workflows import (
    _elements_to_test_def,
    _job_inputs_template_from_invocation,
)

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


SLASH_HDA_ELEMENT = dict(HDA_ELEMENT, element_identifier="Video/Audio File")


def test_hda_with_slash_in_identifier_to_input_test_def():
    element_def = _elements_to_test_def(
        elements=[SLASH_HDA_ELEMENT],
        test_data_base_path="test-data/label",
        download_function=lambda *args, **kwargs: None,
    )
    assert element_def == [
        {
            "class": "File",
            "identifier": "Video/Audio File",
            "path": "test-data/label_Video_Audio File.fastqsanger",
        }
    ]


def test_hda_with_slash_in_identifier_to_output_test_def():
    element_def = _elements_to_test_def(
        elements=[SLASH_HDA_ELEMENT],
        test_data_base_path="test-data/label",
        download_function=lambda *args, **kwargs: None,
        definition_style="output",
    )
    assert element_def == {"Video/Audio File": {"path": "test-data/label_Video_Audio File.fastqsanger"}}


def _hda_element(identifier, dataset_id):
    return dict(HDA_ELEMENT, element_identifier=identifier, object=dict(HDA_ELEMENT["object"], id=dataset_id))


def test_labels_sanitizing_to_the_same_filename_get_distinct_paths():
    downloads = []
    element_def = _elements_to_test_def(
        elements=[_hda_element("a/b", "first"), _hda_element("a_b", "second"), _hda_element("a:b", "third")],
        test_data_base_path="test-data/label",
        download_function=lambda dataset_id, **kwargs: downloads.append((dataset_id, kwargs["file_path"])),
    )
    paths = [element["path"] for element in element_def]
    assert paths == [
        "test-data/label_a_b.fastqsanger",
        "test-data/label_a_b_1.fastqsanger",
        "test-data/label_a_b_2.fastqsanger",
    ]
    assert downloads == list(zip(["first", "second", "third"], paths))


def test_same_element_identifier_in_nested_collections_gets_distinct_paths():
    collections = [dict(COLLECTION_ELEMENT, element_identifier=name) for name in ("sample1", "sample2")]
    element_def = _elements_to_test_def(
        elements=collections,
        test_data_base_path="test-data/label",
        download_function=lambda *args, **kwargs: None,
    )
    paths = [nested["elements"][0]["path"] for nested in element_def]
    assert paths == ["test-data/label_forward.fastqsanger", "test-data/label_forward_1.fastqsanger"]


def test_disambiguated_path_does_not_take_a_path_that_is_already_used():
    element_def = _elements_to_test_def(
        elements=[_hda_element("a/b", "1"), _hda_element("a_b_1", "2"), _hda_element("a_b", "3")],
        test_data_base_path="test-data/label",
        download_function=lambda *args, **kwargs: None,
    )
    paths = [element["path"] for element in element_def]
    assert len(set(paths)) == 3


def test_input_datasets_sanitizing_to_the_same_filename_are_not_overwritten():
    user_gi = mock.MagicMock()
    user_gi.invocations.show_invocation.return_value = {
        "inputs": {
            "0": {"label": "a/b", "src": "hda", "id": "first"},
            "1": {"label": "a_b", "src": "hda", "id": "second"},
        },
        "input_step_parameters": {},
    }
    user_gi.datasets.show_dataset.return_value = {"extension": "txt"}
    with mock.patch("planemo.galaxy.workflows.gi", return_value=user_gi):
        template = _job_inputs_template_from_invocation("invocation_id", "http://galaxy.example", "key")
    assert template["a/b"]["path"] == "test-data/a_b.txt"
    assert template["a_b"]["path"] == "test-data/a_b_1.txt"
    downloads = {call.args[0]: call.kwargs["file_path"] for call in user_gi.datasets.download_dataset.call_args_list}
    assert downloads == {"first": "test-data/a_b.txt", "second": "test-data/a_b_1.txt"}
