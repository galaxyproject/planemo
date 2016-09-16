""" Test app to emulate planemo-relevant portions of the
the ToolShed API... for now :).
"""
import json
import os
import tarfile

from uuid import uuid4
from xml.etree import ElementTree

from flask import (
    Flask,
    request,
    send_file,
)

app = Flask(__name__)


@app.route("/api/repositories")
def get_repositories():
    repos = app.config["model"].get_repositories()
    return json.dumps(list(repos.values()))


@app.route('/api/repositories', methods=['POST'])
def create_repository():
    repo = _request_post_message()
    # TODO: rework things to not need to hardcode this
    # i.e. simulate key stuff also.
    repo["owner"] = "iuc"
    id = str(uuid4())
    model = app.config["model"]
    model.add_repository(id, **repo)
    return json.dumps(model.get_repository(id))


@app.route('/api/repositories/<id>', methods=['PUT'])
def update_repository(id):
    repo = _request_post_message()
    repo["owner"] = "iuc"
    model = app.config["model"]
    model.update_repository(id, **repo)
    return json.dumps(model.get_repository(id))


@app.route('/api/repositories/<id>/changeset_revision', methods=['POST'])
def update_repository_contents(id):
    updated_tar = request.files['file']
    message = request.form["commit_message"]
    repo_path = app.config["model"].repository_path_for_update(id, message)
    repo_tar_path = repo_path + ".tar.gz"
    if not os.path.exists(repo_path):
        os.makedirs(repo_path)
    updated_tar.save(repo_tar_path)
    tar = tarfile.open(repo_tar_path, "r:gz")
    try:
        tar.extractall(repo_path)
    finally:
        tar.close()
    _modify_repository(repo_path)
    return json.dumps({"id": id})


@app.route("/api/categories")
def get_categories():
    categories = app.config["model"].get_categories()
    return json.dumps(categories)


@app.route("/api/users")
def get_users():
    return json.dumps([{
        "username": "iuc"
    }])


@app.route("/repository/download")
def repository_download():
    id = request.args.get("repository_id", None)
    repo_path = app.config["model"].repository_path(id)
    if not os.path.exists(repo_path):
        os.makedirs(repo_path)
    repo_tar_download_path = repo_path + "downlaod.tar.gz"
    tar = tarfile.open(repo_tar_download_path, "w:gz")
    try:
        base_name = os.path.basename(repo_path)
        tar.add(
            os.path.join(repo_path),
            arcname=base_name,
            recursive=True,
        )
    finally:
        tar.close()
    return send_file(repo_tar_download_path)


@app.route('/shutdown', methods=['POST'])
def shutdown():
    # Used to shutdown test server.
    _shutdown_server()
    return ''


def _request_post_message():
    return json.loads(request.data.decode("utf-8"))


class InMemoryShedDataModel(object):

    def __init__(self, directory):
        self.directory = directory
        self._repositories = {}
        self._repositories_msg = {}
        self._categories = []

    def add_category(self, id, name):
        self._categories.append({"id": id, "name": name})
        return self

    def add_repository(self, id, **kwds):
        repo_metadata = kwds.copy()
        repo_metadata["id"] = id
        self._repositories[id] = repo_metadata
        return self

    def update_repository(self, id, **kwds):
        repo_metadata = kwds.copy()
        repo_metadata["id"] = id
        self._repositories[id] = repo_metadata
        return self

    def get_categories(self):
        return self._categories

    def get_repositories(self):
        return self._repositories

    def get_repository(self, id):
        return self._repositories[id]

    def repository_path(self, id):
        return os.path.join(self.directory, id)

    def repository_path_for_update(self, id, message):
        if id not in self._repositories_msg:
            self._repositories_msg[id] = []
        self._repositories_msg[id].append(message)
        return self.repository_path(id)


def _modify_repository(path):
    arch = os.path.join(path, ".hg_archival.txt")
    with open(arch, "w") as f:
        f.write("Hello World!")
    deps = os.path.join(path, "tool_dependencies.xml")
    suite = os.path.join(path, "repository_dependencies.xml")
    _modify_xml(deps)
    _modify_xml(suite)


def _modify_xml(path):
    if not os.path.exists(path):
        return
    element = ElementTree.parse(path)
    _modify_attributes(element.getroot())
    as_str = ElementTree.tostring(element.getroot()).decode("UTF-8")
    with open(path, "w") as f:
        f.write(as_str)


def _modify_attributes(xml_element):
    if xml_element.tag == "repository":
        xml_element.attrib["toolshed"] = "localhost:9012"
        xml_element.attrib["changeset_revision"] = "12345"

    children = xml_element.getchildren()
    if len(children) > 0:
        for child in children:
            _modify_attributes(child)


def _shutdown_server():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()
