"""Lint the data tables referenced by a tool against the sample files of its repository.

The Tool Shed validates every tool of a repository revision before it defines
metadata for it. Tools whose dynamically generated select lists are populated
from a data table (``<options from_data_table="..."/>``) or from a data file
(``<options from_file="..."/>``) are only accepted together with sample files
in the repository. Otherwise uploading the revision fails with a message like::

    Unexpected HTTP status code: 500: {"err_msg": "Metadata may have been
    defined for some items in revision 'ee910cdc8874'.  Correct the following
    problems if necessary and reset metadata.<br/><b>tool.xml</b> - This file
    refers to a file named gene_ontology.loc.  Upload a file named
    gene_ontology.loc.sample to the repository to correct this error.<br/>",
    "err_code": 500001}

The linters in this module reproduce the corresponding checks of Galaxy's
``ToolValidator`` (``lib/galaxy/tool_shed/tools/tool_validator.py``), which is
what runs when metadata is set for a revision:

- a tool that uses a data table requires a ``tool_data_table_conf.xml.sample``
  file anywhere in the repository. The Tool Shed loads the tables declared by
  that file against the data path of its own Galaxy instance, where their data
  files are never present, so
- every data file of a table that is declared in that file and used by the tool
  requires a ``<index file>.sample`` file in the repository.

``ToolDataTableConfRefersToSample`` and ``ToolDataTableFromFileSampleMissing``
are warnings, not errors: the Tool Shed accepts both, but the table or the file
will usually not be found at run time. Galaxy installs repository tool data
files with the ``.sample`` suffix stripped, so a ``<file path=...>`` entry
naming ``gene_ontology.loc.sample`` keeps looking for a file that installing
the repository does not create; ``from_file`` paths are resolved in the tool
data directory, which only a data manager or an administrator populates.

Sample files are matched by base name anywhere in the repository, exactly like
the Tool Shed does. The other reasons for the Tool Shed to fail setting
metadata (an unparsable ``tool_data_table_conf.xml.sample``, broken
``tool_dependencies.xml`` / ``repository_dependencies.xml`` files, files
included by a tool that are not in the repository) are covered by the
``ToolDataTableConfSampleInvalid`` linter and by the repository linters in
:mod:`planemo.shed_lint`.
"""

import os
from dataclasses import dataclass
from typing import (
    Dict,
    FrozenSet,
    Optional,
    Tuple,
    TYPE_CHECKING,
)
from xml.etree import ElementTree

from galaxy.tool_util.lint import Linter

from planemo.linters.util import cached_lint_result

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource
    from galaxy.util import Element

DATA_TABLE_CONF_SAMPLE = "tool_data_table_conf.xml.sample"
TOOL_DATA_DIRECTORY = "tool-data"
SAMPLE_SUFFIX = ".sample"
# Markers identifying the root of a Tool Shed repository. The nearest ancestor
# carrying one of these is used, since the Tool Shed looks for sample files in
# the whole repository and not only in the directory of the tool.
REPOSITORY_MARKERS = (".shed.yml", ".shed.yaml", ".git", ".hg")
EXCLUDED_DIRECTORIES = (".git", ".hg")

MESSAGE_ERROR_CONF_SAMPLE_MISSING = (
    "No {conf} was found in the repository but this tool uses data table [{name}]. "
    "The Tool Shed will not set metadata for this revision - upload a {conf} that declares table "
    "[{name}] together with the sample file of the data file it refers to (usually "
    "{tool_data}/{name}.loc{suffix})."
)
MESSAGE_ERROR_CONF_SAMPLE_INVALID = (
    "{conf} could not be parsed ({error}). The Tool Shed will not set metadata for this revision."
)
MESSAGE_ERROR_INDEX_SAMPLE_MISSING = (
    "Data table [{name}] refers to {index_file} but the repository does not include {index_file}{suffix}. "
    "The Tool Shed will not set metadata for this revision - upload a file named {index_file}{suffix} "
    "(usually to {tool_data}/)."
)
MESSAGE_WARN_CONF_REFERS_TO_SAMPLE = (
    "The data table [{name}] in {conf} names the sample file {index_file} instead of the data file it "
    "samples. Galaxy installs repository tool data files with the {suffix} suffix stripped, so installing "
    "this repository creates {data_file} while the entry keeps looking for {index_file} - the table only "
    "has entries if a data manager or the administrator creates that file. Name the data file in the "
    "entry, for example {tool_data}/{name}.loc."
)
MESSAGE_WARN_FROM_FILE_SAMPLE_MISSING = (
    "The tool reads {index_file} with from_file but the repository does not ship {index_file}{suffix}. "
    "Galaxy resolves from_file paths in the tool data directory, so administrators need a file named "
    "{index_file} there - ship {index_file}{suffix} (usually in {tool_data}/) as a template for them to "
    "create it, or declare the file as a data table."
)

# These linters only add messages for tools that actually reference data tables.
lint_tool_types = ["*"]


@dataclass(frozen=True)
class DataTableReference:
    """A ``<options from_data_table="..."/>`` element of a tool."""

    node: "Element"
    name: str


@dataclass(frozen=True)
class IndexFileReference:
    """A ``<options from_file="..."/>`` element of a tool."""

    node: "Element"
    index_file: str


@dataclass(frozen=True)
class DataTableFile:
    """A ``<file path="..."/>`` element of a ``<table>`` in a data table config."""

    path: str
    index_file: str


@dataclass(frozen=True)
class RepositorySamples:
    """The ``*.sample`` files and the data table config shipped by a repository."""

    root: str
    conf_path: Optional[str]
    conf_error: Optional[str]
    sample_names: FrozenSet[str]
    tables: Dict[str, Tuple[DataTableFile, ...]]


@dataclass(frozen=True)
class DataTableAnalysis:
    """Data table and data file references of a tool plus the repository they live in."""

    data_tables: Tuple[DataTableReference, ...]
    index_files: Tuple[IndexFileReference, ...]
    samples: RepositorySamples


def _data_table_analysis(tool_source: "ToolSource") -> Optional[DataTableAnalysis]:
    def calculate() -> Optional[DataTableAnalysis]:
        tool_xml = getattr(tool_source, "xml_tree", None)
        source_path = tool_source.source_path
        if tool_xml is None or not source_path:
            # Without the path of the tool XML there is no repository to check against.
            return None
        data_tables = []
        index_files = []
        for options in tool_xml.findall(".//options"):
            table_name = options.get("from_data_table")
            if table_name and "$" not in table_name:
                data_tables.append(DataTableReference(options, table_name))
            index_file = options.get("from_file")
            if index_file and "$" not in index_file:
                index_files.append(IndexFileReference(options, _base_name(index_file)))
        return DataTableAnalysis(
            data_tables=tuple(data_tables),
            index_files=tuple(index_files),
            samples=_repository_samples(_repository_root(source_path)),
        )

    return cached_lint_result(tool_source, "tool_data_tables", calculate)


def _base_name(path: str) -> str:
    # The Tool Shed compares sample files by base name, see basic_util.strip_path.
    return os.path.normpath(path).replace("\\", "/").rsplit("/", 1)[-1]


def _sample_name(index_file: str) -> str:
    return f"{index_file}{SAMPLE_SUFFIX}"


def _repository_root(source_path: str) -> str:
    """Return the directory the Tool Shed would treat as repository root."""
    tool_dir = os.path.dirname(os.path.abspath(source_path))
    directory = tool_dir
    while True:
        if _is_repository_root(directory):
            return directory
        parent = os.path.dirname(directory)
        if parent == directory:
            return tool_dir
        directory = parent


def _is_repository_root(directory: str) -> bool:
    if any(os.path.exists(os.path.join(directory, marker)) for marker in REPOSITORY_MARKERS):
        return True
    return any(os.path.exists(os.path.join(directory, name)) for name in (DATA_TABLE_CONF_SAMPLE, TOOL_DATA_DIRECTORY))


def _repository_samples(root: str) -> RepositorySamples:
    sample_names = set()
    conf_path = None
    for current_root, dirs, files in os.walk(root):
        dirs[:] = [directory for directory in dirs if directory not in EXCLUDED_DIRECTORIES]
        for name in files:
            if not name.endswith(SAMPLE_SUFFIX):
                continue
            sample_names.add(name)
            # Like hg_util.get_config_from_disk the first match in the repository is used.
            if conf_path is None and name == DATA_TABLE_CONF_SAMPLE:
                conf_path = os.path.join(current_root, name)
    tables, conf_error = _parse_conf(conf_path)
    return RepositorySamples(
        root=root,
        conf_path=conf_path,
        conf_error=conf_error,
        sample_names=frozenset(sample_names),
        tables=tables,
    )


def _parse_conf(conf_path: Optional[str]) -> Tuple[Dict[str, Tuple[DataTableFile, ...]], Optional[str]]:
    """Return the data files of each ``<table>`` declared by a data table config."""
    if conf_path is None:
        return {}, None
    try:
        conf_root = ElementTree.parse(conf_path).getroot()
    except (ElementTree.ParseError, OSError) as e:
        return {}, str(e)
    tables: Dict[str, Tuple[DataTableFile, ...]] = {}
    # Galaxy reads the tables of a config file with ``getroot().findall("table")``.
    for table in conf_root.findall("table"):
        name = table.get("name")
        if not name:
            continue
        tables[name] = tuple(
            DataTableFile(path=path, index_file=_base_name(path))
            for path in (file.get("path") for file in table.findall("file"))
            if path
        )
    return tables, None


def _relative_conf_path(samples: RepositorySamples) -> str:
    if samples.conf_path is None:
        return DATA_TABLE_CONF_SAMPLE
    return os.path.relpath(samples.conf_path, samples.root)


def _data_files_without_sample(samples: RepositorySamples, data_tables: Tuple[DataTableReference, ...]):
    """Yield the data files of the tables a tool uses which have no sample file in the repository."""
    reported = set()
    for reference in data_tables:
        for data_file in samples.tables.get(reference.name, ()):
            if _sample_name(data_file.index_file) in samples.sample_names:
                continue
            if (reference.name, data_file.index_file) in reported:
                continue
            reported.add((reference.name, data_file.index_file))
            yield reference, data_file


class ToolDataTableConfSampleMissing(Linter):
    """Check that the repository declares the data tables used by the tool."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        analysis = _data_table_analysis(tool_source)
        if analysis is None or analysis.samples.conf_path is not None:
            return
        reported = set()
        for reference in analysis.data_tables:
            if reference.name in reported:
                continue
            reported.add(reference.name)
            lint_ctx.error(
                MESSAGE_ERROR_CONF_SAMPLE_MISSING.format(
                    name=reference.name,
                    conf=DATA_TABLE_CONF_SAMPLE,
                    tool_data=TOOL_DATA_DIRECTORY,
                    suffix=SAMPLE_SUFFIX,
                ),
                linter=cls.name(),
                node=reference.node,
            )


class ToolDataTableConfSampleInvalid(Linter):
    """Check that the data table config of the repository can be parsed."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        analysis = _data_table_analysis(tool_source)
        if analysis is None or not analysis.data_tables or analysis.samples.conf_error is None:
            return
        lint_ctx.error(
            MESSAGE_ERROR_CONF_SAMPLE_INVALID.format(
                conf=_relative_conf_path(analysis.samples), error=analysis.samples.conf_error
            ),
            linter=cls.name(),
            node=analysis.data_tables[0].node,
        )


class ToolDataTableIndexSampleMissing(Linter):
    """Check that the repository ships a sample file for the data files of the used tables."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        analysis = _data_table_analysis(tool_source)
        if analysis is None or analysis.samples.conf_path is None or analysis.samples.conf_error:
            return
        for reference, data_file in _data_files_without_sample(analysis.samples, analysis.data_tables):
            if data_file.index_file.endswith(SAMPLE_SUFFIX):
                continue
            lint_ctx.error(
                MESSAGE_ERROR_INDEX_SAMPLE_MISSING.format(
                    name=reference.name,
                    index_file=data_file.index_file,
                    suffix=SAMPLE_SUFFIX,
                    tool_data=TOOL_DATA_DIRECTORY,
                ),
                linter=cls.name(),
                node=reference.node,
            )


class ToolDataTableConfRefersToSample(Linter):
    """Warn about ``<file path=...>`` entries naming a sample file instead of a data file."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        analysis = _data_table_analysis(tool_source)
        if analysis is None or analysis.samples.conf_path is None or analysis.samples.conf_error:
            return
        conf = _relative_conf_path(analysis.samples)
        for reference, data_file in _data_files_without_sample(analysis.samples, analysis.data_tables):
            if not data_file.index_file.endswith(SAMPLE_SUFFIX):
                continue
            lint_ctx.warn(
                MESSAGE_WARN_CONF_REFERS_TO_SAMPLE.format(
                    name=reference.name,
                    conf=conf,
                    index_file=data_file.index_file,
                    data_file=data_file.index_file[: -len(SAMPLE_SUFFIX)],
                    suffix=SAMPLE_SUFFIX,
                    tool_data=TOOL_DATA_DIRECTORY,
                ),
                linter=cls.name(),
                node=reference.node,
            )


class ToolDataTableFromFileSampleMissing(Linter):
    """Warn about ``from_file`` data files without a sample file in the repository."""

    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        analysis = _data_table_analysis(tool_source)
        if analysis is None:
            return
        reported = set()
        for reference in analysis.index_files:
            index_file = reference.index_file
            if _sample_name(index_file) in analysis.samples.sample_names or index_file in reported:
                continue
            reported.add(index_file)
            lint_ctx.warn(
                MESSAGE_WARN_FROM_FILE_SAMPLE_MISSING.format(
                    index_file=index_file, suffix=SAMPLE_SUFFIX, tool_data=TOOL_DATA_DIRECTORY
                ),
                linter=cls.name(),
                node=reference.node,
            )
