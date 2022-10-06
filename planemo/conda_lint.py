"""Logic for linting conda recipes."""


from functools import wraps
from typing import (
    Callable,
    Dict,
    List,
    Tuple,
    TYPE_CHECKING,
    Union,
)

from galaxy.tool_util.deps.conda_compat import raw_metadata
from galaxy.tool_util.lint import LintContext
from galaxy.util import unicodify

from planemo.conda_verify.recipe import (
    check_build_number,
    check_dir_content,
    check_license_family,
    check_name,
    check_requirements,
    check_source,
    check_url,
    check_version,
    FIELDS,
    get_field,
    RecipeError,
)
from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_OK,
)
from planemo.io import (
    coalesce_return_codes,
    find_matching_directories,
    info,
)
from planemo.lint import (
    handle_lint_complete,
    setup_lint,
)

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext


def lint_recipes_on_paths(ctx: "PlanemoCliContext", paths: Tuple[str], **kwds) -> int:
    """Apply conda linting procedure to recipes on supplied paths."""
    assert_tools = kwds.get("assert_recipes", True)
    recursive = kwds.get("recursive", False)
    exit_codes = []
    for recipe_dir in yield_recipes_on_paths(ctx, paths, recursive):
        if lint_conda_recipe(ctx, recipe_dir, **kwds) != 0:
            exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
        else:
            exit_codes.append(EXIT_CODE_OK)
    return coalesce_return_codes(exit_codes, assert_at_least_one=assert_tools)


def lint_conda_recipe(ctx: "PlanemoCliContext", recipe_dir: str, **kwds) -> int:
    info("Linting conda recipe %s" % recipe_dir)
    lint_args, lint_ctx = setup_lint(ctx, **kwds)

    def apply(f):
        lint_ctx.lint(f.__name__, f, recipe_dir)

    apply(lint_name)
    apply(lint_version)
    apply(lint_summary)
    apply(lint_build_number)
    apply(lint_directory_content)
    apply(lint_license_family)
    apply(lint_about_urls)
    apply(lint_source)
    apply(lint_fields)
    apply(lint_requirements)

    return handle_lint_complete(lint_ctx, lint_args)


def wraps_recipe_error(is_error: bool = True) -> Callable:
    def outer_wrapper(f):
        @wraps(f)
        def wrapper(recipe_dir, lint_ctx):
            try:
                f(recipe_dir, lint_ctx)
            except RecipeError as e:
                if is_error:
                    lint_ctx.error(unicodify(e))
                else:
                    lint_ctx.warn(unicodify(e))
            except TypeError as e:  # Errors in recipe checking code from YAML.
                lint_ctx.error(unicodify(e))

        return wrapper

    return outer_wrapper


def lints_metadata(f: Callable) -> Callable:
    @wraps(f)
    def wrapper(recipe_dir, lint_ctx):
        meta = raw_metadata(recipe_dir)
        f(meta, lint_ctx)

    return wrapper


@wraps_recipe_error(is_error=False)
def lint_directory_content(recipe_dir: str, lint_ctx: LintContext) -> None:
    check_dir_content(recipe_dir)
    lint_ctx.info("Directory content seems okay.")


@lints_metadata
@wraps_recipe_error(is_error=False)
def lint_license_family(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    check_license_family(meta)
    lint_ctx.info("License from vaild license family.")


@lints_metadata
def lint_summary(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    summary = get_field(meta, "about/summary")

    if not summary:
        lint_ctx.warn("No summary supplied in about metadata.")

    if summary and len(summary) > 80:
        msg = "summary exceeds 80 characters"
        lint_ctx.warn(msg)


@lints_metadata
@wraps_recipe_error(is_error=False)
def lint_about_urls(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    for field in ("about/home", "about/dev_url", "about/doc_url", "about/license_url"):
        url = get_field(meta, field)
        if url:
            check_url(url)
    lint_ctx.info("About urls (if present) are valid")


@lints_metadata
@wraps_recipe_error(is_error=True)
def lint_source(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    check_source(meta)
    lint_ctx.info("Source (if present) is valid")


@lints_metadata
@wraps_recipe_error(is_error=True)
def lint_build_number(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    build_number = get_field(meta, "build/number", 0)
    check_build_number(build_number)
    lint_ctx.info("Valid build number [%s]" % build_number)


@lints_metadata
@wraps_recipe_error(is_error=True)
def lint_version(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    version = get_field(meta, "package/version")
    check_version(version)
    lint_ctx.info("Valid version number [%s]" % version)


@lints_metadata
@wraps_recipe_error(is_error=True)
def lint_name(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    name = get_field(meta, "package/name")
    check_name(name)
    lint_ctx.info("Valid recipe name [%s]" % name)


@lints_metadata
@wraps_recipe_error(is_error=False)
def lint_fields(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    # Taken from validate_meta
    for section in meta:
        if section not in FIELDS:
            raise RecipeError("Unknown section: %s" % section)
        submeta = meta.get(section)
        if submeta is None:
            submeta = {}
        for key in submeta:
            # Next two lines added for planemo since we don't do the
            # select lines thing.
            if key == "skip":
                continue

            if key not in FIELDS[section]:
                raise RecipeError("in section %r: unknown key %r" % (section, key))


@lints_metadata
@wraps_recipe_error(is_error=False)
def lint_requirements(
    meta: Dict[str, Union[Dict[str, str], Dict[str, List[str]], Dict[str, Union[str, int]], Dict[str, int], str]],
    lint_ctx: LintContext,
) -> None:
    check_requirements(meta)
    lint_ctx.info("Reference recipe files appear valid")


def yield_recipes_on_paths(ctx: "PlanemoCliContext", paths: Tuple[str], recursive: bool) -> None:
    for path in paths:
        recipe_dirs = find_matching_directories(path, "meta.yaml", recursive=recursive)
        yield from recipe_dirs


__all__ = ("lint_recipes_on_paths",)
