"""Ensure requirements are matched in best practice conda channels."""

from galaxy.tools.deps.conda_util import best_search_result

from planemo.conda import tool_source_conda_targets

BEST_PRACTICE_CHANNELS = ["conda-forge", "anaconda", "r", "bioconda"]


def lint_requirements_in_conda(tool_source, lint_ctx):
    """Check requirements of tool source against best practice Conda channels."""
    conda_targets = tool_source_conda_targets(tool_source)
    for conda_target in conda_targets:
        (best_hit, exact) = best_search_result(conda_target, channels_override=BEST_PRACTICE_CHANNELS)
        conda_target_str = conda_target.package
        if conda_target.version:
            conda_target_str += "@%s" % (conda_target.version)
        if best_hit and exact:
            template = "Requirement [%s] matches target in best practice Conda channel [%s]."
            message = template % (conda_target_str, best_hit.get("channel"))
            lint_ctx.info(message)
        elif best_hit:
            template = "Requirement [%s] doesn't exactly match available version [%s] in best practice Conda channel [%s]."
            message = template % (conda_target_str, best_hit['version'], best_hit.get("channel"))
            lint_ctx.warn(message)
        else:
            template = "Requirement [%s] doesn't match any recipe in a best practice conda channel [%s]."
            message = template % (conda_target_str, BEST_PRACTICE_CHANNELS)
            lint_ctx.warn(message)
