"""Ensure best-practice biocontainer registered for this tool."""

from galaxy.tools.deps.mulled.util import (
    build_target,
    image_name,
    mulled_tags_for,
    split_tag,
)

from planemo.conda import tool_source_conda_targets

MESSAGE_WARN_NO_REQUIREMENTS = "No valid package requirement tags found to infer BioContainer from."
MESSAGE_WARN_NO_CONTAINER = "Failed to find a BioContainer registered for these requirements."
MESSAGE_INFO_FOUND_BIOCONTAINER = "BioContainer best-practice container found [%s]."


def lint_biocontainer_registered(tool_source, lint_ctx):
    conda_targets = tool_source_conda_targets(tool_source)
    if not conda_targets:
        lint_ctx.warn(MESSAGE_WARN_NO_REQUIREMENTS)
        return

    mulled_targets = map(lambda c: build_target(c.package, c.version), conda_targets)
    name = mulled_container_name("biocontainers", mulled_targets)
    if name:
        lint_ctx.info(MESSAGE_INFO_FOUND_BIOCONTAINER % name)
    else:
        lint_ctx.warn(MESSAGE_WARN_NO_CONTAINER)


# TODO: Refactor following method into mulled util and then refactor same code out of mulled
# container resolver.
def mulled_container_name(namespace, targets):
    name = None

    if len(targets) == 1:
        target = targets[0]
        target_version = target.version
        tags = mulled_tags_for(namespace, target.package_name)

        if not tags:
            return None

        if target_version:
            for tag in tags:
                version, build = split_tag(tag)
                if version == target_version:
                    name = "%s:%s--%s" % (target.package_name, version, build)
                    break
        else:
            version, build = split_tag(tags[0])
            name = "%s:%s--%s" % (target.package_name, version, build)
    else:
        base_image_name = image_name(targets)
        tags = mulled_tags_for(namespace, base_image_name)
        if tags:
            name = "%s:%s" % (base_image_name, tags[0])

    if name:
        return "quay.io/%s/%s" % (namespace, name)
