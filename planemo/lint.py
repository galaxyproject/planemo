import os

from planemo.xml import validation


def lint_xsd(lint_ctx, schema_path, path):
    name = os.path.basename(path)
    validator = validation.get_validator(require=True)
    validation_result = validator.validate(schema_path, path)
    if not validation_result.passed:
        msg = "Invalid %s found. Errors [%s]"
        msg = msg % (name, validation_result.output)
        lint_ctx.error(msg)
    else:
        lint_ctx.info("%s found and appears to be valid XML" % name)
