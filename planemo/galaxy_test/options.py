import os

OUTPUT_DFEAULTS = {
    "output": "tool_test_output.html",
    "output_json": "tool_test_output.json",
    "output_xunit": None,
}


def process_defaults(ctx, kwds):
    for name, default in OUTPUT_DFEAULTS.items():
        _populate_default_output(ctx, name, kwds, default)


def _populate_default_output(ctx, type, kwds, default):
    kwd_key = "test_%s" % type
    kwd_value = kwds.get(kwd_key, None)
    if kwd_value is None:
        global_config = ctx.global_config
        global_config_key = "default_test_%s" % type
        if global_config_key in global_config:
            default_value = global_config[global_config_key]
        else:
            default_value = default

        if default_value:
            default_value = os.path.abspath(default_value)
        kwds[kwd_key] = default_value
