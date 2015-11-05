

def run_galaxy(ctx, path, job_path, **kwds):
    kwds["cwl"] = True
    conformance_test = kwds.get("conformance_test", False)
    with conditionally_captured_io(conformance_test):
        with galaxy_serve.serve_daemon(ctx, [path], **kwds) as config:
            try:
                cwl_run = cwl.run_cwl_tool(path, job_path, config, **kwds)
            except Exception:
                io.warn("Problem running cwl tool...")
                print(config.log_contents)
                raise

    print(cwl_run.cwl_command_state)
    return 0
