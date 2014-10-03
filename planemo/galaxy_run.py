
# planemo environment contains a copy a deeply restricted subset of Galaxy, so
# for Galaxy to function we need to deactivate any virtualenv we are in.
DEACTIVATE_COMMAND = "type deactivate >/dev/null 2>&1 && deactivate"

# Activate galaxy's virtualenv if present (needed for tests say but not for
# server because run.sh does this).
ACTIVATE_COMMAND = "[ -e .venv ] && . .venv/bin/activate"


DOWNLOAD_GALAXY = (
    "wget https://codeload.github.com/jmchilton/galaxy-central/tar.gz/master"
)
