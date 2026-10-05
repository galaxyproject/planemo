import importlib.metadata

try:
    planemo_metadata = importlib.metadata.metadata("planemo")
except importlib.metadata.PackageNotFoundError:
    planemo_metadata = importlib.metadata.metadata("planemo-cli")

__version__ = "0.75.48.dev0"

PROJECT_NAME = "planemo"
PROJECT_EMAIL = planemo_metadata["Author-email"].split(" ")[-1]
PROJECT_AUTHOR = PROJECT_USERNAME = "galaxyproject"

PROJECT_URL = "https://github.com/galaxyproject/planemo"
RAW_CONTENT_URL = f"https://raw.github.com/{PROJECT_USERNAME}/{PROJECT_NAME}/master/"
