import os
import string
import sys


SCRIPTS_DIRECTORY = os.path.dirname(__file__)
TEMPLATE_PATH = os.path.join(SCRIPTS_DIRECTORY, "slideshow_template.html")

with open(TEMPLATE_PATH, "r") as tfh:
    TEMPLATE = string.Template(tfh.read())


def main(argv=None):
    if argv is None:
        argv = sys.argv
    title = argv[1]
    markdown_source = argv[2]
    if len(argv) >= 4:
        output = argv[3]
    else:
        output = os.path.splitext(markdown_source)[0] + '.html'
    with open(markdown_source, "r") as s:
        content = s.read()
    html = TEMPLATE.safe_substitute(**{
        'title': title,
        'content': content,
    })
    print(html)
    with open(output, "w") as ofh:
        ofh.write(html)


if __name__ == "__main__":
    main()
