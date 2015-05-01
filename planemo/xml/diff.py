
def diff(x1, x2, reporter=None):
    return 0 if xml_compare(x1, x2, reporter) else 1


# From
# bitbucket.org/ianb/formencode/src/tip/formencode/doctest_xml_compare.py
# with (PSF license)
def xml_compare(x1, x2, reporter=None):
    if reporter is None:
        def reporter(x):
            return None

    if x1.tag != x2.tag:
        reporter('Tags do not match: %s and %s' % (x1.tag, x2.tag))
        return False
    for name, value in x1.attrib.items():
        if x2.attrib.get(name) != value:
            reporter('Attributes do not match: %s=%r, %s=%r'
                     % (name, value, name, x2.attrib.get(name)))
            return False
    for name in x2.attrib.keys():
        if name not in x1.attrib:
            reporter('x2 has an attribute x1 is missing: %s'
                     % name)
            return False
    if not text_compare(x1.text, x2.text):
        reporter('text: %r != %r' % (x1.text, x2.text))
        return False
    if not text_compare(x1.tail, x2.tail):
        reporter('tail: %r != %r' % (x1.tail, x2.tail))
        return False
    return _compare_children(x1, x2, reporter)


def _compare_children(x1, x2, reporter):
    cl1 = x1.getchildren()
    cl2 = x2.getchildren()
    if len(cl1) != len(cl2):
        reporter('children length differs, %i != %i'
                 % (len(cl1), len(cl2)))
        return False
    i = 0
    for c1, c2 in zip(cl1, cl2):
        i += 1
        if not xml_compare(c1, c2, reporter=reporter):
            reporter('children %i do not match: %s'
                     % (i, c1.tag))
            return False
    return True


def text_compare(t1, t2):
    if not t1 and not t2:
        return True
    return (t1 or '').strip() == (t2 or '').strip()
