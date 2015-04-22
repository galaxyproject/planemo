"""
This file is a mess, it is a merge of random stuff that is in galaxy.util and
stuff that was in pulsar.util. This should be reworked to only contain stuff in
galaxy.util and the rest should be moved into galaxy.util.pulsar_io or something
like that.
"""
import os
import stat
try:
    import grp
except ImportError:
    grp = None
try:
    import docutils.core
    import docutils.writers.html4css1
except ImportError:
    pass
import errno
from six.moves.urllib import parse
from six import text_type
from tempfile import NamedTemporaryFile
import logging
log = logging.getLogger(__name__)

BUFFER_SIZE = 4096
DEFAULT_ENCODING = os.environ.get('GALAXY_DEFAULT_ENCODING', 'utf-8')


def enum(**enums):
    """
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    return type('Enum', (), enums)


def copy_to_path(object, path):
    """
    Copy file-like object to path.
    """
    output = open(path, 'wb')
    _copy_and_close(object, output)


def _copy_and_close(object, output):
    try:
        while True:
            buffer = object.read(BUFFER_SIZE)
            if not buffer:
                break
            output.write(buffer)
    finally:
        output.close()


def copy_to_temp(object):
    """
    Copy file-like object to temp file and return
    path.
    """
    temp_file = NamedTemporaryFile(delete=False)
    _copy_and_close(object, temp_file)
    return temp_file.name


def is_in_directory(file, directory, local_path_module=os.path):
    """
    Return true, if the common prefix of both is equal to directory
    e.g. /a/b/c/d.rst and directory is /a/b, the common prefix is /a/b

    Heavily inspired by similar method in from Galaxy's BaseJobRunner class.
    """

    # Make both absolute.
    directory = local_path_module.abspath(directory)
    file = local_path_module.abspath(file)
    return local_path_module.commonprefix([file, directory]) == directory


in_directory = is_in_directory  # For compat. w/Galaxy.


def umask_fix_perms(path, umask, unmasked_perms, gid=None):
    """
    umask-friendly permissions fixing
    """
    perms = unmasked_perms & ~umask
    try:
        st = os.stat(path)
    except OSError as e:
        log.exception('Unable to set permissions or group on %s' % path)
        return
    # fix modes
    if stat.S_IMODE(st.st_mode) != perms:
        try:
            os.chmod(path, perms)
        except Exception as e:
            log.warning('Unable to honor umask (%s) for %s, tried to set: %s but mode remains %s, error was: %s' % (oct(umask),
                                                                                                                    path,
                                                                                                                    oct(perms),
                                                                                                                    oct(stat.S_IMODE(st.st_mode)),
                                                                                                                    e))
    # fix group
    if gid is not None and st.st_gid != gid:
        try:
            os.chown(path, -1, gid)
        except Exception as e:
            try:
                desired_group = grp.getgrgid(gid)
                current_group = grp.getgrgid(st.st_gid)
            except:
                desired_group = gid
                current_group = st.st_gid
            log.warning('Unable to honor primary group (%s) for %s, group remains %s, error was: %s' % (desired_group,
                                                                                                        path,
                                                                                                        current_group,
                                                                                                        e))


def xml_text(root, name=None):
    """Returns the text inside an element"""
    if name is not None:
        # Try attribute first
        val = root.get(name)
        if val:
            return val
        # Then try as element
        elem = root.find(name)
    else:
        elem = root
    if elem is not None and elem.text:
        text = ''.join(elem.text.splitlines())
        return text.strip()
    # No luck, return empty string
    return ''


# asbool implementation pulled from PasteDeploy
truthy = frozenset(['true', 'yes', 'on', 'y', 't', '1'])
falsy = frozenset(['false', 'no', 'off', 'n', 'f', '0'])


def asbool(obj):
    if isinstance(obj, basestring):
        obj = obj.strip().lower()
        if obj in truthy:
            return True
        elif obj in falsy:
            return False
        else:
            raise ValueError("String is not true/false: %r" % obj)
    return bool(obj)


string_as_bool = asbool


def force_symlink(source, link_name):
    try:
        os.symlink(source, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(source, link_name)
        else:
            raise e


def listify( item, do_strip=False ):
    """
    Make a single item a single item list, or return a list if passed a
    list.  Passing a None returns an empty list.
    """
    if not item:
        return []
    elif isinstance( item, list ):
        return item
    elif isinstance( item, basestring ) and item.count( ',' ):
        if do_strip:
            return [token.strip() for token in item.split( ',' )]
        else:
            return item.split( ',' )
    else:
        return [ item ]


def mask_password_from_url( url ):
    """
    Masks out passwords from connection urls like the database connection in galaxy.ini

    >>> mask_password_from_url( 'sqlite+postgresql://user:password@localhost/' )
    'sqlite+postgresql://user:********@localhost/'
    >>> mask_password_from_url( 'amqp://user:amqp@localhost' )
    'amqp://user:********@localhost'
    >>> mask_password_from_url( 'amqp://localhost')
    'amqp://localhost'
    """
    split = parse.urlsplit(url)
    if split.password:
        if url.count(split.password) == 1:
            url = url.replace(split.password, "********")
        else:
            # This can manipulate the input other than just masking password,
            # so the previous string replace method is preferred when the
            # password doesn't appear twice in the url
            split = split._replace(netloc=split.netloc.replace("%s:%s" % (split.username, split.password), '%s:********' % split.username))
            url = parse.urlunsplit(split)
    return url


def unicodify( value, encoding=DEFAULT_ENCODING, error='replace', default=None ):
    """
    Returns a unicode string or None
    """

    if isinstance( value, text_type ):
        return value
    try:
        return unicode( str( value ), encoding, error )
    except:
        return default


def rst_to_html( s ):
    """Convert a blob of reStructuredText to HTML"""
    log = logging.getLogger( "docutils" )

    class FakeStream( object ):
        def write( self, str ):
            if len( str ) > 0 and not str.isspace():
                log.warn( str )
    return unicodify( docutils.core.publish_string( s,
                      writer=docutils.writers.html4css1.Writer(),
                      settings_overrides={ "embed_stylesheet": False, "template": os.path.join(os.path.dirname(__file__), "docutils_template.txt"), "warning_stream": FakeStream() } ) )
