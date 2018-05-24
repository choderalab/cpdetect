from pkg_resources import resource_filename
import os


def get_fn(filename, written=False):
    """Get the full path to one of the reference files shipped for testing

        These files are in torsionfit/testing/reference

    :param
        name: str
            Name of file to load

    :returns
        fn : str
            full path to file
    """
    if written:
        fn = resource_filename('cpdetect', os.path.join('tests', 'files', 'writes', filename))
    else:
        fn = resource_filename('cpdetect', os.path.join('tests', 'files', filename))

    #if not os.path.exists(fn):
    #    raise ValueError('%s does not exist. If you just added it you will have to re install' % fn)

    return fn
