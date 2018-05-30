import collections, os, platform


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, str):
            yield from flatten(el)
        else:
            yield el


def open_with_os(path):
    sys = platform.system()

    if sys == "Darwin":
        os.system('open "%s"' % path)
    elif sys == "Windows":
        os.startfile(path)
    elif sys == "Linux":
        os.system('evince "%s"' % path)
    else:
        raise NotImplementedError("Unable to open files in this particular system: %s" % sys)
