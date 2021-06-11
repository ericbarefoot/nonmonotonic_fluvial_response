def _get_version():
    """
    Extract version number from single file, and make it availabe everywhere.
    """
    from . import _version
    return _version.__version__()
