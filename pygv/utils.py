"""Utility functions for PyGV."""

import os


def check_accessibility(
    file_path: str, allow_remote: bool = False, raise_except: bool = True
):
    """Check if a file is accessible.

    Parameters
    ----------
    file_path : str
        File path
    allow_remote : bool
        Allow remote file access via http/https/ftp protocols.
    raise_except : bool
        If set as True and the file is not accessible,
        this function raises a ValueError

    Returns
    -------
    accessibility : bool
        Returns True for an accessible file
    """
    accessibility = False
    if os.path.exists(file_path):
        accessibility = True
    elif allow_remote:
        if file_path.startswith("http") or file_path.startswith("ftp"):
            accessibility = True

    if raise_except and not accessibility:
        raise ValueError(f"File {file_path} is not accessible.")
    return accessibility
