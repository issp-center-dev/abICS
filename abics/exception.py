class Error(Exception):
    """Base class for exceptions in abICS"""
    pass


class InputError(Error):
    """Exception raised for errors in the input

    Attributes
    ----------
    message: str
        explanation
    """
    def __init__(self, message):
        self.message = message


class MatrixParseError(Error):
    """Exception raised for matrix parse errors

    Attributes
    ----------
    message: str
        explanation
    """
    def __init__(self, message):
        self.message = message
