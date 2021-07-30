"""
    A collection of custom error classes.

"""


class BaseGenonetsError(Exception):
    """ Base class for all error classes in the module. """


class FileError(BaseGenonetsError):
    """ Raised in case of a file/directory related error. """


class MissingMandatoryArgumentError(BaseGenonetsError):
    """ Raised if a mandatory command-line argument is missing. """


class InternalError(BaseGenonetsError):
    """ Raised if the error is not directly related to user input. """


class GenonetsError(Exception):
    def __init__(self, errId, info=""):
        self.errId = errId
        self.info = info

    def __str__(self):
        return repr(self.errId)
