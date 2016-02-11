"""
    genonets_exceptions
    ~~~~~~~~~~~~~~~~~~~

    Defines custom exceptions.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""


class GenonetsError(Exception):
    def __init__(self, errId, info=""):
        self.errId = errId
        self.info = info

    def __str__(self):
        return repr(self.errId)
