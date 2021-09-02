"""All of the exceptions which are raised in the project package"""


class InputError(Exception):
    """Raises an error if the user does not put in the correct input files"""
    def __init__(self, message="Input is unacceptable."):
        self.message = message
        super().__init__(message)


class QueryNotFoundError(Exception):
    """Raises an error if the user's query cannot be found"""
    def __init__(self, message="Query is unacceptable."):
        self.message = message
        super().__init__(message)


class FileMissingError(Exception):
    """Raises an error if the user asks for a file which is not there"""
    def __init__(self, message="This file is not found. Try running parse() once again."):
        self.message = message
        super().__init__(message)
