import functools
import inspect
import warnings

__version__ = '2.1.0'


class deprecated(object):
    """
    Decorator which, when used, will result in automatic reporting to user of the 
    deprecation and indicate replacement object.
    """

    def __init__(self, reason):
        if inspect.isclass(reason) or inspect.isfunction(reason):
            raise TypeError('Reason or deprecation must be supplied.')
        self.reason = reason

    def __call__(self, class_or_function):
        if inspect.isfunction(class_or_function):
            if hasattr(class_or_function, 'func_code'):
                _code = class_or_function.func_code
            else:
                _code = class_or_function.__code__
            fmt = 'Call to deprecated function {name} ({reason}).'
            filename = _code.co_filename
            line_number = _code.co_firstlineno + 1

        elif inspect.isclass(class_or_function):
            fmt = 'Call to deprecated class {name} ({reason}).'
            filename = class_or_function.__module__
            line_number = 1
        else:
            raise TypeError(type(class_or_function))

        message = fmt.format(name=class_or_function.__name__, reason=self.reason)

        @functools.wraps(class_or_function)
        def new_func(*args, **kwargs):
            warnings.simplefilter('always', DeprecationWarning)
            warnings.warn_explicit(message, category=DeprecationWarning, filename=filename, lineno=line_number)
            warnings.simplefilter('default', DeprecationWarning)
            return class_or_function

        return new_func
