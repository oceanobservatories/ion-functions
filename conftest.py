"""
Stub out the legacy `nose` test-tagging import so pytest can collect test
modules under this Python version, where the installed `nose` package is
broken (it imports the `imp` module, removed in Python 3.12). This does not
touch nose itself or any test file; it only provides a no-op `attr` decorator
for `from nose.plugins.attrib import attr`.
"""
import sys
import types


def attr(*args, **kwargs):
    def decorator(func):
        return func
    return decorator


_nose = types.ModuleType('nose')
_nose_plugins = types.ModuleType('nose.plugins')
_nose_plugins_attrib = types.ModuleType('nose.plugins.attrib')
_nose_plugins_attrib.attr = attr
_nose.plugins = _nose_plugins
_nose_plugins.attrib = _nose_plugins_attrib

sys.modules['nose'] = _nose
sys.modules['nose.plugins'] = _nose_plugins
sys.modules['nose.plugins.attrib'] = _nose_plugins_attrib
