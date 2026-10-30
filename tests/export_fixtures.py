"""Export precomputed extension fixtures through the real export functions."""
from aptaswitch_studio.extension_runtime import _finish_exports, _write_data_exports


def export_extension_result(result, *, log=()):
    return _finish_exports(result, _write_data_exports(result), log)
