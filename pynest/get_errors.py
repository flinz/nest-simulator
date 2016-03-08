import re

out = """\
# -*- coding: utf-8 -*-
#
# hl_api_connections.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

\"\"\"
Error classes
\"\"\"

import sys
import inspect


class NESTError(Exception):
    pass


# SLI Errors
class NESTSLIError(NESTError):
    nest_class = "SLIException"


# Kernel Errors
class NESTKernelError(NESTError):
    nest_class = "KernelException"


"""


m = 'class ([A-Za-z]+\w) : public ([A-Za-z]+)'


def process_exceptions(filename, superclass_str):
    with open(filename) as f:
        out = ""
        for l in f:
            names = re.match(m, l)
            if names:
                groups = names.groups()
                if groups[0] in ["KernelException", "SLIException"]:
                    continue
                out += """\
class NEST%sError(NESTKernelError):
    nest_class = "%s"


""" % (groups[0], groups[0])
        return out

out += process_exceptions('../nestkernel/exceptions.h', 'NESTKernelError')
out += process_exceptions('../sli/sliexceptions.h', 'NESTSLIError')
out += """\
def build_error_map():
    out = {None: NESTError}
    clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)

    for c in clsmembers:
        if hasattr(c[1], 'nest_class'):
            out[c[1].nest_class] = c[1]
    return out

ERROR_MAP = build_error_map()
"""

print out,