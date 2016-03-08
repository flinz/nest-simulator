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

"""
Error classes
"""

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


class NESTUnknownModelNameError(NESTKernelError):
    nest_class = "UnknownModelName"


class NESTNewModelNameExistsError(NESTKernelError):
    nest_class = "NewModelNameExists"


class NESTUnknownModelIDError(NESTKernelError):
    nest_class = "UnknownModelID"


class NESTModelInUseError(NESTKernelError):
    nest_class = "ModelInUse"


class NESTUnknownSynapseTypeError(NESTKernelError):
    nest_class = "UnknownSynapseType"


class NESTUnknownNodeError(NESTKernelError):
    nest_class = "UnknownNode"


class NESTNoThreadSiblingsAvailableError(NESTKernelError):
    nest_class = "NoThreadSiblingsAvailable"


class NESTLocalNodeExpectedError(NESTKernelError):
    nest_class = "LocalNodeExpected"


class NESTNodeWithProxiesExpectedError(NESTKernelError):
    nest_class = "NodeWithProxiesExpected"


class NESTUnknownReceptorTypeError(NESTKernelError):
    nest_class = "UnknownReceptorType"


class NESTIncompatibleReceptorTypeError(NESTKernelError):
    nest_class = "IncompatibleReceptorType"


class NESTUnknownPortError(NESTKernelError):
    nest_class = "UnknownPort"


class NESTIllegalConnectionError(NESTKernelError):
    nest_class = "IllegalConnection"


class NESTUnknownThreadError(NESTKernelError):
    nest_class = "UnknownThread"


class NESTBadDelayError(NESTKernelError):
    nest_class = "BadDelay"


class NESTUnexpectedEventError(NESTKernelError):
    nest_class = "UnexpectedEvent"


class NESTUnsupportedEventError(NESTKernelError):
    nest_class = "UnsupportedEvent"


class NESTBadPropertyError(NESTKernelError):
    nest_class = "BadProperty"


class NESTBadParameterError(NESTKernelError):
    nest_class = "BadParameter"


class NESTDimensionMismatchError(NESTKernelError):
    nest_class = "DimensionMismatch"


class NESTDistributionErrorError(NESTKernelError):
    nest_class = "DistributionError"


class NESTSubnetExpectedError(NESTKernelError):
    nest_class = "SubnetExpected"


class NESTSimulationErrorError(NESTKernelError):
    nest_class = "SimulationError"


class NESTInvalidDefaultResolutionError(NESTKernelError):
    nest_class = "InvalidDefaultResolution"


class NESTInvalidTimeInModelError(NESTKernelError):
    nest_class = "InvalidTimeInModel"


class NESTStepMultipleRequiredError(NESTKernelError):
    nest_class = "StepMultipleRequired"


class NESTTimeMultipleRequiredError(NESTKernelError):
    nest_class = "TimeMultipleRequired"


class NESTGSLSolverFailureError(NESTKernelError):
    nest_class = "GSLSolverFailure"


class NESTNumericalInstabilityError(NESTKernelError):
    nest_class = "NumericalInstability"


class NESTMUSICPortUnconnectedError(NESTKernelError):
    nest_class = "MUSICPortUnconnected"


class NESTMUSICPortHasNoWidthError(NESTKernelError):
    nest_class = "MUSICPortHasNoWidth"


class NESTMUSICPortAlreadyPublishedError(NESTKernelError):
    nest_class = "MUSICPortAlreadyPublished"


class NESTMUSICSimulationHasRunError(NESTKernelError):
    nest_class = "MUSICSimulationHasRun"


class NESTMUSICChannelUnknownError(NESTKernelError):
    nest_class = "MUSICChannelUnknown"


class NESTMUSICPortUnknownError(NESTKernelError):
    nest_class = "MUSICPortUnknown"


class NESTMUSICChannelAlreadyMappedError(NESTKernelError):
    nest_class = "MUSICChannelAlreadyMapped"


class NESTInterpreterErrorError(NESTKernelError):
    nest_class = "InterpreterError"


class NESTWrappedThreadExceptionError(NESTKernelError):
    nest_class = "WrappedThreadException"


class NESTDivisionByZeroError(NESTKernelError):
    nest_class = "DivisionByZero"


class NESTTypeMismatchError(NESTKernelError):
    nest_class = "TypeMismatch"


class NESTSystemSignalError(NESTKernelError):
    nest_class = "SystemSignal"


class NESTRangeCheckError(NESTKernelError):
    nest_class = "RangeCheck"


class NESTArgumentTypeError(NESTKernelError):
    nest_class = "ArgumentType"


class NESTBadParameterValueError(NESTKernelError):
    nest_class = "BadParameterValue"


class NESTDictErrorError(NESTKernelError):
    nest_class = "DictError"


class NESTUndefinedNameError(NESTKernelError):
    nest_class = "UndefinedName"


class NESTEntryTypeMismatchError(NESTKernelError):
    nest_class = "EntryTypeMismatch"


class NESTStackUnderflowError(NESTKernelError):
    nest_class = "StackUnderflow"


class NESTIOErrorError(NESTKernelError):
    nest_class = "IOError"


class NESTUnaccessedDictionaryEntryError(NESTKernelError):
    nest_class = "UnaccessedDictionaryEntry"


class NESTDynamicModuleManagementErrorError(NESTKernelError):
    nest_class = "DynamicModuleManagementError"


class NESTNamingConflictError(NESTKernelError):
    nest_class = "NamingConflict"


class NESTNotImplementedError(NESTKernelError):
    nest_class = "NotImplemented"


def build_error_map():
    out = {None: NESTError}
    clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)

    for c in clsmembers:
        if hasattr(c[1], 'nest_class'):
            out[c[1].nest_class] = c[1]
    return out

ERROR_MAP = build_error_map()
