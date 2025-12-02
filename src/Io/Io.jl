"""

Module for the reading and writing data using JLD2.
"""
module Io


using JLD2
using Printf
using EllipsisNotation


include("FromFile.jl")

export inputs_from_file, evals_from_file, efunc_from_file, surfaces_from_file


include("ToFile.jl")

export inputs_to_file, evals_to_file, efuncs_to_file, surfaces_to_file


end
